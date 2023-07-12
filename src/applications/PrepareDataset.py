import os, argparse, sys, csv, platform, subprocess, shutil, posixpath, yaml
from pathlib import Path
from datetime import date
import pandas as pd
import SimpleITK as sitk
from tqdm import tqdm
import numpy as np
from skimage.measure import label
from copy import deepcopy

from FigureGenerator.screenshot_maker import figure_generator

# check against all these modality ID strings with extensions
modality_id_dict = {
    "T1": ["t1", "t1pre", "t1precontrast"],
    "T1GD": ["t1ce", "t1gd", "t1post", "t1postcontrast", "t1gallodinium", "t1c"],
    "T2": ["t2"],
    "FLAIR": ["flair", "fl", "t2flair"],
}


def setup_parser():
    copyrightMessage = (
        "Contact: admin@fets.ai\n\n"
        + "This program is NOT FDA/CE approved and NOT intended for clinical use.\nCopyright (c) "
        + str(date.today().year)
        + " University of Pennsylvania. All rights reserved."
    )
    parser = argparse.ArgumentParser(
        prog="PrepareDataset",
        formatter_class=argparse.RawTextHelpFormatter,
        description="This application calls the BraTSPipeline for all input images and stores the final and intermediate files separately.\n\n"
        + copyrightMessage,
    )
    parser.add_argument(
        "-inputCSV",
        type=str,
        help="The absolute path of the input CSV file containing the list of subjects and their corresponding images",
        required=True,
    )
    parser.add_argument(
        "-outputDir",
        type=str,
        help="The output dir to write the results",
        required=True,
    )

    return parser


def _get_relevant_dicom_tags(filename: str) -> dict:
    """
    This function reads the relevant DICOM tags from the input DICOM directory.

    Args:
        filename (str): The input DICOM filename.

    Returns:
        dict: The relevant DICOM tags.
    """
    input_dicom_dir = filename
    if os.path.isfile(filename):
        input_dicom_dir = os.path.dirname(filename)

    output_dict = {}
    try:
        series_IDs = sitk.ImageSeriesReader.GetGDCMSeriesIDs(input_dicom_dir)
        # if len(series_IDs) > 1:
        #     print(
        #         f"WARNING: Multiple series IDs detected in {input_dicom_dir}.",
        #         file=sys.stderr,
        #     )

        series_file_names = sitk.ImageSeriesReader.GetGDCMSeriesFileNames(
            input_dicom_dir, series_IDs[0]
        )
        series_reader = sitk.ImageSeriesReader()
        series_reader.SetFileNames(series_file_names)
        series_reader.MetaDataDictionaryArrayUpdateOn()
        series_reader.LoadPrivateTagsOn()
        itk_image = series_reader.Execute()
        output_dict = {
            "Resolution": str(itk_image.GetSpacing()).replace(" ", ""),
        }
        # although _technically_ the metadata is different for each slice, we'll just use the first slice's metadata, since the rest is not relevant for our purposes
        ## reference: https://simpleitk.readthedocs.io/en/master/link_DicomSeriesReadModifyWrite_docs.html
        keys_to_extract = [
            "0008|0070",  # Manufacturer
            "0008|1090",  # Manufacturer's Model Name
            "0008|103e",  # Series Description
            "0008|0021",  # Series Date
            "0008|0031",  # Series Time
        ]
        keys_to_extract = {
            "0008|0070": "Manufacturer",
            "0008|1090": "Manufacturer's Model Name",
            "0008|0022": "Acquisition Date",
            "0008|0032": "Acquisition Time",
            "0018|0087": "Magnetic Field Strength",
            "0018|1050": "Slice Thickness",
            "0018|0088": "Spacing Between Slices",
            "0010|1010": "Patient's Age",
            "0010|0040": "Patient's Sex",
        }
        for key in keys_to_extract:
            output_dict[keys_to_extract[key]] = series_reader.GetMetaData(0, key)
    except RuntimeError as e:
        # print(
        #     f"WARNING: Could not read DICOM tags from {input_dicom_dir}.",
        # )
        pass

    return output_dict


def _read_image_with_min_check(filename):
    """
    This function fixes negatives by scaling the image according to the following logic:
    if min(input) < 0:
    for all x in image:
        if x != 0:
        x -= min

    Args:
        filename (str): The input filename.

    Returns:
        sitk.Image: The read image.
        int: The negative count.
    """
    input_image = sitk.ReadImage(filename)
    input_image_array = sitk.GetArrayFromImage(input_image)
    min = np.min(input_image_array)

    # the threshold above which an error is displayed, otherwise, the intensities are scaled
    max_negative_count_threshold = 5000

    if min < 0:
        blobs = input_image_array < 0
        all_labels_nonZero = np.nonzero(label(blobs))
        _, counts = np.unique(all_labels_nonZero, return_counts=True)

        if np.max(counts) < max_negative_count_threshold:
            output_array = deepcopy(input_image_array)
            mask = output_array != 0
            output_array[mask] = output_array[mask] - min
            output_image = sitk.GetImageFromArray(output_array)
            output_image.CopyInformation(input_image)
            sitk.WriteImage(output_image, filename)
            return 0
        else:
            return counts.astype(int)

    return 0


def parse_csv_header(filename):
    """
    Read filename and return the parsed headers.

    Args:
        filename (str): The input filename.

    Returns:
        dict: The parsed headers.
    """
    with open(filename, "r") as csvfile:
        datareader = csv.reader(csvfile)

        headers = {}  # save headers
        for row in datareader:
            for col in row:
                temp = col.lower()  # convert to lower case
                temp = temp.replace(" ", "")  # remove spaces
                temp = temp.replace("_", "")  # remove underscores
                temp = temp.replace("-", "")  # remove dashes
                if (
                    (temp == "patientid")
                    or (temp == "subjectid")
                    or (temp == "subject")
                    or (temp == "subid")
                ):
                    headers["ID"] = col
                elif (
                    (temp == "timepoint")
                    or (temp == "tp")
                    or (temp == "time")
                    or (temp == "series")
                    or (temp == "subseries")
                ):
                    headers["Timepoint"] = col
                else:
                    for key in modality_id_dict.keys():
                        if temp in modality_id_dict[key]:
                            headers[key] = col
                            break

    if "Timepoint" not in headers:
        headers["Timepoint"] = None
    return headers


def copyFilesToCorrectLocation(interimOutputDir, finalSubjectOutputDir, subjectID):
    """
    This function copies the intermediate files and final outputs to correct location and if these are absent, returns a bool flag stating that brats pipeline needs to run again
    """

    # copy files to correct location for inference and training
    runBratsPipeline = False
    input_files = {
        "T1": os.path.join(interimOutputDir, "T1_to_SRI.nii.gz"),
        "T1GD": os.path.join(interimOutputDir, "T1CE_to_SRI.nii.gz"),
        "T2": os.path.join(interimOutputDir, "T2_to_SRI.nii.gz"),
        "FLAIR": os.path.join(interimOutputDir, "FL_to_SRI.nii.gz"),
    }
    expected_outputs = {
        "ID": subjectID,
        "T1": posixpath.join(finalSubjectOutputDir, subjectID + "_t1.nii.gz"),
        "T1GD": posixpath.join(finalSubjectOutputDir, subjectID + "_t1ce.nii.gz"),
        "T2": posixpath.join(finalSubjectOutputDir, subjectID + "_t2.nii.gz"),
        "FLAIR": posixpath.join(finalSubjectOutputDir, subjectID + "_flair.nii.gz"),
    }

    for key in input_files.keys():
        if not os.path.exists(expected_outputs[key]):
            if os.path.exists(input_files[key]):
                shutil.copyfile(input_files[key], expected_outputs[key])
            else:
                runBratsPipeline = True

    return runBratsPipeline, expected_outputs


class Preparator:
    def __init__(self, input_csv: str, output_dir: str):
        self.input_csv = input_csv
        self.input_dir = str(Path(input_csv).parent)
        self.output_dir = os.path.normpath(output_dir)
        self.interim_output_dir = os.path.join(self.output_dir, "DataForQC")
        self.final_output_dir = os.path.join(self.output_dir, "DataForFeTS")
        self.subjects_file = os.path.join(self.final_output_dir, "processed_data.csv")
        self.neg_subjects_file = os.path.join(
            self.final_output_dir, "QC_subjects_with_negative_intensities.csv"
        )
        self.failing_subjects_file = os.path.join(
            self.final_output_dir, "QC_subjects_with_bratspipeline_error.csv"
        )
        self.__init_out_dfs()
        self.stdout_log = os.path.join(self.output_dir, "preparedataset_stdout.txt")
        self.stderr_log = os.path.join(self.output_dir, "preparedataset_stderr.txt")
        self.dicom_tag_information_to_write_collab = {}
        self.dicom_tag_information_to_write_anon = {}
        self.brats_pipeline_exe = "BraTSPipeline"

        if platform.system() == "Windows":
            self.brats_pipeline_exe += ".exe"

    def __init_out_dfs(self):
        self.subjects = pd.DataFrame(
            columns=["SubjectID", "Timepoint", "T1", "T1GD", "T2", "FLAIR"]
        )

        self.neg_subjects = pd.DataFrame(
            columns=["SubjectID", "Timepoint", "Modality", "Count"]
        )
        self.failing_subjects = pd.DataFrame(columns=["SubjectID", "Timepoint"])

    def validate(self):
        assert os.path.exists(self.input_csv), "Input CSV file not found"

        assert shutil.which(
            self.brats_pipeline_exe
        ) is not None, "BraTS Pipeline executable not found, please contact admin@fets.ai for help."

    def process_data(self):
        items = self.subjects_df.iterrows()
        total = self.subjects_df.shape[0]
        desc = "Preparing Dataset (1-10 min per subject)"
        for idx, row in tqdm(items, total=total, desc=desc):
            self.process_row(idx, row)

    def process_row(self, idx: int, row: pd.Series):
        parsed_headers = self.parsed_headers
        bratsPipeline_exe = self.brats_pipeline_exe

        subject_id = row[parsed_headers["ID"]]
        subject_id_timepoint = subject_id

        # create QC and Final output dirs for each subject
        interimOutputDir_actual = posixpath.join(
            self.interim_output_dir, subject_id_timepoint
        )
        finalSubjectOutputDir_actual = posixpath.join(
            self.final_output_dir, subject_id_timepoint
        )

        # per the data ingestion step, we are creating a new folder called timepoint, can join timepoint to subjectid if needed
        if parsed_headers["Timepoint"] is not None:
            timepoint = row[parsed_headers["Timepoint"]]
            subject_id_timepoint += "_" + timepoint
            interimOutputDir_actual = posixpath.join(interimOutputDir_actual, timepoint)
            finalSubjectOutputDir_actual = posixpath.join(
                finalSubjectOutputDir_actual, timepoint
            )

        # get the relevant dicom tags
        self.dicom_tag_information_to_write_collab[subject_id_timepoint] = {}
        self.dicom_tag_information_to_write_anon[str(idx)] = {}
        for modality in ["T1", "T1GD", "T2", "FLAIR"]:
            tags_from_modality = _get_relevant_dicom_tags(row[parsed_headers[modality]])
            self.dicom_tag_information_to_write_collab[subject_id_timepoint][
                modality
            ] = tags_from_modality
            self.dicom_tag_information_to_write_anon[str(idx)][modality] = tags_from_modality

        Path(interimOutputDir_actual).mkdir(parents=True, exist_ok=True)
        Path(finalSubjectOutputDir_actual).mkdir(parents=True, exist_ok=True)
        # if files already exist in DataForQC, then copy to DataForFeTS, and if files exist in DataForFeTS, then skip
        runBratsPipeline, _ = copyFilesToCorrectLocation(
            interimOutputDir_actual, finalSubjectOutputDir_actual, subject_id_timepoint
        )

        # check if the files exist already, if so, skip
        if not runBratsPipeline:
            return

        command = (
            bratsPipeline_exe
            + " -t1 "
            + row[parsed_headers["T1"]]
            + " -t1c "
            + row[parsed_headers["T1GD"]]
            + " -t2 "
            + row[parsed_headers["T2"]]
            + " -fl "
            + row[parsed_headers["FLAIR"]]
            + " -o "
            + interimOutputDir_actual
        )

        with open(self.stdout_log, "a+") as out, open(self.stderr_log, "a+") as err:
            out.write(f"***\n{command}\n***")
            err.write(f"***\n{command}\n***")
            subprocess.Popen(command, stdout=out, stderr=err, shell=True).wait()

        runBratsPipeline, outputs = copyFilesToCorrectLocation(
            interimOutputDir_actual,
            finalSubjectOutputDir_actual,
            subject_id_timepoint,
        )

        if runBratsPipeline:
            # The BraTS command failed, and no files were found
            # flag this subject as failing
            failing_data = {"SubjectID": subject_id, "Timepoint": timepoint}
            failing_subject = pd.DataFrame(failing_data, index=[0])
            self.failing_subjects = pd.concat([self.failing_subjects, failing_subject])
            return

        # store the outputs in a dictionary when there are no errors
        negatives_detected = False
        for modality in ["T1", "T1GD", "T2", "FLAIR"]:
            count = _read_image_with_min_check(outputs[modality])
            # if there are any negative values, then store the subjectid, timepoint, modality and count of negative values
            if count == 0:
                continue
            neg_data = {
                "SubjectID": subject_id,
                "Timepoint": timepoint,
                "Modality": modality,
                "Count": count,
            }
            neg_subject = pd.DataFrame(neg_data, index=[0])
            self.neg_subjects = pd.concat([self.neg_subjects, neg_subject])
            negatives_detected = True

        # store the outputs in a dictionary when there are no errors
        if negatives_detected:
            return

        subject_data = {
            "SubjectID": subject_id,
            "Timepoint": timepoint,
            "T1": outputs["T1"],
            "T1GD": outputs["T1GD"],
            "T2": outputs["T2"],
            "FLAIR": outputs["FLAIR"],
        }
        subject = pd.DataFrame(subject_data, index=[0])
        self.subjects = pd.concat([self.subjects,subject,])

        # save the screenshot
        images = (",").join(
            [
                outputs["T1"],
                outputs["T1GD"],
                outputs["T2"],
                outputs["FLAIR"],
            ]
        )
        ylabels = (",").join(
            [
                "T1",
                "T1GD",
                "T2",
                "FLAIR",
            ]
        )
        figure_generator(
            images,
            ylabels,
            os.path.join(interimOutputDir_actual, "screenshot.png"),
            flip_sagittal=True,
            flip_coronal=True,
        )

    def write(self):
        if self.subjects.shape[0]:
            self.subjects.to_csv(self.subjects_file, index=False)
        if self.neg_subjects.shape[0]:
            self.neg_subjects.to_csv(self.neg_subjects_file, index=False)
        if self.failing_subjects.shape[0]:
            self.failing_subjects.to_csv(self.failing_subjects_file, index=False)

    def read(self):
        self.parsed_headers = parse_csv_header(self.input_csv)
        self.subjects_df = pd.read_csv(self.input_csv)
        if os.path.exists(self.subjects_file):
            self.subjects = pd.read_csv(self.subjects_file)
        if os.path.exists(self.neg_subjects_file):
            self.neg_subjects = pd.read_csv(self.neg_subjects_file)
        if os.path.exists(self.failing_subjects_file):
            self.failing_subjects = pd.read_csv(self.failing_subjects_file)


def main():
    parser = setup_parser()
    args = parser.parse_args()

    prep = Preparator(args.inputCSV, args.outputDir)
    prep.validate()
    prep.read()
    prep.process_data()
    prep.write()


if __name__ == "__main__":
    if platform.system() == "Darwin":
        sys.exit("macOS is not supported")
    else:
        main()
