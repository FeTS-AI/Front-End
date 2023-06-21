import os, argparse, sys, csv, platform, subprocess, shutil
from pathlib import Path
from datetime import date
import pandas as pd
import SimpleITK as sitk
from tqdm import tqdm
import numpy as np
from skimage.measure import label
from copy import deepcopy

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
    Read filename and return a list of dictionaries that have the csv contents
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
            #     elif (
            #         (temp == "t1gd")
            #         or (temp == "t1ce")
            #         or (temp == "t1post")
            #         or (temp == "t1postcontrast")
            #         or (temp == "t1gallodinium")
            #         or (temp == "t1c")
            #     ):
            #         headers["T1GD"] = col
            #     elif (
            #         (temp == "t1")
            #         or (temp == "t1pre")
            #         or (temp == "t1precontrast")
            #         or (temp == "t1p")
            #     ):
            #         headers["T1"] = col
            #     elif temp == "t2":
            #         headers["T2"] = col
            #     elif (
            #         (temp == "t2flair")
            #         or (temp == "flair")
            #         or (temp == "fl")
            #         or ("fl" in temp)
            #         or ("t2fl" in temp)
            #     ):
            #         headers["FLAIR"] = col
            # break

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
        "T1": os.path.join(finalSubjectOutputDir, subjectID + "_t1.nii.gz"),
        "T1GD": os.path.join(finalSubjectOutputDir, subjectID + "_t1ce.nii.gz"),
        "T2": os.path.join(finalSubjectOutputDir, subjectID + "_t2.nii.gz"),
        "FLAIR": os.path.join(finalSubjectOutputDir, subjectID + "_flair.nii.gz"),
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
        script_path = os.path.dirname(os.path.abspath(__file__))

        self.input_csv = input_csv
        self.output_dir = os.path.normpath(output_dir)
        self.interim_output_dir = os.path.join(self.output_dir, "DataForQC")
        self.final_output_dir = os.path.join(self.output_dir, "DataForFeTS")
        self.parsed_headers = parse_csv_header(input_csv)
        self.subjects_df = pd.read_csv(input_csv)
        self.__init_out_dfs()
        self.stdout_log = os.path.join(self.output_dir, "preparedataset_stdout.txt")
        self.stderr_log = os.path.join(self.output_dir, "preparedataset_stderr.txt")
        self.brats_pipeline_exe = os.path.join(script_path, "BraTSPipeline")

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

        assert os.path.exists(
            self.brats_pipeline_exe
        ), "BraTS Pipeline executable not found, please contact admin@fets.ai for help."

    def process_data(self):
        total = self.subjects_df.shape[0]
        for _, row in tqdm(self.subjects_df.iterrows(), total=total):
            self.process_row(row)

    def process_row(self, row: pd.Series):
        parsed_headers = self.parsed_headers
        bratsPipeline_exe = self.brats_pipeline_exe

        subject_id = row[parsed_headers["ID"]]
        subject_id_timepoint = subject_id

        # create QC and Final output dirs for each subject
        interimOutputDir_actual = os.path.join(self.interim_output_dir, subject_id_timepoint)
        finalSubjectOutputDir_actual = os.path.join(self.final_output_dir, subject_id_timepoint)

        # per the data ingestion step, we are creating a new folder called timepoint, can join timepoint to subjectid if needed
        string_to_write_to_logs = f"Processing {subject_id}"
        if parsed_headers["Timepoint"] is not None:
            timepoint = row[parsed_headers["Timepoint"]]
            subject_id_timepoint += "_" + timepoint
            interimOutputDir_actual = os.path.join(interimOutputDir_actual, timepoint)
            finalSubjectOutputDir_actual = os.path.join(
                finalSubjectOutputDir_actual, timepoint
            )
            # print(f"Processing {subject_id} timepoint {timepoint}")
            string_to_write_to_logs = f"Processing {subject_id} timepoint {timepoint}"

        with open(self.stdout_log, "a+") as out:
            out.write(string_to_write_to_logs)

        Path(interimOutputDir_actual).mkdir(parents=True, exist_ok=True)
        Path(finalSubjectOutputDir_actual).mkdir(parents=True, exist_ok=True)
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
            + " -s 0 -b 0 -o "
            + interimOutputDir_actual
        )

        with open(self.stdout_log, "a+") as out, open(
            self.stderr_log, "a+"
        ) as err:
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

        if not negatives_detected:
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

    def write(self):
        self.__write_csv(self.subjects, "processed_data.csv")
        self.__write_csv(self.neg_subjects, "QC_subjects_with_negative_intensities.csv")
        self.__write_csv(self.failing_subjects, "QC_subjects_with_bratspipeline_error.csv")

    def __write_csv(self, df: pd.DataFrame, output_file: str):
        if df.shape[0] > 0:
            df.to_csv(
                os.path.join(self.final_output_dir, output_file), index=False
            )


def main():
    parser = setup_parser()
    args = parser.parse_args()

    prep = Preparator(args.inputCSV, args.outputDir)
    prep.validate()
    prep.process_data()
    prep.write()


if __name__ == "__main__":
    if platform.system() == "Darwin":
        sys.exit("macOS is not supported")
    else:
        main()
