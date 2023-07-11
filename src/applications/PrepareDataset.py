import os, argparse, sys, csv, platform, subprocess, shutil, posixpath, yaml
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


def main():
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

    args = parser.parse_args()

    assert os.path.exists(args.inputCSV), "Input CSV file not found"

    outputDir_qc = posixpath.join(args.outputDir, "DataForQC")
    outputDir_final = posixpath.join(args.outputDir, "DataForFeTS")

    for dir_to_create in [args.outputDir, outputDir_qc, outputDir_final]:
        Path(dir_to_create).mkdir(parents=True, exist_ok=True)

    bratsPipeline_exe = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "BraTSPipeline"
    )
    if platform.system().lower() == "windows":
        bratsPipeline_exe += ".exe"

    assert os.path.exists(
        bratsPipeline_exe
    ), "BraTS Pipeline executable not found, please contact admin@fets.ai for help."

    # only parse the headers here
    parsed_headers = parse_csv_header(args.inputCSV)

    # use pandas for this
    subjects_df = pd.read_csv(args.inputCSV)

    output_df_for_csv = pd.DataFrame(
        columns=["SubjectID", "Timepoint", "T1", "T1GD", "T2", "FLAIR"]
    )

    subjects_with_negatives = pd.DataFrame(
        columns=["SubjectID", "Timepoint", "Modality", "Count"]
    )
    subjects_with_bratspipeline_error = pd.DataFrame(columns=["SubjectID", "Timepoint"])

    # files to store the logs from BraTSPipeline
    preparedataset_stdout_log = os.path.join(
        args.outputDir, "preparedataset_stdout.txt"
    )
    preparedataset_stderr_log = os.path.join(
        args.outputDir, "preparedataset_stderr.txt"
    )

    dicom_tag_information_to_write_collaborator, dicom_tag_information_to_write_anon = (
        {},
        {},
    )

    for idx, row in tqdm(
        subjects_df.iterrows(),
        total=subjects_df.shape[0],
        desc="Preparing Dataset (1-10 min per subject)",
    ):
        subject_id = row[parsed_headers["ID"]]
        subject_id_timepoint = subject_id
        # create QC and Final output dirs for each subject
        interimOutputDir_subject = posixpath.join(outputDir_qc, subject_id_timepoint)
        Path(interimOutputDir_subject).mkdir(parents=True, exist_ok=True)
        finalSubjectOutputDir_subject = posixpath.join(
            outputDir_final, subject_id_timepoint
        )
        Path(finalSubjectOutputDir_subject).mkdir(parents=True, exist_ok=True)
        interimOutputDir_actual = interimOutputDir_subject
        finalSubjectOutputDir_actual = finalSubjectOutputDir_subject
        # per the data ingestion step, we are creating a new folder called timepoint, can join timepoint to subjectid if needed
        if parsed_headers["Timepoint"] is not None:
            timepoint = row[parsed_headers["Timepoint"]]
            subject_id_timepoint += "_" + timepoint
            interimOutputDir_actual = posixpath.join(
                interimOutputDir_subject, timepoint
            )
            finalSubjectOutputDir_actual = posixpath.join(
                finalSubjectOutputDir_subject, timepoint
            )

        dicom_tag_information_to_write_collaborator[subject_id_timepoint] = {}
        dicom_tag_information_to_write_anon[str(idx)] = {}

        for modality in ["T1", "T1GD", "T2", "FLAIR"]:
            tags_from_modality = _get_relevant_dicom_tags(row[parsed_headers[modality]])
            dicom_tag_information_to_write_collaborator[subject_id_timepoint][
                modality
            ] = tags_from_modality
            dicom_tag_information_to_write_anon[str(idx)][modality] = tags_from_modality

        Path(interimOutputDir_actual).mkdir(parents=True, exist_ok=True)
        Path(finalSubjectOutputDir_actual).mkdir(parents=True, exist_ok=True)
        # if files already exist in DataForQC, then copy to DataForFeTS, and if files exist in DataForFeTS, then skip
        runBratsPipeline, outputs = copyFilesToCorrectLocation(
            interimOutputDir_actual, finalSubjectOutputDir_actual, subject_id_timepoint
        )

        negatives_detected = False
        # run the pipeline if needed
        if runBratsPipeline:
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
            with open(preparedataset_stdout_log, "a+") as out, open(
                preparedataset_stderr_log, "a+"
            ) as err:
                # save the command for debugging
                out.write(f"***\n{command}\n***")
                err.write(f"***\n{command}\n***")
                subprocess.Popen(command, stdout=out, stderr=err, shell=True).wait()

            runBratsPipeline, outputs = copyFilesToCorrectLocation(
                interimOutputDir_actual,
                finalSubjectOutputDir_actual,
                subject_id_timepoint,
            )
            # if there are any errors, then store the subjectid and timepoint
            if runBratsPipeline:
                subjects_with_bratspipeline_error = pd.concat(
                    [
                        subjects_with_bratspipeline_error,
                        pd.DataFrame(
                            {
                                "SubjectID": subject_id,
                                "Timepoint": timepoint,
                            },
                            index=[0],
                        ),
                    ]
                )
            else:
                for modality in ["T1", "T1GD", "T2", "FLAIR"]:
                    count = _read_image_with_min_check(outputs[modality])
                    # if there are any negative values, then store the subjectid, timepoint, modality and count of negative values
                    if count > 0:
                        subjects_with_negatives = pd.concat(
                            [
                                subjects_with_negatives,
                                pd.DataFrame(
                                    {
                                        "SubjectID": subject_id,
                                        "Timepoint": timepoint,
                                        "Modality": modality,
                                        "Count": count,
                                    },
                                    index=[0],
                                ),
                            ]
                        )
                        negatives_detected = True
        # store the outputs in a dictionary when there are no errors
        if not negatives_detected:
            output_df_for_csv = pd.concat(
                [
                    output_df_for_csv,
                    pd.DataFrame(
                        {
                            "SubjectID": subject_id,
                            "Timepoint": timepoint,
                            "T1": outputs["T1"],
                            "T1GD": outputs["T1GD"],
                            "T2": outputs["T2"],
                            "FLAIR": outputs["FLAIR"],
                        },
                        index=[0],
                    ),
                ]
            )

    # write the output file
    if output_df_for_csv.shape[0] > 0:
        output_df_for_csv.to_csv(
            os.path.join(outputDir_final, "processed_data.csv"), index=False
        )

    # write the QC files
    if subjects_with_negatives.shape[0] > 0:
        subjects_with_negatives.to_csv(
            os.path.join(outputDir_final, "QC_subjects_with_negative_intensities.csv"),
            index=False,
        )
    if subjects_with_bratspipeline_error.shape[0] > 0:
        subjects_with_bratspipeline_error.to_csv(
            os.path.join(outputDir_final, "QC_subjects_with_bratspipeline_error.csv"),
            index=False,
        )

    # write the dicom tag information
    with open(
        os.path.join(outputDir_final, "dicom_tag_information_collaborator.yaml"), "w"
    ) as f:
        yaml.safe_dump(
            dicom_tag_information_to_write_collaborator, f, allow_unicode=True
        )

    with open(
        os.path.join(outputDir_final, "dicom_tag_information_anon.yaml"), "w"
    ) as f:
        yaml.safe_dump(dicom_tag_information_to_write_anon, f, allow_unicode=True)


if __name__ == "__main__":
    if platform.system() == "Darwin":
        sys.exit("macOS is not supported")
    else:
        main()
