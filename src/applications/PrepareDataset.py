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
from GANDLF.cli import main_run
from LabelFusion.wrapper import fuse_images

# check against all these modality ID strings with extensions
MODALITY_ID_DICT = {
    "T1": ["t1", "t1pre", "t1precontrast"],
    "T1GD": ["t1ce", "t1gd", "t1post", "t1postcontrast", "t1gallodinium", "t1c"],
    "T2": ["t2"],
    "FLAIR": ["flair", "fl", "t2flair"],
}
MODALITIES_LIST = list(MODALITY_ID_DICT.keys())
SUBJECT_NAMES = {"patientid", "subjectid", "subject", "subid"}
TIMEPOINT_NAMES = {"timepoint", "tp", "time", "series", "subseries"}
INPUT_FILENAMES = {
    "T1": "T1_to_SRI.nii.gz",
    "T1GD": "T1CE_to_SRI.nii.gz",
    "T2": "T2_to_SRI.nii.gz",
    "FLAIR": "FL_to_SRI.nii.gz",
}

GANDLF_DF_COLUMNS = ["SubjectID", "Channel_0"]

INTERIM_FOLDER = "DataForQC"
FINAL_FOLDER = "DataForFeTS"
TUMOR_MASK_FOLDER = "TumorMasksForQC"
TESTING_FOLDER = "testing"
REORIENTED_FOLDER = "reoriented"

BRAIN_FILENAME = "gandlf_brain_extraction.csv"
TUMOR_FILENAME = "gandlf_tumor_segmentation.csv"
SUBJECTS_FILENAME = "processed_data.csv"
NEG_SUBJECTS_FILENAME = "QC_subjects_with_negative_intensities.csv"
FAIL_SUBJECTS_FILENAME = "QC_subjects_with_bratspipeline_error.csv"
DICOM_ANON_FILENAME = "dicom_tag_information_to_write_anon.yaml"
DICOM_COLLAB_FILENAME = "dicom_tag_information_to_write_collab.yaml"
STDOUT_FILENAME = "preparedataset_stdout.txt"
STDERR_FILENAME = "preparedataset_stderr.txt"

EXEC_NAME = "BraTSPipeline"


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
    parser.add_argument(
        "-executablePath",
        type=str,
        help="The path to the BraTSPipeline executable. If not given, will infer from current script's location.",
        nargs="?",
        const=None,
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


def _save_screenshot(
    input_images: dict, output_filename: str = None, input_mask: str = None
) -> None:
    """
    This function saves the screenshot of the input images and mask.

    Args:
        input_images (dict): The input multi-modal images.
        output_filename (str, optional): The output filename to save the screenshot. Defaults to None.
        input_mask (str, optional): The input mask filename. Defaults to None.
    """
    # save the screenshot
    images = (",").join(
        [
            input_images["T1"],
            input_images["T1GD"],
            input_images["T2"],
            input_images["FLAIR"],
        ]
    )
    ylabels = (",").join(MODALITIES_LIST)

    figure_generator(
        input_images=images,
        ylabels=ylabels,
        output=output_filename,
        input_mask=input_mask,
        flip_sagittal=True,
        flip_coronal=True,
    )


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


def _parse_csv_header(filename):
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
                if temp in SUBJECT_NAMES:
                    headers["ID"] = col
                elif temp in TIMEPOINT_NAMES:
                    headers["Timepoint"] = col
                else:
                    for key in MODALITY_ID_DICT.keys():
                        if temp in MODALITY_ID_DICT[key]:
                            headers[key] = col
                            break

    if "Timepoint" not in headers:
        headers["Timepoint"] = None
    return headers


def _copy_files_to_correct_location(interimOutputDir, finalSubjectOutputDir, subjectID):
    """
    This function copies the intermediate files and final outputs to correct location and if these are absent, returns a bool flag stating that brats pipeline needs to run again

    Args:
        interimOutputDir (str): The interim output directory.
        finalSubjectOutputDir (str): The final subject output directory.
        subjectID (str): The subject ID.

    Returns:
        bool, dict: The flag stating whether brats pipeline needs to run again and the output files in the expected location.
    """

    # copy files to correct location for inference and training
    runBratsPipeline = False
    input_files = {
        k: posixpath.join(interimOutputDir, v) for k, v in INPUT_FILENAMES.items()
    }
    expected_outputs = get_expected_outputs(subjectID, finalSubjectOutputDir)

    for key in input_files.keys():
        if not os.path.exists(expected_outputs[key]):
            if os.path.exists(input_files[key]):
                shutil.copyfile(input_files[key], expected_outputs[key])
            else:
                runBratsPipeline = True

    return runBratsPipeline, expected_outputs


def get_expected_outputs(subjectID: str, output_dir: str) -> dict:
    expected_outputs = {
        "ID": subjectID,
        "T1": posixpath.join(output_dir, subjectID + "_t1.nii.gz"),
        "T1GD": posixpath.join(output_dir, subjectID + "_t1ce.nii.gz"),
        "T2": posixpath.join(output_dir, subjectID + "_t2.nii.gz"),
        "FLAIR": posixpath.join(output_dir, subjectID + "_flair.nii.gz"),
    }
    return expected_outputs


def get_brain_mask_files(subject_id, output_dir) -> dict:
    files = {}
    for modality in MODALITIES_LIST:
        files[modality] = posixpath.join(
            output_dir,
            f"{subject_id}_brain_{modality}.nii.gz",
        )
    return files


def _run_brain_extraction_using_gandlf(
    subject_id: str,
    input_oriented_images: dict,
    models_to_infer: str,
    base_output_dir: str,
) -> sitk.Image:
    """
    This function runs brain extraction using gandlf.

    Args:
        subject_id (str): The subject ID.
        input_oriented_images (dict): The input oriented images.
        models_to_infer (str): The models to infer, comma-separated.
        base_output_dir (str): The base output directory.

    Returns:
        sitk.Image: The fused brain mask.
    """
    df_for_gandlf = pd.DataFrame(columns=GANDLF_DF_COLUMNS)
    for key in MODALITIES_LIST:
        current_modality = {
            "SubjectID": subject_id + "_" + key,
            "Channel_0": input_oriented_images[key],
        }
        df_for_gandlf = pd.concat(
            [df_for_gandlf, pd.DataFrame(current_modality, index=[0])]
        )
    data_path = posixpath.join(base_output_dir, BRAIN_FILENAME)
    df_for_gandlf.to_csv(
        data_path,
        index=False,
    )

    models_to_run = models_to_infer.split(",")

    model_counter = 0
    images_for_fusion = []
    for model_dir in models_to_run:
        model_output_dir = posixpath.join(
            base_output_dir, "model_" + str(model_counter)
        )
        file_list = os.listdir(model_dir)
        for file in file_list:
            if file.endswith(".yaml") or file.endswith(".yml"):
                config_file = posixpath.join(model_dir, file)
                break

        # ensure the openvino version is used
        parameters = yaml.safe_load(open(config_file, "r"))
        parameters["model"]["type"] = "openvino"
        yaml.safe_dump(parameters, open(config_file, "w"))

        main_run(
            data_csv=data_path,
            config_file=config_file,
            model_dir=model_dir,
            train_mode=False,
            device="cpu",
            resume=False,
            reset=False,
            output_dir=model_output_dir,
        )

        modality_outputs = os.listdir(posixpath.join(model_output_dir, TESTING_FOLDER))
        for modality in modality_outputs:
            modality_output_dir = posixpath.join(modality_outputs, modality)
            files_in_modality = os.listdir(modality_output_dir)
            for file in files_in_modality:
                if file.endswith(".nii.gz"):
                    file_path = posixpath.join(modality_output_dir, file)
                    shutil.copyfile(
                        file_path,
                        posixpath.join(
                            base_output_dir,
                            f"brainMask_{model_counter}_{modality}.nii.gz",
                        ),
                    )
                    images_for_fusion.append(sitk.ReadImage(file_path, sitk.sitkInt16))
        model_counter += 1

    return fuse_images(images_for_fusion, "staple", [0, 1])


def _run_tumor_segmentation_using_gandlf(
    subject_id: str,
    input_oriented_brain_images: dict,
    models_to_infer: str,
    base_output_dir: str,
) -> sitk.Image:
    """
    This function runs tumor segmentation using gandlf.

    Args:
        subject_id (str): The subject ID.
        input_oriented_brain_images (dict): The input oriented brain images.
        models_to_infer (str): The models to infer, comma-separated.
        base_output_dir (str): The base output directory.

    Returns:
        sitk.Image: The fused tumor mask.
    """
    df_for_gandlf = pd.DataFrame(columns=GANDLF_DF_COLUMNS)
    current_subject = {"SubjectID": subject_id}
    channel_idx = 0
    # todo: confirm the order for modalities
    for key in MODALITIES_LIST:
        current_subject = {
            f"Channel_{channel_idx}": input_oriented_brain_images[key],
        }
        channel_idx += 1
    df_for_gandlf = pd.DataFrame(current_subject, index=[0])
    data_path = posixpath.join(base_output_dir, TUMOR_FILENAME)
    df_for_gandlf.to_csv(
        data_path,
        index=False,
    )

    models_to_run = models_to_infer.split(",")

    model_counter = 0
    images_for_fusion = []
    mask_output_dir = posixpath.join(base_output_dir, TUMOR_MASK_FOLDER)
    for model_dir in models_to_run:
        model_output_dir = posixpath.join(
            base_output_dir, "model_" + str(model_counter)
        )
        file_list = os.listdir(model_dir)
        for file in file_list:
            if file.endswith(".yaml") or file.endswith(".yml"):
                config_file = posixpath.join(model_dir, file)
                break

        # ensure the openvino version is used
        parameters = yaml.safe_load(open(config_file, "r"))
        parameters["model"]["type"] = "openvino"
        yaml.safe_dump(parameters, open(config_file, "w"))

        main_run(
            data_csv=data_path,
            config_file=config_file,
            model_dir=model_dir,
            train_mode=False,
            device="cpu",
            resume=False,
            reset=False,
            output_dir=model_output_dir,
        )

        subject_model_output_dir = os.listdir(
            posixpath.join(model_output_dir, TESTING_FOLDER)
        )
        for subject in subject_model_output_dir:
            subject_output_dir = posixpath.join(subject_model_output_dir, subject)
            files_in_modality = os.listdir(subject_output_dir)
            for file in files_in_modality:
                if file.endswith(".nii.gz"):
                    file_path = posixpath.join(subject_output_dir, file)
                    shutil.copyfile(
                        file_path,
                        posixpath.join(
                            mask_output_dir,
                            f"{subject_id}_tumorMask_model-{model_counter}.nii.gz",
                        ),
                    )
                    images_for_fusion.append(sitk.ReadImage(file_path, sitk.sitkInt16))
        model_counter += 1

    tumor_class_list = [0, 1, 2, 3, 4]

    tumor_masks_to_return = images_for_fusion

    if len(images_for_fusion) > 1:
        for fusion_type in ["staple", "simple", "voting"]:
            fused_mask = fuse_images(images_for_fusion, fusion_type, tumor_class_list)
            fused_mask_file = posixpath.join(
                mask_output_dir,
                f"{subject_id}_tumorMask_fused-{fusion_type}.nii.gz",
            )
            sitk.WriteImage(fused_mask, fused_mask_file)
            tumor_masks_to_return.append(fused_mask_file)

    return tumor_masks_to_return


class Preparator:
    def __init__(self, input_csv: str, output_dir: str, executablePath: str):
        self.input_csv = input_csv
        self.input_dir = str(Path(input_csv).parent)
        self.output_dir = os.path.normpath(output_dir)
        self.interim_output_dir = posixpath.join(self.output_dir, INTERIM_FOLDER)
        self.final_output_dir = posixpath.join(self.output_dir, FINAL_FOLDER)
        self.subjects_file = posixpath.join(self.final_output_dir, SUBJECTS_FILENAME)
        self.neg_subjects_file = posixpath.join(
            self.final_output_dir, NEG_SUBJECTS_FILENAME
        )
        self.failing_subjects_file = posixpath.join(
            self.final_output_dir, FAIL_SUBJECTS_FILENAME
        )
        self.dicom_tag_information_to_write_anon_file = posixpath.join(
            self.final_output_dir, DICOM_ANON_FILENAME
        )
        self.dicom_tag_information_to_write_collab_file = posixpath.join(
            self.final_output_dir, DICOM_COLLAB_FILENAME
        )
        self.__init_out_dfs()
        self.stdout_log = posixpath.join(self.output_dir, STDOUT_FILENAME)
        self.stderr_log = posixpath.join(self.output_dir, STDERR_FILENAME)
        self.dicom_tag_information_to_write_collab = {}
        self.dicom_tag_information_to_write_anon = {}
        self.brats_pipeline_exe = executablePath
        if self.brats_pipeline_exe is None:
            self.brats_pipeline_exe = os.path.join(
                Path(__file__).parent.resolve(), EXEC_NAME
            )

        if platform.system() == "Windows":
            if not self.brats_pipeline_exe.endswith(".exe"):
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

        assert (
            shutil.which(self.brats_pipeline_exe) is not None
        ), "BraTS Pipeline executable not found, please contact admin@fets.ai for help."

    def process_data(self):
        items = self.subjects_df.iterrows()
        total = self.subjects_df.shape[0]
        pbar = tqdm(range(total), desc="Preparing Dataset (1-10 min per subject)")
        for idx, row in items:
            self.process_row(idx, row, pbar)

    def process_row(self, idx: int, row: pd.Series, pbar: tqdm):
        self.convert_to_dicom(idx, row, pbar)
        self.extract_brain(row, pbar)
        self.extract_tumor(row, pbar)

    def __get_row_information(self, row: pd.Series):
        subject_id = row[self.parsed_headers["ID"]]
        subject_id_timepoint = subject_id

        # create QC and Final output dirs for each subject
        interimOutputDir_actual = posixpath.join(
            self.interim_output_dir, subject_id_timepoint
        )
        finalSubjectOutputDir_actual = posixpath.join(
            self.final_output_dir, subject_id_timepoint
        )

        return subject_id, interimOutputDir_actual, finalSubjectOutputDir_actual

    def convert_to_dicom(self, idx: int, row: pd.Series, pbar: tqdm):
        parsed_headers = self.parsed_headers
        bratsPipeline_exe = self.brats_pipeline_exe

        (
            subject_id,
            interimOutputDir_actual,
            finalSubjectOutputDir_actual,
        ) = self.__get_row_information(row)
        subject_id_timepoint = subject_id

        # per the data ingestion step, we are creating a new folder called timepoint, can join timepoint to subjectid if needed
        if parsed_headers["Timepoint"] is not None:
            timepoint = row[parsed_headers["Timepoint"]]
            subject_id_timepoint += "_" + timepoint
            interimOutputDir_actual = posixpath.join(interimOutputDir_actual, timepoint)
            finalSubjectOutputDir_actual = posixpath.join(
                finalSubjectOutputDir_actual, timepoint
            )

        # create QC and Final output dirs for each subject
        Path(interimOutputDir_actual).mkdir(parents=True, exist_ok=True)
        Path(finalSubjectOutputDir_actual).mkdir(parents=True, exist_ok=True)

        pbar.set_description(f"Processing {subject_id_timepoint}")

        # get the relevant dicom tags
        self.dicom_tag_information_to_write_collab[subject_id_timepoint] = {}
        self.dicom_tag_information_to_write_anon[str(idx)] = {}
        for modality in MODALITIES_LIST:
            tags_from_modality = _get_relevant_dicom_tags(row[parsed_headers[modality]])
            self.dicom_tag_information_to_write_collab[subject_id_timepoint][
                modality
            ] = tags_from_modality
            with open(
                os.path.join(
                    interimOutputDir_actual, f"dicom_tag_information_{modality}.yaml"
                ),
                "w",
            ) as f:
                yaml.safe_dump(tags_from_modality, f, allow_unicode=True)
            self.dicom_tag_information_to_write_anon[str(idx)][
                modality
            ] = tags_from_modality

        interimOutputDir_actual_reoriented = posixpath.join(
            interimOutputDir_actual, REORIENTED_FOLDER
        )
        Path(interimOutputDir_actual_reoriented).mkdir(parents=True, exist_ok=True)
        # if files already exist in DataForQC, then copy to "reorient" folder, and if files exist in "reorient" folder, then skip
        runBratsPipeline, _ = _copy_files_to_correct_location(
            interimOutputDir_actual,
            interimOutputDir_actual_reoriented,
            subject_id_timepoint,
        )

        # check if the files exist already, if so, skip
        if runBratsPipeline:
            pbar.set_description(f"Running BraTSPipeline")

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
                + " -s 0 -o "
                + interimOutputDir_actual
            )

            with open(self.stdout_log, "a+") as out, open(self.stderr_log, "a+") as err:
                out.write(f"***\n{command}\n***")
                err.write(f"***\n{command}\n***")
                subprocess.Popen(command, stdout=out, stderr=err, shell=True).wait()

        runBratsPipeline, outputs_reoriented = _copy_files_to_correct_location(
            interimOutputDir_actual,
            interimOutputDir_actual_reoriented,
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
        for modality in MODALITIES_LIST:
            count = _read_image_with_min_check(outputs_reoriented[modality])
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
            "T1": outputs_reoriented["T1"],
            "T1GD": outputs_reoriented["T1GD"],
            "T2": outputs_reoriented["T2"],
            "FLAIR": outputs_reoriented["FLAIR"],
        }
        subject = pd.DataFrame(subject_data, index=[0])
        self.subjects = pd.concat(
            [
                self.subjects,
                subject,
            ]
        )

        pbar.set_description(f"Saving screenshot")

        screenshot_path = posixpath.join(
            interimOutputDir_actual_reoriented,
            f"{subject_id_timepoint}_summary_coregistration.png",
        )
        # save the screenshot
        _save_screenshot(outputs_reoriented, screenshot_path)

        if os.path.exists(screenshot_path):
            shutil.copyfile(
                screenshot_path,
                posixpath.join(
                    interimOutputDir_actual,
                    f"{subject_id_timepoint}_summary_coregistration.png",
                ),
            )

    def extract_brain(self, row: pd.Series, pbar: tqdm):
        (
            subject_id_timepoint,
            interimOutputDir_actual,
            finalSubjectOutputDir_actual,
        ) = self.__get_row_information(row)
        interimOutputDir_actual_reoriented = posixpath.join(
            interimOutputDir_actual, REORIENTED_FOLDER
        )
        outputs_reoriented = get_expected_outputs(
            subject_id_timepoint, interimOutputDir_actual_reoriented
        )

        pbar.set_description(f"Brain Extraction")

        brain_mask = _run_brain_extraction_using_gandlf(
            subject_id_timepoint,
            outputs_reoriented,
            interimOutputDir_actual
            + ","
            + interimOutputDir_actual,  # todo: this needs to be changed appropriately
            interimOutputDir_actual,
        )
        brain_mask_path = posixpath.join(
            interimOutputDir_actual, "brainMask_fused.nii.gz"
        )
        sitk.WriteImage(brain_mask, brain_mask_path)

        # this is to ensure that the mask and reoriented images are in the same byte order
        brain_mask = sitk.Cast(brain_mask, sitk.sitkFloat32)
        input_for_tumor_models = get_brain_mask_files(
            subject_id_timepoint, finalSubjectOutputDir_actual
        )
        for modality in MODALITIES_LIST:
            image = sitk.ReadImage(outputs_reoriented[modality])
            masked_image = sitk.Mask(image, brain_mask)
            file_to_save = input_for_tumor_models[modality]
            sitk.WriteImage(masked_image, file_to_save)

        # save the screenshot
        _save_screenshot(
            input_for_tumor_models,
            posixpath.join(
                interimOutputDir_actual,
                f"{subject_id_timepoint}_summary_brain-extraction.png",
            ),
            brain_mask_path,
        )

    def extract_tumor(self, row: pd.Series, pbar: tqdm):
        (
            subject_id_timepoint,
            interimOutputDir_actual,
            finalSubjectOutputDir_actual,
        ) = self.__get_row_information(row)
        input_for_tumor_models = get_brain_mask_files(
            subject_id_timepoint, finalSubjectOutputDir_actual
        )

        pbar.set_description(f"Brain Tumor Segmentation")

        tumor_masks_for_qc = _run_tumor_segmentation_using_gandlf(
            subject_id_timepoint,
            input_for_tumor_models,
            interimOutputDir_actual
            + ","
            + interimOutputDir_actual,  # todo: this needs to be changed appropriately
            interimOutputDir_actual,
        )

        tumor_mask_idx = 0
        for tumor_mask in tumor_masks_for_qc:
            # save the screenshot
            _save_screenshot(
                input_for_tumor_models,
                posixpath.join(
                    interimOutputDir_actual,
                    f"{subject_id_timepoint}_summary_tumor-segmentation_model-{tumor_mask_idx}.png",
                ),
                tumor_mask,
            )
            tumor_mask_idx += 1

        with open(self.stdout_log, "a+") as f:
            f.write(f"***\nTumor Masks For QC:\n{tumor_masks_for_qc}\n***")

    def write(self):
        if self.subjects.shape[0]:
            self.subjects.to_csv(self.subjects_file, index=False)
        if self.neg_subjects.shape[0]:
            self.neg_subjects.to_csv(self.neg_subjects_file, index=False)
        if self.failing_subjects.shape[0]:
            self.failing_subjects.to_csv(self.failing_subjects_file, index=False)
        with open(self.dicom_tag_information_to_write_collab_file, "w") as f:
            yaml.safe_dump(
                self.dicom_tag_information_to_write_collab, f, allow_unicode=True
            )
        with open(self.dicom_tag_information_to_write_anon_file, "w") as f:
            yaml.safe_dump(
                self.dicom_tag_information_to_write_anon, f, allow_unicode=True
            )

    def read(self):
        self.parsed_headers = _parse_csv_header(self.input_csv)
        self.subjects_df = pd.read_csv(self.input_csv)
        if os.path.exists(self.subjects_file):
            self.subjects = pd.read_csv(self.subjects_file)
        if os.path.exists(self.neg_subjects_file):
            self.neg_subjects = pd.read_csv(self.neg_subjects_file)
        if os.path.exists(self.failing_subjects_file):
            self.failing_subjects = pd.read_csv(self.failing_subjects_file)
        if os.path.exists(self.dicom_tag_information_to_write_collab_file):
            with open(self.dicom_tag_information_to_write_collab_file, "r") as f:
                self.dicom_tag_information_to_write_collab = yaml.safe_load(f)
        if os.path.exists(self.dicom_tag_information_to_write_anon_file):
            with open(self.dicom_tag_information_to_write_anon_file, "r") as f:
                self.dicom_tag_information_to_write_anon = yaml.safe_load(f)


def main():
    parser = setup_parser()
    args = parser.parse_args()

    prep = Preparator(args.inputCSV, args.outputDir, args.executablePath)
    prep.validate()
    prep.read()
    prep.process_data()
    prep.write()


if __name__ == "__main__":
    if platform.system().lower() == "darwin":
        sys.exit("macOS is not supported")
    else:
        main()
