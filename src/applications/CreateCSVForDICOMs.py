import os, argparse, sys, platform, posixpath, re
from pathlib import Path
from datetime import date
from tqdm import tqdm
import pandas as pd
import SimpleITK as sitk

from .constants import MODALITY_ID_DICT


def verify_dicom_folder(dicom_folder: str) -> (bool, str):
    """
    This function verifies that the folder is a valid DICOM folder. In the case of NIfTI file input, it will just verify if a 3D NIfTI is being passed in.

    Args:
        dicom_folder (str): The path to the DICOM folder or NIfTI file.

    Returns:
        bool: True if the folder is a valid DICOM folder or a 3D NIfTI file, False otherwise.
        str: The path to the first DICOM file in the folder if the folder is a valid DICOM folder, the NIfTI file itself otherwise.
        str: If verification failed, a message explaining the cause of failure is returned. Otherwise an empty string is returned.
    """

    if os.path.isdir(dicom_folder):
        series_IDs = sitk.ImageSeriesReader.GetGDCMSeriesIDs(dicom_folder)
        if not series_IDs:
            msg = (
                f"No valid or readable DICOM files were found in `{dicom_folder}`. "
                "Please convert your images using an external tool (such as dcm2niix) "
                "and try again by passing the NIfTI images in the CSV files instead of the DICOM images."
            )
            return False, None, msg

        if len(series_IDs) > 1:
            msg = (
                f"More than 1 DICOM series was found in `{dicom_folder}`, "
                "which usually means that multiple modalities or timepoints "
                "are present in that folder. Please ensure that only 1 modality is present per folder and try again."
            )
            return False, None, msg

        series_file_names = sitk.ImageSeriesReader.GetGDCMSeriesFileNames(
            dicom_folder, series_IDs[0]
        )
        series_reader = sitk.ImageSeriesReader()
        series_reader.SetFileNames(series_file_names)
        series_reader.MetaDataDictionaryArrayUpdateOn()
        series_reader.LoadPrivateTagsOn()
        image = series_reader.Execute()
    else:
        image = sitk.ReadImage(dicom_folder)
        series_file_names = [dicom_folder]

    if image.GetDimension() != 3:
        msg = (
            f"Only 3D volumes are supported, but the image in `{dicom_folder}` "
            f"appears to have {int(image.GetDimension())} dimensions. If your are scans "
            "have been encoded differently, please contact your vendor to convert them to 3D volumes."
        )
        return False, None, msg

    return True, series_file_names[0], ""


def setup_argparser():
    copyrightMessage = (
        "Contact: admin@fets.ai\n\n"
        + "This program is NOT FDA/CE approved and NOT intended for clinical use.\nCopyright (c) "
        + str(date.today().year)
        + " University of Pennsylvania. All rights reserved."
    )
    parser = argparse.ArgumentParser(
        prog="CreateCSVForDICOMS",
        formatter_class=argparse.RawTextHelpFormatter,
        description="This application creates the CSV for the DICOM folder structure.\n\n"
        + copyrightMessage,
    )
    parser.add_argument(
        "-inputDir",
        type=str,
        help="The absolute path to the input directory",
        required=True,
    )
    parser.add_argument(
        "-outputCSV",
        type=str,
        help="The output csv file name",
        required=True,
    )

    return parser.parse_args()


class CSVCreator:
    def __init__(self, inputDir: str, outputCSV: str):
        self.inputDir = inputDir
        self.outputCSV = outputCSV
        self.subject_timepoint_missing_modalities = []
        self.subject_timepoint_extra_modalities = []
        self.output_df_for_csv = pd.DataFrame(
            columns=["SubjectID", "Timepoint", "T1", "T1GD", "T2", "FLAIR"]
        )

    def process_data(self):
        for subject in tqdm(os.listdir(self.inputDir)):
            self.process_row(subject)

    def process_row(self, subject):
        inputDir = posixpath.normpath(self.inputDir)
        current_subject_dir = posixpath.join(inputDir, subject)

        if not os.path.isdir(current_subject_dir):
            return

        for timepoint in os.listdir(current_subject_dir):
            self.process_timepoint(timepoint, subject, current_subject_dir)

    def process_timepoint(self, timepoint, subject, subject_dir):
        timepoint_dir = posixpath.join(subject_dir, timepoint)
        if not os.path.isdir(timepoint_dir):
            return

        modality_folders = os.listdir(timepoint_dir)
        modality_folders = [folder for folder in modality_folders if not folder.startswith(".")]
        # check if there are missing modalities
        subject_tp = subject + "_" + timepoint
        if len(modality_folders) < 4:
            msg = (
                f"Less than 4 modalities where identified: {modality_folders}. "
                "Please ensure all modalities are present."
            )
            self.subject_timepoint_missing_modalities.append((subject_tp, msg))
            return
        # check if there are extra modalities
        if len(modality_folders) > 4:
            msg = (
                f"More than 4 modalities where identified: {modality_folders}. "
                "Please review your data and remove any unwanted files/folders."
            )
            self.subject_timepoint_extra_modalities.append((subject_tp, msg))
            return

        # goldilocks zone
        detected_modalities = {
            "T1": None,
            "T1GD": None,
            "T2": None,
            "FLAIR": None,
        }
        for modality in modality_folders:
            modality_path = posixpath.join(timepoint_dir, modality)
            modality_lower = modality.lower()
            modality_norm = re.sub('\.nii\.gz', '', modality_lower)
            for modality_to_check in MODALITY_ID_DICT:
                if detected_modalities[modality_to_check] is not None:
                    continue

                for modality_id in MODALITY_ID_DICT[modality_to_check]:
                    if modality_id != modality_norm:
                        continue

                    valid_dicom, first_dicom_file, msg = verify_dicom_folder(modality_path)
                    if valid_dicom:
                        detected_modalities[modality_to_check] = first_dicom_file
                        break
                    else:
                        subject_tp_mod = subject + "_" + timepoint + "_" + modality
                        self.subject_timepoint_missing_modalities.append((subject_tp_mod, msg))

        # check if any modalities are missing
        modalities_missing = False
        for modality in detected_modalities:
            if detected_modalities[modality] is None:
                modalities_missing = True
                subject_tp_mod = subject + "_" + timepoint + "_" + modality
                msg = (
                    f"A valid file/directory corresponding to modality {modality} could not be found. "
                    "Please ensure all modalities are included and are valid."
                )
                self.subject_timepoint_missing_modalities.append((subject_tp_mod, msg))

        if modalities_missing:
            return

        # if no modalities are missing, then add to the output csv
        dict_to_append = {
            "SubjectID": subject,
            "Timepoint": timepoint,
            "T1": detected_modalities["T1"],
            "T1GD": detected_modalities["T1GD"],
            "T2": detected_modalities["T2"],
            "FLAIR": detected_modalities["FLAIR"],
        }
        self.output_df_for_csv = pd.concat(
            [
                self.output_df_for_csv,
                pd.DataFrame(
                    [dict_to_append],
                    columns=[
                        "SubjectID",
                        "Timepoint",
                        "T1",
                        "T1GD",
                        "T2",
                        "FLAIR",
                    ],
                ),
            ],
        )

    def write(self):
        if self.output_df_for_csv.shape[0] > 0:
            if not (self.outputCSV.endswith(".csv")):
                self.outputCSV += ".csv"
            self.output_df_for_csv.to_csv(self.outputCSV, index=False)


def main(inputDir: str, outputCSV: str):
    inputDir = str(Path(inputDir).resolve())
    csv_creator = CSVCreator(inputDir, outputCSV)
    csv_creator.process_data()
    csv_creator.write()

    # print out the missing modalities
    missing = csv_creator.subject_timepoint_missing_modalities
    extra = csv_creator.subject_timepoint_extra_modalities
    if len(missing) > 0:
        print(
            "WARNING: The following subject timepoints are missing modalities: ",
            [subject for subject, _ in missing],
        )
    if len(extra) > 0:
        print(
            "WARNING: The following subject timepoints have extra modalities: ",
            [subject for subject, _ in extra],
        )

    print("Done!")


if __name__ == "__main__":
    args = setup_argparser()
    if platform.system().lower() == "darwin":
        sys.exit("macOS is not supported")
    else:
        main(args.inputDir, args.outputCSV)
