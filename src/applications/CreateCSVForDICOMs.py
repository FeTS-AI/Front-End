import os, argparse, sys, platform, posixpath
from pathlib import Path
from datetime import date
from tqdm import tqdm
import pandas as pd
import SimpleITK as sitk

from .constants import MODALITY_ID_DICT


def verify_dicom_folder(dicom_folder):
    series_IDs = sitk.ImageSeriesReader.GetGDCMSeriesIDs(dicom_folder)
    if not series_IDs:
        return False, None

    if len(series_IDs) > 1:
        return False, None

    series_file_names = sitk.ImageSeriesReader.GetGDCMSeriesFileNames(
        dicom_folder, series_IDs[0]
    )
    series_reader = sitk.ImageSeriesReader()
    series_reader.SetFileNames(series_file_names)
    series_reader.MetaDataDictionaryArrayUpdateOn()
    series_reader.LoadPrivateTagsOn()
    image_dicom = series_reader.Execute()
    if image_dicom.GetDimension() != 3:
        return False, None

    return True, series_file_names[0]


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
        # check if there are missing modalities
        if len(modality_folders) < 4:
            self.subject_timepoint_missing_modalities.append(subject + "_" + timepoint)
            return
        # check if there are extra modalities
        if len(modality_folders) > 4:
            self.subject_timepoint_extra_modalities.append(subject + "_" + timepoint)
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
            for modality_to_check in MODALITY_ID_DICT:
                if detected_modalities[modality_to_check] is not None:
                    continue

                for modality_id in MODALITY_ID_DICT[modality_to_check]:
                    if modality_id != modality_lower:
                        continue

                    valid_dicom, first_dicom_file = verify_dicom_folder(modality_path)
                    if valid_dicom:
                        detected_modalities[modality_to_check] = first_dicom_file
                        break
                    else:
                        self.subject_timepoint_missing_modalities.append(
                            subject + "_" + timepoint + "_" + modality
                        )

        # check if any modalities are missing
        modalities_missing = False
        for modality in detected_modalities:
            if detected_modalities[modality] is None:
                modalities_missing = True
                self.subject_timepoint_missing_modalities.append(
                    subject + "_" + timepoint + "_" + modality
                )

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
            missing,
        )
    if len(extra) > 0:
        print(
            "WARNING: The following subject timepoints have extra modalities: ",
            extra,
        )

    print("Done!")


if __name__ == "__main__":
    args = setup_argparser()
    if platform.system().lower() == "darwin":
        sys.exit("macOS is not supported")
    else:
        main(args.inputDir, args.outputCSV)
