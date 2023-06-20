import os, argparse, sys, platform
from datetime import date
from tqdm import tqdm
import pandas as pd
import SimpleITK as sitk


# check against all these modality ID strings with extensions
modality_id_dict = {
    "T1": ["t1", "t1pre", "t1precontrast"],
    "T1GD": ["t1ce", "t1gd", "t1post", "t1postcontrast", "t1gallodinium", "t1c"],
    "T2": ["t2"],
    "FLAIR": ["flair", "fl", "t2flair"],
}


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


def main():
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

    args = parser.parse_args()

    output_df_for_csv = pd.DataFrame(
        columns=["SubjectID", "Timepoint", "T1", "T1GD", "T2", "FLAIR"]
    )

    subject_timepoint_missing_modalities, subject_timepoint_extra_modalities = [], []

    for subject in tqdm(os.listdir(args.inputDir)):
        current_subject_dir = os.path.join(args.inputDir, subject)
        if os.path.isdir(current_subject_dir):
            for timepoint in os.listdir(current_subject_dir):
                current_subject_timepoint_dir = os.path.join(
                    current_subject_dir, timepoint
                )
                if os.path.isdir(current_subject_timepoint_dir):
                    modality_folders = os.listdir(current_subject_timepoint_dir)
                    # check if there are missing modalities
                    if len(modality_folders) < 4:
                        subject_timepoint_missing_modalities.append(
                            subject + "_" + timepoint
                        )
                    # check if there are extra modalities
                    elif len(modality_folders) > 4:
                        subject_timepoint_extra_modalities.append(
                            subject + "_" + timepoint
                        )
                    # goldilocks zone
                    else:
                        detected_modalities = {
                            "T1": None,
                            "T1GD": None,
                            "T2": None,
                            "FLAIR": None,
                        }
                        for modality in modality_folders:
                            modality_path = os.path.join(
                                current_subject_timepoint_dir, modality
                            )
                            modality_lower = modality.lower()
                            for modality_to_check in modality_id_dict:
                                if detected_modalities[modality_to_check] is None:
                                    for modality_id in modality_id_dict[
                                        modality_to_check
                                    ]:
                                        if modality_id in modality_lower:
                                            (
                                                valid_dicom,
                                                first_dicom_file,
                                            ) = verify_dicom_folder(modality_path)
                                            if valid_dicom:
                                                detected_modalities[
                                                    modality_to_check
                                                ] = first_dicom_file
                                                break
                                            else:
                                                subject_timepoint_missing_modalities.append(
                                                    subject
                                                    + "_"
                                                    + timepoint
                                                    + "_"
                                                    + modality
                                                )

                        # check if any modalities are missing
                        modalities_missing = False
                        for modality in detected_modalities:
                            if detected_modalities[modality] is None:
                                modalities_missing = True
                                subject_timepoint_missing_modalities.append(
                                    subject + "_" + timepoint + "_" + modality
                                )

                        # if no modalities are missing, then add to the output csv
                        if not modalities_missing:
                            dict_to_append = {
                                "SubjectID": subject,
                                "Timepoint": timepoint,
                                "T1": detected_modalities["T1"],
                                "T1GD": detected_modalities["T1GD"],
                                "T2": detected_modalities["T2"],
                                "FLAIR": detected_modalities["FLAIR"],
                            }
                            # output_df_for_csv = pd.concat(
                            #     [
                            #         dict_to_append,
                            #         pd.DataFrame(
                            #             columns=[
                            #                 "SubjectID",
                            #                 "Timepoint",
                            #                 "T1",
                            #                 "T1GD",
                            #                 "T2",
                            #                 "FLAIR",
                            #             ],
                            #         ),
                            #     ],
                            #     ignore_index=True,
                            # )
                            output_df_for_csv = pd.concat(
                                [
                                    output_df_for_csv,
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

    # write the output csv
    if output_df_for_csv.shape[0] > 0:
        if not(args.outputCSV.endswith(".csv")):
            args.outputCSV += ".csv"
        output_df_for_csv.to_csv(args.outputCSV, index=False)

    # print out the missing modalities
    if len(subject_timepoint_missing_modalities) > 0:
        print(
            "WARNING: The following subject timepoints are missing modalities: ",
            subject_timepoint_missing_modalities,
        )
    if len(subject_timepoint_extra_modalities) > 0:
        print(
            "WARNING: The following subject timepoints have extra modalities: ",
            subject_timepoint_extra_modalities,
        )

    print("Done!")


if __name__ == "__main__":
    if platform.system() == "Darwin":
        sys.exit("macOS is not supported")
    else:
        main()
