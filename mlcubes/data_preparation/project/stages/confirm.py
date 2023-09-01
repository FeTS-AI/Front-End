from typing import Union
import os
import yaml
import shutil

import pandas as pd
from pandas import DataFrame

from .dset_stage import DatasetStage
from .utils import get_id_tp
from .constants import TUMOR_MASK_FOLDER, INTERIM_FOLDER, FINAL_FOLDER


class ConfirmStage(DatasetStage):
    def __init__(
        self,
        data_csv: str,
        out_data_path: str,
        out_labels_path: str,
        prev_stage_path: str,
        backup_path: str,
        staging_folders: str,
    ):
        self.data_csv = data_csv
        self.out_data_path = out_data_path
        self.out_labels_path = out_labels_path
        self.prev_stage_path = prev_stage_path
        self.backup_path = backup_path
        self.staging_folders = staging_folders
        self.status_code = 7

    def get_name(self):
        return "Annotations Confirmation"

    def __get_input_data_path(self, index: Union[str, int]):
        id, tp = get_id_tp(index)
        path = os.path.join(self.prev_stage_path, FINAL_FOLDER, id, tp)
        return path

    def __get_input_label_path(self, index: Union[str, int]):
        id, tp = get_id_tp(index)
        path = os.path.join(self.prev_stage_path, INTERIM_FOLDER, id, tp, "reviewed")
        case = os.listdir(path)[0]

        return os.path.join(path, case)

    def __get_output_data_path(self, index: Union[str, int]):
        id, tp = get_id_tp(index)
        path = os.path.join(self.out_data_path, id, tp)
        return path

    def __get_output_label_path(self, index: Union[str, int]):
        id, tp = get_id_tp(index)
        path = os.path.join(self.out_labels_path, id, tp)
        filename = f"{id}_{tp}_final_seg.nii.gz"
        return os.path.join(path, filename)

    def __confirm(self, exact_match_percent: float) -> bool:
        exact_match_percent = round(exact_match_percent * 100, 2)
        msg = (
            f"We've identified {exact_match_percent}% of cases have not been modified "
            + "with respect to the baseline segmentation. Do you confirm this is intended? "
            + "[Y]/n"
        )

        # user_input = input(msg).lower()
        user_input = "y"
        return user_input == "y" or user_input == ""

    def __report_failure(self, report: DataFrame) -> DataFrame:
        # For this stage, failure is done when the user doesn't confirm
        # This means he probably wants to keep working on the data
        # And needs to know which rows are exact matches.
        # Because of this, failing this stage keeps the report intact
        return report

    def __process_row(self, row: pd.Series) -> pd.Series:
        """process a row by moving the required files
        to their respective locations, and removing any extra files

        Args:
            report (DataFrame): data preparation report

        Returns:
            DataFrame: modified data preparation report
        """
        print(row)
        index = row.name
        input_data_path = self.__get_input_data_path(index)
        input_label_path = self.__get_input_label_path(index)
        output_data_path = self.__get_output_data_path(index)
        output_label_path = self.__get_output_label_path(index)

        shutil.rmtree(output_data_path, ignore_errors=True)
        shutil.copytree(input_data_path, output_data_path)
        shutil.copy(input_label_path, output_label_path)

        row["status"] = self.status_code
        row["status_name"] = "ANNOTATION_CONFIRMED"
        row["data_path"] = output_data_path
        row["labels_path"] = output_label_path
        row["comment"] = ""
        return row

    def __cleanup_storage(self):
        for folder in self.staging_folders:
            staging_path = os.path.join(self.out_data_path, folder)
            shutil.rmtree(staging_path, ignore_errors=True)

    def should_run(self, report: DataFrame) -> bool:
        # Should run once all cases have been compared to the ground truth
        missing_voxels = report["num_changed_voxels"].isnull().values.any()

        return not missing_voxels

    def execute(self, report: DataFrame) -> DataFrame:
        exact_match_percent = (report["num_changed_voxels"] == 0).sum() / len(report)
        confirmed = self.__confirm(exact_match_percent)

        if not confirmed:
            return self.__report_failure(report)

        report = report.apply(self.__process_row, axis=1)
        # Remove all intermediary steps
        self.__cleanup_storage()

        return report
