from typing import Union, Tuple
import pandas as pd
import os
import shutil

from .row_stage import RowStage
from .constants import TUMOR_MASK_FOLDER, INTERIM_FOLDER
from .utils import get_id_tp, update_row_with_dict, set_files_read_only


class ManualStage(RowStage):
    def __init__(
        self, data_csv: str, out_path: str, prev_stage_path: str, backup_path: str
    ):
        self.data_csv = data_csv
        self.out_path = out_path
        self.prev_stage_path = prev_stage_path
        self.backup_path = backup_path
        self.status_code = 5

    def get_name(self):
        return "Manual review"

    def __get_input_path(self, index: Union[str, int]):
        id, tp = get_id_tp(index)
        path = os.path.join(
            self.prev_stage_path, INTERIM_FOLDER, id, tp, TUMOR_MASK_FOLDER
        )
        return path

    def __get_under_review_path(self, index: Union[str, int]):
        id, tp = get_id_tp(index)
        path = os.path.join(self.out_path, INTERIM_FOLDER, id, tp, "under_review")
        return path

    def __get_output_path(self, index: Union[str, int]):
        id, tp = get_id_tp(index)
        path = os.path.join(self.out_path, INTERIM_FOLDER, id, tp, "reviewed")
        return path

    def __get_backup_path(self, index: Union[str, int]):
        id, tp = get_id_tp(index)
        path = os.path.join(self.backup_path, id, tp, TUMOR_MASK_FOLDER)
        return path

    def __report_success(
        self, index: Union[str, int], report: pd.DataFrame
    ) -> pd.DataFrame:
        out_path = self.__get_output_path(index)
        report_data = {
            "status": 5,
            "status_name": "MANUAL_REVIEW_COMPLETED",
            "comment": "",
            "data_path": out_path,
            "labels_path": "",
        }
        update_row_with_dict(report, report_data, index)
        return report

    def __report_step_missing(
        self, index: Union[str, int], report: pd.DataFrame
    ) -> pd.DataFrame:
        in_path = self.__get_input_path(index)
        out_path = self.__get_output_path(index)
        under_review_path = self.__get_under_review_path(index)
        msg = (
            f"You may find baseline segmentations inside {in_path}. "
            + f"Please inspect those segmentations and move the best one to {under_review_path}. "
            + "Make the necessary corrections to the generated segmentations with your desired tool, "
            + f"and once you're done, move the finalized file to {out_path}."
        )

        report_data = {
            "status": -self.status_code,
            "status_name": "MANUAL_REVIEW_REQUIRED",
            "comment": msg,
            "data_path": in_path,
            "labels_path": "",
        }
        update_row_with_dict(report, report_data, index)
        return report

    def __report_multiple_cases_error(
        self, index: Union[str, int], report: pd.DataFrame
    ) -> pd.DataFrame:
        path = self.__get_output_path(index)
        msg = "More than one reviewed segmentation was identified. Please ensure there's only one NIfTI file present"

        report_data = {
            "status": -self.status_code,
            "status_name": "MULTIPLE_ANNOTATIONS_ERROR",
            "comment": msg,
            "data_path": path,
            "labels_path": "",
        }
        update_row_with_dict(report, report_data, index)
        return report

    def could_run(self, index: Union[str, int], report: pd.DataFrame) -> bool:
        out_path = self.__get_output_path(index)
        cases = []
        if os.path.exists(out_path):
            cases = os.listdir(out_path)

        segmentation_exists = os.path.exists(self.__get_input_path(index))
        annotation_exists = len(cases) == 1
        return segmentation_exists and not annotation_exists

    def execute(self, index: Union[str, int], report: pd.DataFrame) -> Tuple[pd.DataFrame, bool]:
        """Manual steps are by definition not doable by an algorithm. Therefore,
        execution of this step leads to a failed stage message, indicating that
        the manual step has not been done.

        Args:
            index (Union[str, int]): _description_
            report (pd.DataFrame): _description_

        Returns:
            pd.DataFrame: _description_
        """

        # Generate a hidden copy of the baseline segmentations
        in_path = self.__get_input_path(index)
        out_path = self.__get_output_path(index)
        under_review_path = self.__get_under_review_path(index)
        bak_path = self.__get_backup_path(index)
        if not os.path.exists(bak_path):
            shutil.copytree(in_path, bak_path)
            set_files_read_only(bak_path)
        os.makedirs(under_review_path, exist_ok=True)
        os.makedirs(out_path, exist_ok=True)

        cases = os.listdir(out_path)

        if len(cases) > 1:
            # Found more than one reviewed case
            return self.__report_multiple_cases_error(index, report), False
        elif not len(cases):
            # Found no cases yet reviewed
            return self.__report_step_missing(index, report), False
        return self.__report_success(index, report), True
