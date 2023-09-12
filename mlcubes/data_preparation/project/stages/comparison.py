from typing import Union, Tuple
import os
import shutil

import pandas as pd
from pandas import DataFrame
import numpy as np
import nibabel as nib

from .row_stage import RowStage
from .utils import get_id_tp, update_row_with_dict
from .constants import TUMOR_MASK_FOLDER, INTERIM_FOLDER


class SegmentationComparisonStage(RowStage):
    def __init__(
        self,
        data_csv: str,
        out_path: str,
        prev_stage_path,
        backup_path: str,
    ):
        self.data_csv = data_csv
        self.out_path = out_path
        self.prev_stage_path = prev_stage_path
        self.backup_path = backup_path
        self.status_code = 6

    def get_name(self):
        return "Label Segmentation Comparison"

    def __get_input_path(self, index: Union[str, int]) -> str:
        id, tp = get_id_tp(index)
        path = os.path.join(self.prev_stage_path, INTERIM_FOLDER, id, tp)
        return path

    def __get_backup_path(self, index: Union[str, int]) -> str:
        id, tp = get_id_tp(index)
        path = os.path.join(self.backup_path, id, tp, TUMOR_MASK_FOLDER)
        return path

    def __get_output_path(self, index: Union[str, int]) -> str:
        id, tp = get_id_tp(index)
        path = os.path.join(self.out_path, id, tp)
        return path

    def __get_case_path(self, index: Union[str, int]) -> str:
        path = os.path.join(self.__get_input_path(index), "reviewed")
        case = os.listdir(path)[0]

        return os.path.join(path, case)

    def __report_gt_not_found(
        self, index: Union[str, int], report: pd.DataFrame
    ) -> pd.DataFrame:
        case_path = self.__get_case_path(index)
        data_path = report.loc[index, "data_path"]
        msg = (
            f"The original segmentation for file {case_path} was not identified. "
            + "This most probably means the annotation file was renamed. "
            + "Please ensure the reviewed file retains its original name."
        )

        report_data = {
            "status": -self.status_code,
            "status_name": "ANNOTATION_COMPARISON_FAILED",
            "comment": msg,
            "data_path": data_path,
            "labels_path": case_path,
        }
        update_row_with_dict(report, report_data, index)
        return report

    def __report_exact_match(
        self, index: Union[str, int], report: pd.DataFrame
    ) -> pd.DataFrame:
        case_path = self.__get_case_path(index)
        data_path = report.loc[index, "data_path"]
        msg = (
            "The annotation seems to be a copy of a baseline segmentation, "
            + "Please ensure this is intended. Waiting on all cases to be "
            + "reviewed in order to proceed."
        )
        report_data = {
            "status": self.status_code,
            "status_name": "EXACT_MATCH_IDENTIFIED",
            "comment": msg,
            "data_path": data_path,
            "labels_path": case_path,
            "num_changed_voxels": 0,
        }
        update_row_with_dict(report, report_data, index)
        return report

    def __report_success(
        self, index: Union[str, int], report: pd.DataFrame, num_changed_voxels: int
    ) -> pd.DataFrame:
        case_path = self.__get_case_path(index)
        data_path = report.loc[index, "data_path"]
        report_data = {
            "status": self.status_code,
            "status_name": "COMPARISON_COMPLETED",
            "comment": "",
            "data_path": data_path,
            "labels_path": case_path,
            "num_changed_voxels": num_changed_voxels,
        }
        update_row_with_dict(report, report_data, index)
        return report

    def could_run(self, index: Union[str, int], report: DataFrame) -> bool:
        # Ensure a single reviewed segmentation file exists
        path = os.path.join(self.__get_input_path(index), "reviewed")
        gt_path = self.__get_backup_path(index)

        is_valid = True
        path_exists = os.path.exists(path)
        gt_path_exists = os.path.exists(gt_path)
        contains_case = False
        if path_exists:
            num_cases = len(os.listdir(path))
            contains_case = num_cases == 1
        is_valid = path_exists and contains_case and gt_path_exists

        return is_valid

    def execute(
        self, index: Union[str, int], report: DataFrame
    ) -> Tuple[DataFrame, bool]:
        path = os.path.join(self.__get_input_path(index), "reviewed")
        cases = os.listdir(path)

        match_output_path = self.__get_output_path(index)
        os.makedirs(match_output_path, exist_ok=True)
        # Get the necessary files for match check
        # We assume reviewed and gt files have the same name
        reviewed_file = os.path.join(path, cases[0])
        gt_file = os.path.join(self.__get_backup_path(index), cases[0])

        if not os.path.exists(gt_file):
            # Ground truth file not found, reviewed file most probably renamed
            report = self.__report_gt_not_found(index, report, reviewed_file)
            return report, False

        reviewed_img = nib.load(reviewed_file)
        gt_img = nib.load(gt_file)

        reviewed_voxels = np.array(reviewed_img.dataobj)
        gt_voxels = np.array(gt_img.dataobj)

        num_changed_voxels = np.sum(reviewed_voxels != gt_voxels)

        # At this point we know we succeeded, remove the backup files
        shutil.rmtree(self.__get_backup_path(index))

        if num_changed_voxels == 0:
            report = self.__report_exact_match(index, report)
            return report, True

        report = self.__report_success(index, report, num_changed_voxels)
        return report, True
