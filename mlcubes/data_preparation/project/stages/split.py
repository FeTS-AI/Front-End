import os
import yaml
import pandas as pd
from typing import List
import math

from .dset_stage import DatasetStage
from .utils import get_id_tp, cleanup_storage


class SplitStage(DatasetStage):
    def __init__(
        self,
        params: str,
        data_path: str,
        labels_path: str,
        split_csv_path: str,
        staging_folders: List[str],
    ):
        self.params = params
        self.data_path = data_path
        self.labels_path = labels_path
        self.split_csv_path = split_csv_path
        self.staging_folders = staging_folders
        self.status_code = 8

    def get_name(self) -> str:
        return "Generate splits"

    def could_run(self, report: pd.DataFrame) -> bool:
        for index in report.index:
            id, tp = get_id_tp(index)
            case_data_path = os.path.join(self.data_path, id, tp)
            case_labels_path = os.path.join(self.labels_path, id, tp)
            data_exists = os.path.exists(case_data_path)
            labels_exist = os.path.exists(case_labels_path)

            if not data_exists or not labels_exist:
                return False

        return True

    def __report_success(self, report: pd.DataFrame) -> pd.DataFrame:
        report["status"] = self.status_code
        report["status_name"] = "DONE"
        report["comment"] = ""

        return report

    def execute(self, report: pd.DataFrame) -> pd.DataFrame:
        with open(self.params, "r") as f:
            params = yaml.safe_load(f)

        seed = params["seed"]
        train_pct = params["train_percent"]

        split_df = report.copy(deep=True)
        split_df["SubjectID"] = split_df.index.str.split("|").str[0]
        split_df["Timepoint"] = split_df.index.str.split("|").str[1]
        split_df = split_df[["SubjectID", "Timepoint"]].reset_index(drop=True)
        subjects = split_df["SubjectID"].drop_duplicates()
        subjects = subjects.sample(frac=1, random_state=seed)
        train_size = math.floor(len(subjects) * train_pct)

        train_subjects = subjects.iloc[:train_size]
        val_subjects = subjects.iloc[train_size:]

        train_mask = split_df["SubjectID"].isin(train_subjects)
        val_mask = split_df["SubjectID"].isin(val_subjects)

        split_df.loc[train_mask, "Split"] = "Train"
        split_df.loc[val_mask, "Split"] = "Val"

        split_df.to_csv(self.split_csv_path, index=False)

        report = self.__report_success(report)
        cleanup_storage(self.staging_folders)

        return report
