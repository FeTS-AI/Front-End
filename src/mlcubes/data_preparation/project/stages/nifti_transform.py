from typing import Union
import pandas as pd
import os
import shutil

from .row_stage import RowStage
from .PrepareDataset import Preparator
from .utils import update_row_with_dict, get_id_tp


class NIfTITransform(RowStage):
    def __init__(self, data_csv: str, out_path: str, prev_stage_path: str):
        self.data_csv = data_csv
        self.out_path = out_path
        self.prev_stage_path = prev_stage_path
        os.makedirs(self.out_path, exist_ok=True)
        self.prep = Preparator(data_csv, out_path, "BraTSPipeline")

    def should_run(self, index: Union[str, int], report: pd.DataFrame) -> bool:
        """Determine if case at given index needs to be converted to NIfTI

        Args:
            index (Union[str, int]): Case index, as used by the report dataframe
            report (pd.DataFrame): Report Dataframe for providing additional context

        Returns:
            bool: Wether this stage should be executed for the given case
        """
        id, tp = get_id_tp(index)
        prev_case_path = os.path.join(self.prev_stage_path, id, tp)
        return os.path.exists(prev_case_path)

    def execute(self, index: Union[str, int], report: pd.DataFrame) -> pd.DataFrame:
        """Executes the NIfTI transformation stage on the given case

        Args:
            index (Union[str, int]): case index, as used by the report
            report (pd.DataFrame): DataFrame containing the current state of the preparation flow

        Returns:
            pd.DataFrame: Updated report dataframe
        """
        self.__prepare_exec()
        self.__process_case(index)
        report = self.__update_report(index, report)
        self.prep.write()

        return report

    def __prepare_exec(self):
        # Reset the file contents for errors
        open(self.prep.stderr_log, "w").close()

        # Update the out dataframes to current state
        self.prep.read()

    def __process_case(self, index: Union[str, int]):
        id, tp = get_id_tp(index)
        df = self.prep.subjects_df
        row = df[(df["SubjectID"] == id) & (df["Timepoint"] == tp)].iloc[0]
        self.prep.process_row(index, row)

    def __update_prev_stage_state(self, index: Union[str, int], report: pd.DataFrame):
        prev_data_path = report.loc[index]["data_path"]
        shutil.rmtree(prev_data_path)

    def __undo_current_stage_changes(self, index: Union[str, int]):
        id, tp = get_id_tp(index)
        fets_path = os.path.join(self.out_path, "DataForFeTS", id, tp)
        qc_path = os.path.join(self.out_path, "DataForQC", id, tp)
        shutil.rmtree(fets_path, ignore_errors=True)
        shutil.rmtree(qc_path, ignore_errors=True)

    def __update_report(
        self, index: Union[str, int], report: pd.DataFrame
    ) -> pd.DataFrame:
        # TODO: What should be reported? We have the processed data,
        # QC subjects with negative intensities and
        # QC subjects with bratspipeline error
        id, tp = get_id_tp(index)
        failing = self.prep.failing_subjects
        failing_subject = failing[
            (failing["SubjectID"] == id) & (failing["Timepoint"] == tp)
        ]
        if len(failing_subject):
            self.__undo_current_stage_changes(index)
            report = self.__report_failure(index, report)
        else:
            self.__update_prev_stage_state(index, report)
            report = self.__report_success(index, report)

        return report

    def __report_success(
        self, index: Union[str, int], report: pd.DataFrame
    ) -> pd.DataFrame:
        id, tp = get_id_tp(index)
        fets_path = os.path.join(self.out_path, "DataForFeTS", id, tp)
        report_data = {
            "status": 2,
            "status_name": "CONVERTED_TO_NIfTI",
            "comment": "",
            "data_path": fets_path,
            "labels_path": "",
        }
        update_row_with_dict(report, report_data, index)
        return report

    def __report_failure(
        self, index: Union[str, int], report: pd.DataFrame
    ) -> pd.DataFrame:
        id, tp = get_id_tp(index)
        prev_data_path = report.loc[index]["data_path"]

        with open(self.prep.stderr_log, "r") as f:
            msg = f.read()

        report_data = {
            "status": -2,
            "status_name": "NIfTI_CONVERSION_FAILED",
            "comment": msg,
            "data_path": prev_data_path,
            "labels_path": "",
        }
        update_row_with_dict(report, report_data, index)
        return report
