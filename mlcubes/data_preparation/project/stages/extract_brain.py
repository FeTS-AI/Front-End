from typing import Union
from tqdm import tqdm
import pandas as pd
import os
import shutil

from .row_stage import RowStage
from .PrepareDataset import Preparator, INTERIM_FOLDER, FINAL_FOLDER
from .utils import update_row_with_dict, get_id_tp, MockTqdm


class ExtractBrain(RowStage):
    def __init__(self, data_csv: str, out_path: str, prev_stage_path: str, pbar: tqdm):
        self.data_csv = data_csv
        self.out_path = out_path
        self.prev_stage_path = prev_stage_path
        os.makedirs(self.out_path, exist_ok=True)
        self.prep = Preparator(data_csv, out_path, "BraTSPipeline")
        self.pbar = pbar
        self.failed = False
        self.exception = None

    def get_name(self) -> str:
        return "NiFTI Conversion"

    def should_run(self, index: Union[str, int], report: pd.DataFrame) -> bool:
        """Determine if case at given index needs to be converted to NIfTI

        Args:
            index (Union[str, int]): Case index, as used by the report dataframe
            report (pd.DataFrame): Report Dataframe for providing additional context

        Returns:
            bool: Wether this stage should be executed for the given case
        """
        prev_fets_path, prev_qc_path = self.__get_prev_output_paths(index)
        return os.path.exists(prev_fets_path) and os.path.exists(prev_qc_path)

    def execute(self, index: Union[str, int], report: pd.DataFrame) -> pd.DataFrame:
        """Executes the NIfTI transformation stage on the given case

        Args:
            index (Union[str, int]): case index, as used by the report
            report (pd.DataFrame): DataFrame containing the current state of the preparation flow

        Returns:
            pd.DataFrame: Updated report dataframe
        """
        self.__prepare_exec()
        self.__copy_case(index)
        self.__process_case(index)
        report = self.__update_report(index, report)
        self.prep.write()

        return report

    def __prepare_exec(self):
        # Reset the file contents for errors
        open(self.prep.stderr_log, "w").close()

        # Update the out dataframes to current state
        self.prep.read()

    def __get_prev_output_paths(self, index: Union[str, int]):
        id, tp = get_id_tp(index)
        prev_fets_path = os.path.join(self.prev_stage_path, FINAL_FOLDER, id, tp)
        prev_qc_path = os.path.join(self.prev_stage_path, INTERIM_FOLDER, id, tp)
        return prev_fets_path, prev_qc_path

    def __get_output_paths(self, index: Union[str, int]):
        id, tp = get_id_tp(index)
        fets_path = os.path.join(self.prep.final_output_dir, id, tp)
        qc_path = os.path.join(self.prep.interim_output_dir, id, tp)

        return fets_path, qc_path

    def __copy_case(self, index: Union[str, int]):
        prev_fets_path, prev_qc_path = self.__get_prev_output_paths(index)
        fets_path, qc_path = self.__get_output_paths(index)
        shutil.copytree(prev_fets_path, fets_path, dirs_exist_ok=True)
        shutil.copytree(prev_qc_path, qc_path, dirs_exist_ok=True)

    def __process_case(self, index: Union[str, int]):
        id, tp = get_id_tp(index)
        df = self.prep.subjects_df
        row = df[(df["SubjectID"] == id) & (df["Timepoint"] == tp)].iloc[0]
        try:
            self.prep.extract_brain(row, self.pbar)
        except Exception as e:
            self.failed = True
            self.exception = e

    def __update_prev_stage_state(self, index: Union[str, int], report: pd.DataFrame):
        prev_data_path = report.loc[index]["data_path"]
        shutil.rmtree(prev_data_path)

    def __undo_current_stage_changes(self, index: Union[str, int]):
        id, tp = get_id_tp(index)
        fets_path = os.path.join(self.out_path, FINAL_FOLDER, id, tp)
        qc_path = os.path.join(self.out_path, INTERIM_FOLDER, id, tp)
        shutil.rmtree(fets_path, ignore_errors=True)
        shutil.rmtree(qc_path, ignore_errors=True)

    def __update_report(
        self, index: Union[str, int], report: pd.DataFrame
    ) -> pd.DataFrame:
        if self.failed:
            self.__undo_current_stage_changes(index)
            report = self.__report_failure(index, report)
        else:
            self.__update_prev_stage_state(index, report)
            report = self.__report_success(index, report)

        return report

    def __report_success(
        self, index: Union[str, int], report: pd.DataFrame
    ) -> pd.DataFrame:
        paths = self.__get_output_paths(index)
        report_data = {
            "status": 3,
            "status_name": "BRAIN_EXTRACTED",
            "comment": "",
            "data_path": ",".join(paths),
            "labels_path": "",
        }
        update_row_with_dict(report, report_data, index)
        return report

    def __report_failure(
        self, index: Union[str, int], report: pd.DataFrame
    ) -> pd.DataFrame:
        prev_paths = self.__get_prev_output_paths(index)
        msg = str(self.exception)

        report_data = {
            "status": -3,
            "status_name": "BRAIN_EXTRACTION_FAILED",
            "comment": msg,
            "data_path": ",".join(prev_paths),
            "labels_path": "",
        }
        update_row_with_dict(report, report_data, index)
        return report
