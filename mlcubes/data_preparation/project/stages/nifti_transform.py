from typing import Union
from tqdm import tqdm
import pandas as pd
import os
import shutil

from .row_stage import RowStage
from .PrepareDataset import Preparator, INTERIM_FOLDER, FINAL_FOLDER
from .utils import update_row_with_dict, get_id_tp, MockTqdm, unnormalize_path


class NIfTITransform(RowStage):
    def __init__(self, data_csv: str, out_path: str, prev_stage_path: str, pbar: tqdm):
        self.data_csv = data_csv
        self.out_path = out_path
        self.prev_stage_path = prev_stage_path
        os.makedirs(self.out_path, exist_ok=True)
        self.prep = Preparator(data_csv, out_path, "BraTSPipeline")
        # self.pbar = pbar
        self.pbar = tqdm()
        self.status_code = 2

    def get_name(self) -> str:
        return "NiFTI Conversion"

    def could_run(self, index: Union[str, int], report: pd.DataFrame) -> bool:
        """Determine if case at given index needs to be converted to NIfTI

        Args:
            index (Union[str, int]): Case index, as used by the report dataframe
            report (pd.DataFrame): Report Dataframe for providing additional context

        Returns:
            bool: Wether this stage could be executed for the given case
        """
        id, tp = get_id_tp(index)
        prev_case_path = os.path.join(self.prev_stage_path, id, tp)
        if os.path.exists(prev_case_path):
            return len(os.listdir(prev_case_path)) > 0
        return False

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
        report, success = self.__update_report(index, report)
        self.prep.write()

        return report, success

    def __get_output_paths(self, index: Union[str, int]):
        id, tp = get_id_tp(index)
        fets_path = os.path.join(self.prep.final_output_dir, id, tp)
        qc_path = os.path.join(self.prep.interim_output_dir, id, tp)

        return fets_path, qc_path

    def __prepare_exec(self):
        # Reset the file contents for errors
        open(self.prep.stderr_log, "w").close()

        self.prep.read()

    def __process_case(self, index: Union[str, int]):
        id, tp = get_id_tp(index)
        df = self.prep.subjects_df
        row = df[(df["SubjectID"] == id) & (df["Timepoint"] == tp)].iloc[0]
        self.prep.convert_to_dicom(index, row, self.pbar)

    def __update_prev_stage_state(self, index: Union[str, int], report: pd.DataFrame):
        prev_data_path = report.loc[index]["data_path"]
        prev_data_path = unnormalize_path(prev_data_path, "mlcube_io3")
        shutil.rmtree(prev_data_path)

    def __undo_current_stage_changes(self, index: Union[str, int]):
        fets_path, qc_path = self.__get_output_paths(index)
        shutil.rmtree(fets_path, ignore_errors=True)
        shutil.rmtree(qc_path, ignore_errors=True)

    def __update_report(
        self, index: Union[str, int], report: pd.DataFrame
    ) -> pd.DataFrame:
        # TODO: What could be reported? We have the processed data,
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
            success = False
        else:
            self.__update_prev_stage_state(index, report)
            report = self.__report_success(index, report)
            success = True

        return report, success

    def __report_success(
        self, index: Union[str, int], report: pd.DataFrame
    ) -> pd.DataFrame:
        fets_path, qc_path = self.__get_output_paths(index)
        report_data = {
            "status": self.status_code,
            "status_name": "CONVERTED_TO_NIfTI",
            "comment": "",
            "data_path": qc_path,
            "labels_path": fets_path,
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
            "status": -self.status_code,
            "status_name": "NIfTI_CONVERSION_FAILED",
            "comment": msg,
            "data_path": prev_data_path,
            "labels_path": "",
        }
        update_row_with_dict(report, report_data, index)
        return report
