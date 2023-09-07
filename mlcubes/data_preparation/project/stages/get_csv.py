from .row_stage import RowStage
from .CreateCSVForDICOMs import CSVCreator
from .utils import update_row_with_dict, get_id_tp
from pathlib import Path

import pandas as pd
from typing import Union, Tuple
import os
import shutil


class AddToCSV(RowStage):
    def __init__(
        self, input_dir: str, output_csv: str, out_dir: str, prev_stage_path: str
    ):
        self.input_dir = input_dir
        self.output_csv = output_csv
        self.out_dir = out_dir
        self.prev_stage_path = prev_stage_path
        self.status_code = 1
        os.makedirs(self.out_dir, exist_ok=True)
        self.csv_processor = CSVCreator(self.input_dir, self.output_csv)
        if os.path.exists(self.output_csv):
            # Use the updated version of the CSV
            self.contents = pd.read_csv(self.output_csv)
            self.csv_processor.output_df_for_csv = self.contents
        else:
            # Use the default, empty version
            self.contents = self.csv_processor.output_df_for_csv

    def get_name(self) -> str:
        return "Initial Validation"

    def could_run(self, index: Union[str, int], report: pd.DataFrame) -> bool:
        """Determines if getting a new CSV is necessary.
        This is done by checking the existence of the expected file

        Args:
            index (Union[str, int]): case index in the report
            report (pd.DataFrame): Dataframe containing the current state of the preparation flow

        Returns:
            bool: wether this stage could be executed
        """
        id, tp = get_id_tp(index)
        prev_case_path = os.path.join(self.prev_stage_path, id, tp)

        return os.path.exists(prev_case_path)

    def execute(
        self, index: Union[str, int], report: pd.DataFrame
    ) -> Tuple[pd.DataFrame, bool]:
        """Adds valid cases to the data csv that is used for later processing
        Invalid cases are flagged in the report

        Args:
            index (Union[str, int]): case index in the report
            report (pd.DataFrame): DataFrame containing the current state of the preparation flow

        Returns:
            pd.DataFrame: Updated report dataframe
        """
        id, tp = get_id_tp(index)
        subject_path = os.path.join(self.input_dir, id)
        tp_path = os.path.join(subject_path, tp)
        subject_out_path = os.path.join(self.out_dir, id)
        tp_out_path = os.path.join(subject_out_path, tp)
        # We will first copy the timepoint to the out folder
        # This is so, if successful, the csv will point to the data
        # in the next stage, instead of the previous
        shutil.copytree(tp_path, tp_out_path)

        self.csv_processor.process_timepoint(tp, id, subject_out_path)
        report_data = {
            "status": self.status_code,
            "status_name": "VALIDATED",
            "comment": "",
            "data_path": tp_out_path,
            "labels_path": "",
        }
        if f"{id}_{tp}" in self.csv_processor.subject_timepoint_missing_modalities:
            shutil.rmtree(tp_out_path, ignore_errors=True)
            # Differentiate errors by floating point value
            status_code = -self.status_code - 0.1 # -1.1
            comment = "There are missing modalities. Please check the data"
            report_data["status"] = status_code 
            report_data["status_name"] = "MISSING_MODALITIES"
            report_data["data_path"] = tp_path
            report_data["comment"] = comment
            success = False
        elif f"{id}_{tp}" in self.csv_processor.subject_timepoint_extra_modalities:
            shutil.rmtree(tp_out_path, ignore_errors=True)
            # Differentiate errors by floating point value
            status_code = -self.status_code - 0.2 # -1.2
            comment = "There are extra modalities. Please check the data"
            report_data["status"] = status_code
            report_data["status_name"] = "EXTRA MODALITIES"
            report_data["data_path"] = tp_path
            report_data["comment"] = comment
            success = False
        else:
            shutil.rmtree(tp_path)
            success = True

        update_row_with_dict(report, report_data, index)

        self.csv_processor.write()

        return report, success
