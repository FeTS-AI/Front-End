from .dset_stage import DatasetStage
import pandas as pd
import numpy as np
import os
import shutil
from typing import Tuple
from .utils import has_prepared_folder_structure, md5_dir


class GenerateReport(DatasetStage):
    def __init__(
        self,
        input_path: str,
        output_path: str,
        input_labels_path: str,
        output_labels_path,
        done_data_out_path: str,
        done_status: int,
    ):
        self.input_path = input_path
        self.output_path = output_path
        self.input_labels_path = input_labels_path
        self.output_labels_path = output_labels_path
        self.done_data_out_path = done_data_out_path
        self.done_status_code = done_status

    @property
    def name(self) -> str:
        return "Generate Report"

    @property
    def status_code(self) -> int:
        return 0

    def could_run(self, report: pd.DataFrame):
        return True

    def execute(self, report: pd.DataFrame) -> Tuple[pd.DataFrame, bool]:
        # Rewrite the report
        cols = [
            "status",
            "status_name",
            "comment",
            "data_path",
            "labels_path",
            "input_hash",
        ]
        if report is None:
            report = pd.DataFrame(columns=cols)

        input_is_prepared = has_prepared_folder_structure(
            self.input_path, self.input_labels_path
        )
        if input_is_prepared:
            # If prepared, store data directly in the data folder
            self.output_path = self.done_data_out_path

        observed_cases = set()

        for subject in os.listdir(self.input_path):
            in_subject_path = os.path.join(self.input_path, subject)
            out_subject_path = os.path.join(self.output_path, subject)
            in_labels_subject_path = os.path.join(self.input_labels_path, subject)
            out_labels_subject_path = os.path.join(self.output_labels_path, subject)

            if not os.path.isdir(in_subject_path):
                continue

            for timepoint in os.listdir(in_subject_path):
                in_tp_path = os.path.join(in_subject_path, timepoint)
                out_tp_path = os.path.join(out_subject_path, timepoint)
                in_labels_tp_path = os.path.join(in_labels_subject_path, timepoint)
                out_labels_tp_path = os.path.join(out_labels_subject_path, timepoint)

                if not os.path.isdir(in_tp_path):
                    continue

                input_hash = md5_dir(in_tp_path)

                index = f"{subject}|{timepoint}"

                # Keep track of the cases that were found on the input folder
                observed_cases.add(index)

                if index in report.index:
                    # Case has already been identified, see if input hash is different
                    # if so, override the contents and restart the state for that case
                    if report.loc[index]["input_hash"] == input_hash:
                        continue

                    shutil.rmtree(out_tp_path, ignore_errors=True)
                    shutil.copytree(in_tp_path, out_tp_path)
                    report = report.drop(index)
                else:
                    # New case not identified by the report. Add it
                    shutil.copytree(in_tp_path, out_tp_path)

                data = {
                    "status": self.status_code,
                    "data_path": out_tp_path,
                    "labels_path": "",
                    "num_changed_voxels": np.nan,
                    "brain_mask_hash": "",
                    "input_hash": input_hash,
                }
                if input_is_prepared:
                    data["status_name"] = "DONE"
                    data["status_code"] = self.done_status_code
                    shutil.rmtree(out_labels_tp_path, ignore_errors=True)
                    shutil.copytree(in_labels_tp_path, out_labels_tp_path)

                subject_series = pd.Series(data)
                subject_series.name = index
                report = report.append(subject_series)

        reported_cases = set(report.index)
        removed_cases = reported_cases - observed_cases

        # Stop reporting removed cases
        for case_index in removed_cases:
            report = report.drop(case_index)

        return report, True
