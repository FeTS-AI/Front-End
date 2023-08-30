from .dset_stage import DatasetStage
import pandas as pd
import numpy as np
import os


class GenerateReport(DatasetStage):
    def __init__(self, input_path: str):
        self.input_path = input_path

    def should_run(self, report: pd.DataFrame):
        return report is None

    def execute(self, report: pd.DataFrame):
        # Rewrite the report
        cols = ["status", "status_name", "comment", "data_path", "labels_path"]
        report = pd.DataFrame(columns=cols)

        for subject in os.listdir(self.input_path):
            subject_path = os.path.join(self.input_path, subject)
            if not os.path.isdir(subject_path):
                continue

            for timepoint in os.listdir(subject_path):
                tp_path = os.path.join(subject_path, timepoint)
                if not os.path.isdir(tp_path):
                    continue

                data = {
                    "status": 0,
                    "status_name": "IDENTIFIED",
                    "comment": "",
                    "data_path": tp_path,
                    "labels_path": "",
                    "changed_voxels": np.nan,
                }
                subject_series = pd.Series(data)
                subject_series.name = f"{subject}|{timepoint}"
                report = report.append(subject_series)

        return report
