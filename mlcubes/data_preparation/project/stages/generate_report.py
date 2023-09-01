from .dset_stage import DatasetStage
import pandas as pd
import numpy as np
import os
from pathlib import Path
import hashlib
import shutil


# Taken from https://stackoverflow.com/questions/24937495/how-can-i-calculate-a-hash-for-a-filesystem-directory-using-python
def md5_update_from_dir(directory, hash):
    assert Path(directory).is_dir()
    for path in sorted(Path(directory).iterdir(), key=lambda p: str(p).lower()):
        hash.update(path.name.encode())
        if path.is_file():
            with open(path, "rb") as f:
                for chunk in iter(lambda: f.read(4096), b""):
                    hash.update(chunk)
        elif path.is_dir():
            hash = md5_update_from_dir(path, hash)
    return hash


def md5_dir(directory):
    return md5_update_from_dir(directory, hashlib.md5()).hexdigest()


class GenerateReport(DatasetStage):
    def __init__(self, input_path: str, output_path: str):
        self.input_path = input_path
        self.output_path = output_path

    def should_run(self, report: pd.DataFrame):
        return True

    def execute(self, report: pd.DataFrame):
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

        observed_cases = set()

        for subject in os.listdir(self.input_path):
            in_subject_path = os.path.join(self.input_path, subject)
            out_subject_path = os.path.join(self.output_path, subject)

            if not os.path.isdir(in_subject_path):
                continue

            for timepoint in os.listdir(in_subject_path):
                in_tp_path = os.path.join(in_subject_path, timepoint)
                out_tp_path = os.path.join(out_subject_path, timepoint)

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

                data = {
                    "status": 0,
                    "status_name": "IDENTIFIED",
                    "comment": "",
                    "data_path": out_tp_path,
                    "labels_path": "",
                    "num_changed_voxels": np.nan,
                    "input_hash": input_hash,
                }
                subject_series = pd.Series(data)
                subject_series.name = index
                report = report.append(subject_series)

        reported_cases = set(report.index)
        removed_cases = reported_cases - observed_cases

        # Stop reporting removed cases
        for case_index in removed_cases:
            report = report.drop(case_index)

        return report
