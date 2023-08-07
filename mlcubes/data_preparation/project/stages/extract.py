from typing import Union, List
from tqdm import tqdm
import pandas as pd
import os
import shutil
import traceback

from .row_stage import RowStage
from .PrepareDataset import Preparator, INTERIM_FOLDER, FINAL_FOLDER
from .utils import update_row_with_dict, get_id_tp, MockTqdm


class Extract(RowStage):
    def __init__(
        self,
        data_csv: str,
        out_path: str,
        subpaths: List[str],
        prev_stage_path: str,
        prev_subpaths: List[str],
        pbar: tqdm,
        func_name: str,
        status_code: int,
    ):
        self.data_csv = data_csv
        self.out_path = out_path
        self.subpaths = subpaths
        self.prev_path = prev_stage_path
        self.prev_subpaths = prev_subpaths
        os.makedirs(self.out_path, exist_ok=True)
        self.prep = Preparator(data_csv, out_path, "BraTSPipeline")
        self.func_name = func_name
        self.func = getattr(self.prep, func_name)
        self.pbar = pbar
        self.failed = False
        self.exception = None
        self.status_code = status_code

    def get_name(self) -> str:
        return self.func_name.replace("_", " ").capitalize()

    def should_run(self, index: Union[str, int], report: pd.DataFrame) -> bool:
        """Determine if case at given index needs to be converted to NIfTI

        Args:
            index (Union[str, int]): Case index, as used by the report dataframe
            report (pd.DataFrame): Report Dataframe for providing additional context

        Returns:
            bool: Wether this stage should be executed for the given case
        """
        prev_paths = self.__get_paths(index, self.prev_path, self.prev_subpaths)
        return all([os.path.exists(path) for path in prev_paths])

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
        report = self.__update_state(index, report)
        self.prep.write()

        return report

    def __prepare_exec(self):
        # Reset the file contents for errors
        open(self.prep.stderr_log, "w").close()

        # Update the out dataframes to current state
        self.prep.read()

    def __get_paths(self, index: Union[str, int], path: str, subpaths: List[str]):
        id, tp = get_id_tp(index)
        out_paths = []
        for subpath in subpaths:
            out_paths.append(os.path.join(path, subpath, id, tp))
        return out_paths

    def __copy_case(self, index: Union[str, int]):
        prev_paths = self.__get_paths(index, self.prev_path, self.prev_subpaths)
        copy_paths = self.__get_paths(index, self.out_path, self.prev_subpaths)
        for prev, copy in zip(prev_paths, copy_paths):
            shutil.copytree(prev, copy, dirs_exist_ok=True)

    def __process_case(self, index: Union[str, int]):
        id, tp = get_id_tp(index)
        df = self.prep.subjects_df
        row = df[(df["SubjectID"] == id) & (df["Timepoint"] == tp)].iloc[0]
        try:
            self.func(row, self.pbar)
        except Exception as e:
            self.failed = True
            self.exception = e
            self.traceback = traceback.format_exc()

    def __update_state(
        self, index: Union[str, int], report: pd.DataFrame
    ) -> pd.DataFrame:
        if self.failed:
            del_paths = self.__get_paths(index, self.out_path, self.subpaths)
            report = self.__report_failure(index, report)
        else:
            del_paths = self.__get_paths(index, self.prev_path, self.prev_subpaths)
            report = self.__report_success(index, report)

        for path in del_paths:
            shutil.rmtree(path, ignore_errors=True)

        return report

    def __report_success(
        self, index: Union[str, int], report: pd.DataFrame
    ) -> pd.DataFrame:
        paths = self.__get_paths(index, self.out_path, self.subpaths)
        report_data = {
            "status": self.status_code,
            "status_name": f"{self.func_name.upper()}_FINISHED",
            "comment": "",
            "data_path": ",".join(paths),
            "labels_path": "",
        }
        update_row_with_dict(report, report_data, index)
        return report

    def __report_failure(
        self, index: Union[str, int], report: pd.DataFrame
    ) -> pd.DataFrame:
        prev_paths = self.__get_paths(index, self.prev_path, self.prev_subpaths)
        msg = f"{str(self.exception)}: {self.traceback}"

        report_data = {
            "status": -self.status_code,
            "status_name": f"{self.func_name.upper()}_FAILED",
            "comment": msg,
            "data_path": ",".join(prev_paths),
            "labels_path": "",
        }
        update_row_with_dict(report, report_data, index)
        return report
