from pandas import DataFrame
from typing import Union, List, Tuple
from tqdm import tqdm
import traceback
import yaml

from .dset_stage import DatasetStage
from .row_stage import RowStage
from .stage import Stage
from .utils import cleanup_storage


def normalize_report_paths(report: DataFrame) -> DataFrame:
    """Ensures paths are normalized and converts them to relative paths for the local machine

    Args:
        report (DataFrame): report to normalize

    Returns:
        DataFrame: report with transformed paths
    """
    report["data_path"] = report["data_path"].str.split("mlcube_io3").str[-1]
    report["labels_path"] = report["labels_path"].str.split("mlcube_io3").str[-1]
    return report


def write_report(report: DataFrame, filepath: str):
    report = normalize_report_paths(report)
    report_dict = report.to_dict()
    with open(filepath, "w") as f:
        yaml.dump(report_dict, f)


class Pipeline:
    def __init__(
        self,
        init_stage: DatasetStage,
        stages: List[Union[DatasetStage, RowStage]],
        staging_folders: List[str],
        trash_folders: List[str],
    ):
        self.init_stage = init_stage
        self.stages = stages
        self.staging_folders = staging_folders
        self.trash_folders = trash_folders

    def __is_subject_done(self, subject: Union[str, int], report: DataFrame) -> bool:
        """Determines if a subject is considered done

        Args:
            subject (Union[str, int]): subject index
            report (DataFrame): DaraFrame containing the state of the processing

        Returns:
            bool: wether the subject is done or not
        """
        subject_status = report.loc[subject, "status_name"]

        return subject_status == "DONE"

    def __is_done(self, report: DataFrame) -> bool:
        """Determines if the preparation is complete

        Args:
            report (DataFrame): DataFrame containing the state of the processing

        Returns:
            bool: Wether the preparation is complete
        """
        return all(report["status_name"] == "DONE")

    def __get_report_stage_to_run(
        self, subject: Union[str, int], report: DataFrame
    ) -> Union[DatasetStage, RowStage]:
        """Retrieves the stage a subject is in indicated by the report

        Args:
            subject (Union[str, int]): Subject index
            report (DataFrame): Dataframe containing the state of the processing

        Returns:
            Union[DatasetStage, RowStage]: Stage the current subject is in
        """
        report_status_code = int(report.loc[subject, "status"])
        if report_status_code < 0:
            # Error code, rerun the stage specified in the report
            report_status_code = abs(report_status_code)
        else:
            # Success code, reported stage works so move to the next one
            report_status_code += 1
        for stage in self.stages:
            if stage.status_code == report_status_code:
                return stage

        return None

    def determine_next_stage(
        self, subject: Union[str, int], report
    ) -> Tuple[List[Union[DatasetStage, RowStage]], bool]:
        could_run_stages = []
        for i, stage in enumerate(self.stages):
            could_run = False
            if isinstance(stage, RowStage):
                could_run = stage.could_run(subject, report)
            else:
                could_run = stage.could_run(report)

            if could_run:
                runnable_stage = self.stages[i]
                could_run_stages.append(runnable_stage)

        if len(could_run_stages) == 1:
            return could_run_stages[0], False

        # Handle errors
        # Either no stage can be executed (len(could_run_stages == 0))
        # or multiple stages can be executed (len(could_run_stages > 1))
        report_stage = self.__get_report_stage_to_run(subject, report)

        if len(could_run_stages) == 0:
            # Either the case processing was on-going but it's state is broken
            # or the next stage is a dataset stage, which means we're done with this one
            # or the case is done and no stage can nor should run
            # We can tell this by looking at the report
            is_done = self.__is_subject_done(subject, report)
            is_dset_stage = isinstance(report_stage, DatasetStage)
            if is_done or is_dset_stage:
                return None, True
            else:
                return None, False
        else:
            # Multiple stages could run. Remove ambiguity by
            # syncing with the report
            if report_stage in could_run_stages:
                return report_stage, False

            return None, False

    def run(self, report: DataFrame, report_path: str):
        # cleanup the trash at the very beginning
        cleanup_storage(self.trash_folders)

        # The init stage always has to be executed
        report, _ = self.init_stage.execute(report)
        write_report(report, report_path)

        should_loop = True
        while should_loop:
            # Sine we could have row and dataset stages interwoven, we want
            # to make sure we continue processing subjects until nothing new has happened.
            # This means we can resume a given subject and its row stages even after a dataset stage
            prev_status = report["status"].copy()
            subjects = list(report.index)
            subjects_loop = tqdm(subjects)

            for subject in subjects_loop:
                report = self.process_subject(
                    subject, report, report_path, subjects_loop
                )

            # Check for report differences. If there are, rerun the loop
            should_loop = any(report["status"] != prev_status)

        if self.__is_done(report):
            cleanup_folders = self.staging_folders + self.trash_folders
            cleanup_storage(cleanup_folders)

    def process_subject(
        self, subject: Union[int, str], report: DataFrame, report_path: str, pbar: tqdm
    ):
        # TODO: implement general cleanup
        # self.cleanup(subject)
        # self.check_for_errors_that_need_user_attention(subject, report)
        while True:
            # TODO also handle when everything's done
            stage, done = self.determine_next_stage(subject, report)

            if done:
                break

            try:
                report, successful = self.run_stage(
                    stage, subject, report, report_path, pbar
                )
            except Exception:
                report = self.__report_unhandled_exception(
                    stage, subject, report, report_path
                )
                print(traceback.format_exc())
                successful = False

            if not successful:
                print("breaking")
                break

        return report

    def run_stage(self, stage, subject, report, report_path, pbar):
        successful = False
        if isinstance(stage, RowStage):
            pbar.set_description(f"{subject} | {stage.name}")
            report, successful = stage.execute(subject, report)
        elif isinstance(stage, DatasetStage):
            pbar.set_description(f"{stage.name}")
            report, successful = stage.execute(report)

        write_report(report, report_path)

        return report, successful

    def __report_unhandled_exception(
        self,
        stage: Stage,
        subject: Union[int, str],
        report: DataFrame,
        report_path: str,
    ):
        # Assign a special status code for unhandled errors, associated
        # to the stage status code
        status_code = -stage.status_code - 0.101
        name = f"{stage.name.upper().replace(' ', '_')}_UNHANDLED_ERROR"
        comment = traceback.format_exc()
        data_path = report.loc[subject]["data_path"]
        labels_path = report.loc[subject]["labels_path"]

        body = {
            "status": status_code,
            "status_name": name,
            "comment": comment,
            "data_path": data_path,
            "labels_path": labels_path,
        }

        report.loc[subject] = body

        write_report(report, report_path)

        return report
