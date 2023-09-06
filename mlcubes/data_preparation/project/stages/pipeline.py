from pandas import DataFrame
from typing import Union, List, Tuple
from tqdm import tqdm
import yaml

from .dset_stage import DatasetStage
from .row_stage import RowStage

def write_report(report: DataFrame, filepath: str):
    report_dict = report.to_dict()
    with open(filepath, "w") as f:
        yaml.dump(report_dict, f)


class Pipeline:
    def __init__(self, init_stage: DatasetStage, stages: List[Union[DatasetStage, RowStage]]):
        self.init_stage = init_stage
        self.stages = stages

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

    def __get_report_stage_to_run(self, subject: Union[str, int], report: DataFrame) -> Union[DatasetStage, RowStage]:
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

        if len(could_run_stages) == 0:
            # Either the case processing was on-going but it's state is broken
            # or the case is done and no stage can nor should run
            # We can tell this by looking at the report
            if self.__is_subject_done(subject, report):
                return None, True
            else:
                return None, False
        else:
            # Multiple stages could run. Remove ambiguity by
            # syncing with the report
            report_stage = self.__get_report_stage_to_run(subject, report)
            if report_stage in could_run_stages:
                return report_stage, False
            
            return None, False

    def run(self, report: DataFrame, report_path: str):
        # The init stage always has to be executed
        report, _ = self.init_stage.execute(report)
        write_report(report, report_path)

        subjects = list(report.index)
        subjects_loop = tqdm(subjects)

        for subject in subjects_loop:
            report = self.process_subject(subject, report, report_path, subjects_loop)

    def process_subject(self, subject: Union[int, str], report: DataFrame, report_path: str, pbar: tqdm):
        # TODO: implement general cleanup
        # self.cleanup(subject)
        # self.check_for_errors_that_need_user_attention(subject, report)
        while True:
            # TODO also handle when everything's done
            stage, done = self.determine_next_stage(subject, report)

            if done:
                break

            if isinstance(stage, RowStage):
                pbar.set_description(f"{subject} | {stage.get_name()}")
                report, successful = stage.execute(subject, report)
            else:
                pbar.set_description(f"{stage.get_name()}")
                report, successful = stage.execute(report)

            write_report(report, report_path)

            if not successful:
                break


        return report
