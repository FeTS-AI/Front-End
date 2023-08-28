from typing import Union
import os
import yaml

import pandas as pd
from pandas import DataFrame

from .dset_stage import DatasetStage
from .utils import get_id_tp
from .constants import TUMOR_MASK_FOLDER, INTERIM_FOLDER
import GANDLF
from GANDLF.cli import generate_metrics


class MatchStage(DatasetStage):
    def __init__(self, data_csv: str, out_path: str, prev_stage_path, backup_path: str):
        self.data_csv = data_csv
        self.out_path = out_path
        self.prev_stage_path = prev_stage_path
        self.backup_path = backup_path

    def get_name(self):
        return "Exact Match Identification"

    def __get_input_path(self, index: Union[str, int]):
        id, tp = get_id_tp(index)
        path = os.path.join(self.prev_stage_path, INTERIM_FOLDER, id, tp)
        return path

    def __get_backup_path(self, index: Union[str, int]):
        id, tp = get_id_tp(index)
        path = os.path.join(self.backup_path, id, tp, TUMOR_MASK_FOLDER)
        return path

    def __get_output_path(self, index: Union[str, int]):
        id, tp = get_id_tp(index)
        path = os.path.join(self.out_path, id, tp)
        return path

    def should_run(self, report: DataFrame) -> bool:
        # Should run once all cases have a selected and corrected segmentation
        data = pd.read_csv(self.data_csv)

        def get_row_path(row: pd.Series):
            return os.path.join(
                self.prev_stage_path,
                INTERIM_FOLDER,
                row["SubjectID"],
                row["Timepoint"],
                "reviewed",
            )

        def valid_case(path):
            path_exists = os.path.exists(path)
            contains_case = False
            if path_exists:
                contains_case = len(os.listdir(path)) == 1
            return path_exists and contains_case

        data["expected_path"] = data.apply(get_row_path, axis="columns")
        valid_data = data["expected_path"].apply(valid_case)
        return all(valid_data)

    def execute(self, index: Union[str, int], report: DataFrame) -> DataFrame:
        # TODO: Create a CSV with the expected structure
        # TODO: Run the generate metrics command
        # TODO: Read the generated metrics and extract the necessary ones
        # TODO: Check exact matches by hash and metrics
        # TODO: Are we checking exact match against all the models?
        # TODO: Add the percent of unchanged files, as well as voxel changes
        # To the report, as separate columns

        match_output_path = self.__get_output_path(index)
        os.makedirs(match_output_path, exist_ok=True)
        # Get the necessary files for match check
        id, tp = get_id_tp(index)
        reviewed_filename = f"reviewed/{id}_{tp}_final_seg.nii.gz"
        reviewed_file = os.path.join(self.__get_input_path(index), reviewed_filename)
        gt_filename = f"{id}_{tp}_tumorMask_fused-voting.nii.gz"
        # TODO: How do we know which segmentation to compare against?
        # Should we compare against all segmentations?
        # If there's no exact match, which segmentation should we compare metrics with?
        ground_truth = os.path.join(self.__get_backup_path(index), gt_filename)

        # Prepare the assets for metrics generation
        inputdata_file = os.path.join(match_output_path, "inputdata.csv")
        data = {"subjectid": f"{id}_{tp}", "prediction": reviewed_file, "target": ground_truth}
        pd.DataFrame(data, index=[0]).to_csv(inputdata_file, index=False)

        # Read gandlf config file.
        # TODO: what are the requirements of config?
        # TODO: do NOT hardcode the filesystem names used below
        config_file = os.path.join(os.path.dirname(__file__), "data_prep_models/tumor_segmentation/model_0/config.yaml")

        out_file = os.path.join(match_output_path, "out.yaml")

        # Run the metrics generation logic
        generate_metrics.generate_metrics_dict(inputdata_file, config_file, out_file)

        # Open the generated metrics
        with open(out_file, "r") as f:
            metrics = yaml.safe_load(f)
        print(metrics)
