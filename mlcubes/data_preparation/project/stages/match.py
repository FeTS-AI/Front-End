from typing import Union
import os
import yaml
import json

import pandas as pd
from pandas import DataFrame

from .dset_stage import DatasetStage
from .utils import get_id_tp
from .constants import TUMOR_MASK_FOLDER
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
        path = os.path.join(self.prev_stage_path, id, tp)
        return path

    def __get_backup_path(self, index: str | int):
        id, tp = get_id_tp(index)
        path = os.path.join(self.backup_path, id, tp)
        return path

    def __get_output_path(self, index: str | int):
        id, tp = get_id_tp(index)
        path = os.path.join(self.out_path, id, tp)
        return path

    def should_run(self, report: DataFrame) -> bool:
        data = pd.read_csv(self.data_csv)
        print(data)
        return False

    def execute(self, index: str | int, report: DataFrame) -> DataFrame:
        # TODO: Create a CSV with the expected structure
        # TODO: Run the generate metrics command
        # TODO: Read the generated metrics and extract the necessary ones
        # TODO: Check exact matches by hash and metrics
        # TODO: Are we checking exact match against all the models?
        # TODO: Add the percent of unchanged files, as well as voxel changes
        # To the report, as separate columns

        # Get the necessary files for match check
        id, tp = get_id_tp(index)
        reviewed_filename = f"{id}_{tp}_final_seg.nii.gz"
        reviewed_file = os.path.join(self.__get_input_path(index), reviewed_filename)
        gt_filename = ""  # TODO: How do we know which segmentation to compare against?
        # Should we compare against all segmentations?
        # If there's no exact match, which segmentation should we compare metrics with?
        ground_truth = os.path.join(self.__get_backup_path(index), gt_filename)

        # Prepare the assets for metrics generation
        inputdata_file = os.path.join(self.__get_output_path(index), "inputdata.csv")
        config_file = os.path.join(self.__get_output_path(index), "parameters.yaml")
        data = {"subjectid": id, "prediction": reviewed_file, "target": ground_truth}
        pd.DataFrame(data).to_csv(inputdata_file)
        # TODO: Where do we get this config file?
        # From reading the code, it seems to expect an MLCube parameters.yaml
        # file which was used for training/generating inference
        # That concept breaks here, because we have multiple models running
        # without an accompanying MLCube, and we would need to know which config to use
        # for which model

        # config.yaml can be found inside project/data_prep_models/tumor_segmentation/{model_id}/config.yaml
        config = {"problem_type": "segmentation"}
        with open(config_file, "w") as f:
            yaml.dump(config, f)

        out_file = os.path.join(self.__get_output_path(index), "out.json")

        # Run the metrics generation logic
        generate_metrics.generate_metrics_dict(inputdata_file, config_file, out_file)

        # Open the generated metrics
        with open(out_file, "r") as f:
            metrics = json.load(f)
