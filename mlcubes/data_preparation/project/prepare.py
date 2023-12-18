import os
import argparse
import pandas as pd
import yaml
import shutil
from stages.generate_report import GenerateReport
from stages.get_csv import AddToCSV
from stages.nifti_transform import NIfTITransform
from stages.extract import Extract
from stages.extract_nnunet import ExtractNnUNet
from stages.manual import ManualStage
from stages.comparison import SegmentationComparisonStage
from stages.confirm import ConfirmStage
from stages.split import SplitStage
from stages.pipeline import Pipeline
from stages.constants import INTERIM_FOLDER, FINAL_FOLDER, TUMOR_MASK_FOLDER

MODELS_PATH = "/project/models"


def find_csv_filenames(path_to_dir, suffix=".csv"):
    filenames = os.listdir(path_to_dir)
    return [filename for filename in filenames if filename.endswith(suffix)]


def setup_argparser():
    parser = argparse.ArgumentParser("Medperf Data Preparator Example")
    parser.add_argument(
        "--data_path", dest="data", type=str, help="path containing raw data"
    )
    parser.add_argument(
        "--labels_path", dest="labels", type=str, help="path containing labels"
    )
    parser.add_argument(
        "--models_path", dest="models", type=str, help="path to the nnunet models"
    )
    parser.add_argument(
        "--data_out", dest="data_out", type=str, help="path to store prepared data"
    )
    parser.add_argument(
        "--labels_out",
        dest="labels_out",
        type=str,
        help="path to store prepared labels",
    )
    parser.add_argument(
        "--report", dest="report", type=str, help="path to the report csv file to store"
    )
    parser.add_argument(
        "--parameters",
        dest="parameters",
        type=str,
        help="path to the parameters yaml file",
    )
    parser.add_argument(
        "--metadata_path",
        dest="metadata_path",
        type=str,
        help="path to the local metadata folder"
    )

    return parser.parse_args()


def init_pipeline(args):
    # RUN COLUMN-WISE PROCESSING
    out_raw = os.path.join(args.data_out, "raw")
    valid_data_out = os.path.join(args.data_out, "validated")
    nifti_data_out = os.path.join(args.data_out, "prepared")
    brain_data_out = os.path.join(args.data_out, "brain_extracted")
    tumor_data_out = os.path.join(args.data_out, "tumor_extracted")
    match_data_out = args.labels_out
    backup_out = os.path.join(args.labels_out, ".tumor_segmentation_backup")
    staging_folders = [
        out_raw,
        valid_data_out,
        nifti_data_out,
        brain_data_out,
        tumor_data_out,
        backup_out,
    ]
    out_data_csv = os.path.join(args.data_out, "data.csv")
    trash_folder = os.path.join(args.data_out, ".trash")
    invalid_subjects_file = os.path.join(args.metadata_path, ".invalid.txt")

    loop = None
    report_gen = GenerateReport(out_data_csv, args.data, out_raw, args.labels, args.labels_out, args.data_out, 8, brain_data_out, 3, tumor_data_out, 5)
    csv_proc = AddToCSV(out_raw, out_data_csv, valid_data_out, out_raw)
    nifti_proc = NIfTITransform(out_data_csv, nifti_data_out, valid_data_out, args.metadata_path, args.data_out)
    brain_extract_proc = Extract(
        out_data_csv,
        brain_data_out,
        INTERIM_FOLDER,
        nifti_data_out,
        INTERIM_FOLDER,
        # loop,
        "extract_brain",
        3,
    )
    tumor_extract_proc = ExtractNnUNet(
        out_data_csv,
        tumor_data_out,
        INTERIM_FOLDER,
        brain_data_out,
        INTERIM_FOLDER,
        4,
    )
    manual_proc = ManualStage(out_data_csv, tumor_data_out, tumor_data_out, backup_out)
    match_proc = SegmentationComparisonStage(
        out_data_csv,
        match_data_out,
        tumor_data_out,
        backup_out,
    )
    confirm_proc = ConfirmStage(
        out_data_csv,
        args.data_out,
        args.labels_out,
        tumor_data_out,
        backup_out,
        staging_folders,
    )
    split_proc = SplitStage(
        args.parameters, args.data_out, args.labels_out, staging_folders
    )
    stages = [
        csv_proc,
        nifti_proc,
        brain_extract_proc,
        tumor_extract_proc,
        manual_proc,
        match_proc,
        confirm_proc,
        split_proc
    ]
    return Pipeline(report_gen, stages, staging_folders, [trash_folder], invalid_subjects_file)

def init_report(args) -> pd.DataFrame:
    report = None
    if os.path.exists(args.report):
        with open(args.report, "r") as f:
            report_data = yaml.safe_load(f)
        report = pd.DataFrame(report_data)

    return report


def main():
    args = setup_argparser()

    # Move models to the expected location
    if not os.path.exists(MODELS_PATH):
        shutil.copytree(args.models, MODELS_PATH)

    report = init_report(args)
    pipeline = init_pipeline(args)
    pipeline.run(report, args.report)

if __name__ == "__main__":
    main()