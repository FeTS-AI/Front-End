import os
import shutil
from pathlib import Path
import argparse
import pandas as pd
import yaml
from tqdm import tqdm
from stages.generate_report import GenerateReport
from stages.get_csv import AddToCSV
from stages.nifti_transform import NIfTITransform
from stages.extract import Extract
from stages.manual import ManualStage
from stages.comparison import SegmentationComparisonStage
from stages.confirm import ConfirmStage
from stages.split import SplitStage
from stages.pipeline import Pipeline
from stages.constants import INTERIM_FOLDER, FINAL_FOLDER, TUMOR_MASK_FOLDER


def find_csv_filenames(path_to_dir, suffix=".csv"):
    filenames = os.listdir(path_to_dir)
    return [filename for filename in filenames if filename.endswith(suffix)]


def cleanup(path: str):
    walk = list(os.walk(path))
    for path, _, _ in walk[::-1]:
        if len(os.listdir(path)) == 0:
            os.rmdir(path)


def remove_files_in_directory(directory_path):
    files = os.listdir(directory_path)
    for file in files:
        file_path = os.path.join(directory_path, file)
        if os.path.isfile(file_path):  # Check if the path is a file
            os.remove(file_path)


def setup_argparser():
    parser = argparse.ArgumentParser("Medperf Data Preparator Example")
    parser.add_argument(
        "--data_path", dest="data", type=str, help="path containing raw data"
    )
    parser.add_argument(
        "--labels_path", dest="labels", type=str, help="path containing labels"
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

    return parser.parse_args()


def init_pipeline():
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
    split_csv_path = os.path.join(args.data_out, "splits.csv")

    loop = None
    report_gen = GenerateReport(args.data, out_raw)
    csv_proc = AddToCSV(out_raw, out_data_csv, valid_data_out, out_raw)
    nifti_proc = NIfTITransform(out_data_csv, nifti_data_out, valid_data_out, loop)
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
    tumor_extract_proc = Extract(
        out_data_csv,
        tumor_data_out,
        TUMOR_MASK_FOLDER,
        brain_data_out,
        INTERIM_FOLDER,
        # loop,
        "extract_tumor",
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
        args.parameters, args.data_out, args.labels_out, split_csv_path, staging_folders
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
    return Pipeline(report_gen, stages, staging_folders)



if __name__ == "__main__":
    args = setup_argparser()
    out_data_csv = os.path.join(args.data_out, "data.csv")

    report = None
    if os.path.exists(args.report):
        with open(args.report, "r") as f:
            report_data = yaml.safe_load(f)
        report = pd.DataFrame(report_data)

    # 1. If there is a csv file in the input folder
    # always reuse it for the prepared dataset
    csvs = find_csv_filenames(args.data_out)
    if len(csvs) == 1:
        # One csv was found. Assume this is the desired csv
        # move it to the expected location
        # TODO: How to deal with inconsistent paths because of MLCube functionality?
        csv_path = os.path.join(args.data_out, csvs[0])
        os.rename(csv_path, out_data_csv)
        # can we assume the paths inside data.csv to be relative to the csv?
        # TODO: Create some logic to turn the csv paths into the expected paths for the MLCube
        # update_csv_paths(out_data_csv)

    # Generate all paths for all steps

    # TODO: If we identify a folder structure similar to the data-preparation one,
    # copy everything into the folders and start from there
    pipeline = init_pipeline()
    pipeline.run(report, args.report)
