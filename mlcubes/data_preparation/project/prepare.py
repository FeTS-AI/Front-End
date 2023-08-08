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
from stages.constants import INTERIM_FOLDER, FINAL_FOLDER, TUMOR_MASK_FOLDER


def find_csv_filenames(path_to_dir, suffix=".csv"):
    filenames = os.listdir(path_to_dir)
    return [filename for filename in filenames if filename.endswith(suffix)]


def write_report(report: pd.DataFrame, filepath: str):
    report_dict = report.to_dict()
    with open(filepath, "w") as f:
        yaml.dump(report_dict, f)


def cleanup(path: str):
    walk = list(os.walk(path))
    for path, _, _ in walk[::-1]:
        if len(os.listdir(path)) == 0:
            os.rmdir(path)


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

    return parser.parse_args()


if __name__ == "__main__":
    args = setup_argparser()
    out_data_csv = os.path.join(args.data_out, "data.csv")

    # DETERMINE CURRENT STATUS
    # We will be using paths with keywords to determine the status
    # of each patient
    # Assumption:
    # There's a deterministic way in which I can find the patient
    # for every state, e.g. patient1 will always be named as such
    # inside every step folder
    #
    # step1/
    # ├─ patient3/
    # ├─ patient4/
    # step2/
    # ├─ patient1/
    # step3/
    # ├─ patient2/
    #
    # Given this, we can iterate through the folders to find where the patient is
    # And continue with the procedure
    report = None
    if os.path.exists(args.report):
        with open(args.report, "r") as f:
            report_data = yaml.safe_load(f)
        report = pd.DataFrame(report_data)

    # Generate all paths for all steps

    #########################################################
    # TODO: Turn the enclosed logic into a Stage
    # Check if there's any folder in our data out
    # If not, then that means this is a new execution
    start_new = not next(os.walk(args.data_out))[1]

    # 1. Copy all the data to the out raw path
    out_raw = os.path.join(args.data_out, "raw")
    if start_new:
        print("Copying data to staging")
        shutil.copytree(args.data, out_raw)
    #########################################################

    # If there is a csv file in the input folder
    # always reuse it for the prepared dataset
    csvs = find_csv_filenames(args.data_out)
    if len(csvs) == 1:
        # One csv was found. Assume this is the desired csv
        # move it to the expected location
        # TODO: How to deal with inconsistent paths because of MLCube functionality?
        csv_path = os.path.join(args.data_out, csvs[0])
        os.rename(csv_path, out_data_csv)

    # 2. Generate the report for all identified subjects
    report_gen = GenerateReport(out_raw)
    if report_gen.should_run(report):
        print("Generating new report")
        report = report_gen.execute(report)
        write_report(report, args.report)

    # RUN COLUMN-WISE PROCESSING
    csv_data_out = os.path.join(args.data_out, "validated")
    nifti_data_out = os.path.join(args.data_out, "prepared")
    brain_data_out = os.path.join(args.data_out, "brain_extracted")
    tumor_data_out = os.path.join(args.data_out, "tumor_extracted")
    review_data_out = os.path.join(args.data_out, "reviewed")
    backup_out = os.path.join(args.labels_out, ".tumor_segmentation_backup")
    subjects = list(report.index)
    loop = tqdm(subjects)
    brain_subpaths = [INTERIM_FOLDER, FINAL_FOLDER]
    tumor_subpaths = [TUMOR_MASK_FOLDER]

    csv_proc = AddToCSV(out_raw, out_data_csv, csv_data_out, out_raw)
    nifti_proc = NIfTITransform(out_data_csv, nifti_data_out, csv_data_out, loop)
    brain_extract_proc = Extract(
        out_data_csv,
        brain_data_out,
        brain_subpaths,
        nifti_data_out,
        brain_subpaths,
        loop,
        "extract_brain",
        3,
    )
    tumor_extract_proc = Extract(
        out_data_csv,
        tumor_data_out,
        tumor_subpaths,
        brain_data_out,
        brain_subpaths,
        loop,
        "extract_tumor",
        4,
    )
    manual_proc = ManualStage(out_data_csv, review_data_out, tumor_data_out, backup_out)

    stages = [csv_proc, nifti_proc, brain_extract_proc, tumor_extract_proc, manual_proc]

    for subject in loop:
        for stage in stages:
            if stage.should_run(subject, report):
                loop.set_description(f"{subject} | {stage.get_name()}")
                report = stage.execute(subject, report)
                write_report(report, args.report)

    cleanup(out_raw)
    cleanup(csv_data_out)
    cleanup(nifti_data_out)
    cleanup(brain_data_out)
