import os
import yaml
import argparse
import pandas as pd
from functools import reduce

# Taken from https://code.activestate.com/recipes/577879-create-a-nested-dictionary-from-oswalk/
def get_directory_structure(rootdir):
    """
    Creates a nested dictionary that represents the folder structure of rootdir
    """
    dir = {}
    rootdir = rootdir.rstrip(os.sep)
    start = rootdir.rfind(os.sep) + 1
    for path, dirs, files in os.walk(rootdir):
        folders = path[start:].split(os.sep)
        subdir = dict.fromkeys(files)
        parent = reduce(dict.get, folders[:-1], dir)
        parent[folders[-1]] = subdir
    return dir


def validate_folder_structure(data_path, labels_path):
    data_struct = get_directory_structure(data_path)["mlcube_io0"]
    labels_struct = get_directory_structure(labels_path)["mlcube_io1"]
    
    expected_data_files = ["brain_t1c.nii.gz", "brain_t1n.nii.gz", "brain_t2f.nii.gz", "brain_t2w.nii.gz"]
    expected_labels_files = ["final_seg.nii.gz"]

    assert "splits.csv" in data_struct

    for id in data_struct.keys():
        if data_struct[id] is None:
            # This is a file, ignore
            continue
        for tp in data_struct[id].keys():
            expected_subject_data_files = set(["_".join([id, tp, file]) for file in expected_data_files])
            expected_subject_labels_files = set(["_".join([id, tp, file]) for file in expected_labels_files])

            found_data_files = set(data_struct[id][tp].keys())
            found_labels_files = set(labels_struct[id][tp].keys())

            assert len(expected_subject_data_files - found_data_files) == 0
            assert len(expected_subject_labels_files - found_labels_files) == 0


def sanity_check(data_path: str, labels_path: str, report_df: pd.DataFrame):
    """Runs a few checks to ensure data quality and integrity

    Args:
        data_path (str): Path to data.
        labels_path (str): Path to labels.
        report_df (pd.DataFrame): Report DataFrame, containing information about the preparation
    """
    # Here you must add all the checks you consider important regarding the
    # state of the data
    validate_folder_structure(data_path, labels_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Medperf Model Sanity Check Example")
    parser.add_argument(
        "--data_path",
        dest="data",
        type=str,
        help="directory containing the prepared data",
    )
    parser.add_argument(
        "--labels_path",
        dest="labels",
        type=str,
        help="directory containing the prepared labels",
    )
    parser.add_argument(
        "--report", dest="report", type=str, help="path to the report file"
    )

    args = parser.parse_args()

    with open(args.report, "r") as f:
        report_dict = yaml.safe_load(f)
        report_df = pd.DataFrame(data=report_dict)

    sanity_check(args.data, args.labels, report_df)
