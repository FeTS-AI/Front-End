import yaml
import argparse
import pandas as pd

from stages.utils import has_prepared_folder_structure

def sanity_check(data_path: str, labels_path: str, report_df: pd.DataFrame):
    """Runs a few checks to ensure data quality and integrity

    Args:
        data_path (str): Path to data.
        labels_path (str): Path to labels.
        report_df (pd.DataFrame): Report DataFrame, containing information about the preparation
    """
    # Here you must add all the checks you consider important regarding the
    # state of the data
    assert has_prepared_folder_structure(data_path, labels_path), "The contents of the labels and data don't ressemble a prepared dataset"


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
