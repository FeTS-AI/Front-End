import os
import yaml
import argparse
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser("MedPerf Statistics Example")
    parser.add_argument(
        "--data_path",
        dest="data",
        type=str,
        help="directory containing the prepared data",
    )
    parser.add_argument(
        "--labels_path",
        dest="labels",
    )
    parser.add_argument(
        "--out_file", dest="out_file", type=str, help="file to store statistics"
    )
    parser.add_argument(
        "--metadata_path",
        dest="metadata_path",
        type=str,
        help="path to the local metadata folder",
    )

    args = parser.parse_args()

    dicom_info_file = "dicom_tag_information_to_write_anon.yaml"
    dicom_info_filepath = os.path.join(args.metadata_path, dicom_info_file)

    stats = {}
    if os.path.exists(dicom_info_filepath):
        with open(dicom_info_filepath, "r") as f:
            stats = yaml.safe_load(f)

    with open(args.out_file, "w") as f:
        yaml.dump(stats, f)
