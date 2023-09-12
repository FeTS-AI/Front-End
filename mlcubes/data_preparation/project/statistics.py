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

    args = parser.parse_args()

    # TODO: implement statistics
    stats = {}

    with open(args.out_file, "w") as f:
        yaml.dump(stats, f)
