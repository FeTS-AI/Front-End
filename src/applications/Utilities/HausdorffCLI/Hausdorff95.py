#!usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
from medpy.metric.binary import hd95
import SimpleITK as sitk
import pkg_resources

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="Hausdorff95",
        formatter_class=argparse.RawTextHelpFormatter,
        description="\nThis code is used to get the Hausdorff 95th percentile. "
        + "For questions and feedback contact: admin@fets.ai",
    )

    parser.add_argument(
        "-gt",
        dest="groundTruth",
        type=str,
        help="The ground truth image for comparison.\n",
        required=True,
    )

    parser.add_argument(
        "-m",
        dest="maskImage",
        type=str,
        help="The annotated mask to compare against ground truth.\n",
        required=True,
    )

    # parser.add_argument('-v', '--version', action='version',
    #                     version=pkg_resources.require("Hausdorff95")[0].version, help="Show program's version number and exit.") # disabled because pyinstaller doesn't import it properly

    args = parser.parse_args()

    groundTruth = os.path.abspath(args.groundTruth)
    maskImage = os.path.abspath(args.maskImage)

    gt = sitk.GetArrayFromImage(sitk.ReadImage(groundTruth))
    mk = sitk.GetArrayFromImage(sitk.ReadImage(maskImage))

    print(hd95(gt, mk))
