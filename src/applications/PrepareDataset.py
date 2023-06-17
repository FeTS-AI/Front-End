import os, argparse, sys, csv, platform, subprocess, shutil
from pathlib import Path
from datetime import date
from tqdm import tqdm
import pandas as pd


def parse_csv_header(filename):
    """
    Read filename and return a list of dictionaries that have the csv contents
    """
    with open(filename, "r") as csvfile:
        datareader = csv.reader(csvfile)

        headers = {}  # save headers
        for row in datareader:
            for col in row:
                temp = col.lower()  # convert to lower case
                temp = temp.replace(" ", "")  # remove spaces
                temp = temp.replace("_", "")  # remove underscores
                temp = temp.replace("-", "")  # remove dashes
                if (
                    (temp == "patientid")
                    or (temp == "subjectid")
                    or (temp == "subject")
                    or (temp == "subid")
                ):
                    headers["ID"] = col
                elif (
                    (temp == "timepoint")
                    or (temp == "tp")
                    or (temp == "time")
                    or (temp == "series")
                    or (temp == "subseries")
                ):
                    headers["TIMEPOINT"] = col
                elif (temp == "t1gd") or (temp == "t1ce") or (temp == "t1post"):
                    headers["T1GD"] = col
                elif (temp == "t1") or (temp == "t1pre"):
                    headers["T1"] = col
                elif temp == "t2":
                    headers["T2"] = col
                elif (
                    (temp == "t2flair")
                    or (temp == "flair")
                    or (temp == "fl")
                    or ("fl" in temp)
                    or ("t2fl" in temp)
                ):
                    headers["FLAIR"] = col
            break

    if "TIMEPOINT" not in headers:
        headers["TIMEPOINT"] = None
    return headers


def copyFilesToCorrectLocation(interimOutputDir, finalSubjectOutputDir, subjectID):
    """
    This function copies the intermediate files and final outputs to correct location and if these are absent, returns a bool flag stating that brats pipeline needs to run again
    """

    # copy files to correct location for inference and training
    runBratsPipeline = False
    output_dict = {"ID": subjectID}
    output_t1c_brain_file_inter = os.path.join(interimOutputDir, "brain_T1CE.nii.gz")
    output_t1c_brain_file_final = os.path.join(
        finalSubjectOutputDir, subjectID + "_brain_t1ce.nii.gz"
    )
    output_dict["T1GD"] = output_t1c_brain_file_final
    if not os.path.exists(output_t1c_brain_file_final):
        if os.path.exists(output_t1c_brain_file_inter):
            shutil.copyfile(output_t1c_brain_file_inter, output_t1c_brain_file_final)
        else:
            output_t1c_brain_file_inter = os.path.join(
                interimOutputDir, "brain_T1GD.nii.gz"
            )
            if os.path.exists(output_t1c_brain_file_inter):
                shutil.copyfile(
                    output_t1c_brain_file_inter, output_t1c_brain_file_final
                )
            else:
                runBratsPipeline = True

    output_t1_brain_file_inter = os.path.join(interimOutputDir, "brain_T1.nii.gz")
    output_t1_brain_file_final = os.path.join(
        finalSubjectOutputDir, subjectID + "_brain_t1.nii.gz"
    )
    output_dict["T1"] = output_t1_brain_file_final
    if not os.path.exists(output_t1_brain_file_final):
        if os.path.exists(output_t1_brain_file_inter):
            shutil.copyfile(output_t1_brain_file_inter, output_t1_brain_file_final)
        else:
            runBratsPipeline = True

    output_t2_brain_file_inter = os.path.join(interimOutputDir, "brain_T2.nii.gz")
    output_t2_brain_file_final = os.path.join(
        finalSubjectOutputDir, subjectID + "_brain_t2.nii.gz"
    )
    output_dict["T2"] = output_t2_brain_file_final
    if not os.path.exists(output_t2_brain_file_final):
        if os.path.exists(output_t2_brain_file_inter):
            shutil.copyfile(output_t2_brain_file_inter, output_t2_brain_file_final)
        else:
            runBratsPipeline = True

    output_fl_brain_file_inter = os.path.join(interimOutputDir, "brain_FL.nii.gz")
    output_fl_brain_file_final = os.path.join(
        finalSubjectOutputDir, subjectID + "_brain_flair.nii.gz"
    )
    output_dict["FLAIR"] = output_fl_brain_file_final
    if not os.path.exists(output_fl_brain_file_final):
        if os.path.exists(output_fl_brain_file_inter):
            shutil.copyfile(output_fl_brain_file_inter, output_fl_brain_file_final)
        else:
            runBratsPipeline = True

    return runBratsPipeline, output_dict


def main():
    copyrightMessage = (
        "Contact: admin@fets.ai\n\n"
        + "This program is NOT FDA/CE approved and NOT intended for clinical use.\nCopyright (c) "
        + str(date.today().year)
        + " University of Pennsylvania. All rights reserved."
    )
    parser = argparse.ArgumentParser(
        prog="PrepareDataset",
        formatter_class=argparse.RawTextHelpFormatter,
        description="This application calls the BraTSPipeline for all input images and stores the final and intermediate files separately.\n\n"
        + copyrightMessage,
    )
    parser.add_argument(
        "-inputCSV",
        type=str,
        help="The absolute, comma-separated paths of labels that need to be fused",
        required=True,
    )
    parser.add_argument(
        "-outputDir",
        type=str,
        help="The output file to write the results",
        required=True,
    )

    args = parser.parse_args()
    outputDir_qc = os.path.normpath(args.outputDir + "/DataForQC")
    outputDir_final = os.path.normpath(args.outputDir + "/DataForFeTS")

    Path(args.outputDir).mkdir(parents=True, exist_ok=True)
    Path(outputDir_qc).mkdir(parents=True, exist_ok=True)
    Path(outputDir_final).mkdir(parents=True, exist_ok=True)

    bratsPipeline_exe = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "BraTSPipeline"
    )
    if platform.system() == "Windows":
        bratsPipeline_exe += ".exe"

    assert os.path.exists(
        bratsPipeline_exe
    ), "BraTS Pipeline executable not found, please contact admin@fets.ai for help."

    # only parse the headers here
    parsed_headers = parse_csv_header(args.inputCSV)

    # use pandas for this
    subjects_df = pd.read_csv(args.inputCSV)

    output_dict_for_writing_csv = {
        "ID": [],
        "TIMEPOINT": [],
        "T1": [],
        "T1GD": [],
        "T2": [],
        "FLAIR": [],
    }

    for row in tqdm(subjects_df.iterrows(), total=subjects_df.shape[0]):
        subject_id_timepoint = row[parsed_headers["ID"]]
        # joining timepoint to subjectid, but can create a new folder called timepoint if needed
        if parsed_headers["TIMEPOINT"] is not None:
            subject_id_timepoint += "_" + row[parsed_headers["TIMEPOINT"]]
        interimOutputDir_actual = os.path.join(outputDir_qc, subject_id_timepoint)
        finalSubjectOutputDir_actual = os.path.join(
            outputDir_final, subject_id_timepoint
        )
        Path(interimOutputDir_actual).mkdir(parents=True, exist_ok=True)
        Path(finalSubjectOutputDir_actual).mkdir(parents=True, exist_ok=True)
        # check if the files exist already, if so, skip
        runBratsPipeline, outputs = copyFilesToCorrectLocation(
            interimOutputDir_actual, finalSubjectOutputDir_actual, subject_id_timepoint
        )

        output_dict_for_writing_csv["ID"].append(outputs["ID"])
        ## we don't need this if the timepoint is embedded in the subject id
        # if parsed_headers["TIMEPOINT"] is not None:
        #     output_dict_for_writing_csv["TIMEPOINT"].append(
        #         row[parsed_headers["TIMEPOINT"]]
        #     )
        output_dict_for_writing_csv["T1"].append(outputs["T1"])
        output_dict_for_writing_csv["T1GD"].append(outputs["T1GD"])
        output_dict_for_writing_csv["T2"].append(outputs["T2"])
        output_dict_for_writing_csv["FLAIR"].append(outputs["FLAIR"])

        if runBratsPipeline:
            command = (
                bratsPipeline_exe
                + " -t1 "
                + row[parsed_headers["T1"]]
                + " -t1c "
                + row[parsed_headers["T1GD"]]
                + " -t2 "
                + row[parsed_headers["T2"]]
                + " -fl "
                + row[parsed_headers["FLAIR"]]
                + " -o "
                + interimOutputDir_actual
            )
            print("Command: ", command)
            subprocess.Popen(command, shell=True).wait()

        runBratsPipeline, _ = copyFilesToCorrectLocation(
            interimOutputDir_actual, finalSubjectOutputDir_actual, subject_id_timepoint
        )
        if runBratsPipeline:
            print(
                "BraTSPipeline failed for subject ID:",
                subject_id_timepoint,
                file=sys.stderr,
            )

    output_csv_file = os.path.join(outputDir_final, "processed_data.csv")
    output_df = pd.DataFrame.from_dict(output_dict_for_writing_csv)
    output_df.to_csv(output_csv_file, index=False)


if __name__ == "__main__":
    if platform.system() == "Darwin":
        sys.exit("macOS is not supported")
    else:
        main()
