from typing import Union, List, Tuple
from tqdm import tqdm
import pandas as pd
import os
from os.path import realpath, dirname, join
import shutil
import time
import SimpleITK as sitk
import subprocess
import traceback
from LabelFusion.wrapper import fuse_images

from .row_stage import RowStage
from .PrepareDataset import Preparator, FINAL_FOLDER, generate_tumor_segmentation_fused_images, save_screenshot
from .utils import update_row_with_dict, get_id_tp, MockTqdm

MODALITY_MAPPING = {
    "t1c": "t1c",
    "t1ce": "t1c",
    "t1": "t1n",
    "t1n": "t1n",
    "t2": "t2w",
    "t2w": "t2w",
    "t2f": "t2f",
    "flair": "t2f"
}

MODALITY_VARIANTS = {
    "t1c": "T1GD",
    "t1ce": "T1GD",
    "t1": "T1",
    "t1n": "T1",
    "t2": "T2",
    "t2w": "T2",
    "t2f": "FLAIR",
    "flair": "FLAIR"
}

class ExtractNnUNet(RowStage):
    def __init__(
        self,
        data_csv: str,
        out_path: str,
        subpath: str,
        prev_stage_path: str,
        prev_subpath: str,
        status_code: int,
    ):
        self.data_csv = data_csv
        self.out_path = out_path
        self.subpath = subpath
        self.data_subpath = FINAL_FOLDER
        self.prev_path = prev_stage_path
        self.prev_subpath = prev_subpath
        os.makedirs(self.out_path, exist_ok=True)
        self.prep = Preparator(data_csv, out_path, "BraTSPipeline")
        self.pbar = tqdm()
        self.failed = False
        self.exception = None
        self.status_code = status_code

    def get_name(self) -> str:
        return "nnUNet Tumor Extraction"

    def could_run(self, index: Union[str, int], report: pd.DataFrame) -> bool:
        """Determine if case at given index needs to be converted to NIfTI

        Args:
            index (Union[str, int]): Case index, as used by the report dataframe
            report (pd.DataFrame): Report Dataframe for providing additional context

        Returns:
            bool: Wether this stage could be executed for the given case
        """
        prev_paths = self.__get_paths(index, self.prev_path, self.prev_subpath)
        return all([os.path.exists(path) for path in prev_paths])

    def execute(
        self, index: Union[str, int], report: pd.DataFrame
    ) -> Tuple[pd.DataFrame, bool]:
        """Runs the pretrained nnUNet models for tumor segmentation

        Args:
            index (Union[str, int]): case index, as used by the report
            report (pd.DataFrame): DataFrame containing the current state of the preparation flow

        Returns:
            pd.DataFrame: Updated report dataframe
        """
        self.__prepare_exec()
        self.__copy_case(index)
        self.__process_case(index)
        report, success = self.__update_state(index, report)
        self.prep.write()

        return report, success

    def __prepare_exec(self):

        # Reset the file contents for errors
        open(self.prep.stderr_log, "w").close()

        # Update the out dataframes to current state
        self.prep.read()

    def __get_paths(self, index: Union[str, int], path: str, subpath: str):
        id, tp = get_id_tp(index)
        data_path = os.path.join(path, self.data_subpath, id, tp)
        out_path = os.path.join(path, subpath, id, tp)
        return data_path, out_path

    def __copy_case(self, index: Union[str, int]):
        prev_paths = self.__get_paths(index, self.prev_path, self.prev_subpath)
        copy_paths = self.__get_paths(index, self.out_path, self.prev_subpath)
        for prev, copy in zip(prev_paths, copy_paths):
            shutil.copytree(prev, copy, dirs_exist_ok=True)

    def __get_models(self):
        rel_models_path = "../models/nnUNet_trained_models/nnUNet/3d_fullres"
        models_path = realpath(join(dirname(__file__), rel_models_path))
        return os.listdir(models_path)

    def __get_mod_order(self, model):
        rel_orders_path = "../models/nnUNet_modality_order"
        order_path = realpath(join(dirname(__file__), rel_orders_path, model, "order"))
        with open(order_path, "r") as f:
            order_str = f.readline()
        # remove 'order = ' from the splitted list
        modalities = order_str.split()[2:]
        modalities = [MODALITY_MAPPING[mod] for mod in modalities]
        return modalities

    def __prepare_case(self, path, id, tp, order):
        tmp_subject = f"{id}-{tp}"
        tmp_path = os.path.join(path, "tmp-data")
        tmp_subject_path = os.path.join(tmp_path, tmp_subject)
        tmp_out_path = os.path.join(path, "tmp-out")
        shutil.rmtree(tmp_path, ignore_errors=True)
        shutil.rmtree(tmp_out_path, ignore_errors=True)
        os.makedirs(tmp_subject_path)
        os.makedirs(tmp_out_path)
        in_modalities_path = os.path.join(path, "DataForFeTS", id, tp)
        input_modalities = {}
        for modality_file in os.listdir(in_modalities_path):
            if not modality_file.endswith(".nii.gz"):
                continue
            modality = modality_file[:-7].split("_")[-1]
            norm_mod = MODALITY_MAPPING[modality]
            mod_idx = order.index(norm_mod)
            mod_idx = str(mod_idx).zfill(4)

            out_modality_file = f"{tmp_subject}_{mod_idx}.nii.gz"
            in_file = os.path.join(in_modalities_path, modality_file)
            out_file = os.path.join(tmp_subject_path, out_modality_file)
            input_modalities[MODALITY_VARIANTS[modality]] = in_file
            shutil.copyfile(in_file, out_file)

        return tmp_subject_path, tmp_out_path, input_modalities
    
    def __run_model(self, model, data_path, out_path):
        # models are named Task<ID>_..., where <ID> is always 3 numbers
        task_id = model[4:7]
        cmd = f"nnUNet_predict -i {data_path} -o {out_path} -t {task_id} -f all"
        print(cmd)
        print(os.listdir(data_path))
        start = time.time()
        subprocess.call(cmd, shell=True)
        end = time.time()
        total_time = (end - start)
        print(f"Total time elapsed is {total_time} seconds")

    def __finalize_pred(self, tmp_out_path, out_pred_filepath):
        # We assume there's only one file in out_path
        pred = None
        for file in os.listdir(tmp_out_path):
            if file.endswith(".nii.gz"):
                pred = file

        if pred is None:
            raise RuntimeError("No tumor segmentation was found")

        pred_filepath = os.path.join(tmp_out_path, pred)
        shutil.move(pred_filepath, out_pred_filepath)
        return out_pred_filepath

    def __process_case(self, index: Union[str, int]):
        id, tp = get_id_tp(index)
        subject_id = f"{id}_{tp}"
        models = self.__get_models()
        outputs = []
        images_for_fusion = []
        out_path = os.path.join(self.out_path, "DataForQC", id, tp) 
        out_pred_path = os.path.join(out_path, "TumorMasksForQC")
        os.makedirs(out_pred_path, exist_ok=True)
        for i, model in enumerate(models):
            order = self.__get_mod_order(model)
            tmp_data_path, tmp_out_path, input_modalities = self.__prepare_case(self.out_path, id, tp, order)
            out_pred_filepath = os.path.join(out_pred_path, f"{id}_{tp}_tumorMask_model_{i}.nii.gz")
            try:
                self.__run_model(model, tmp_data_path, tmp_out_path)
                output = self.__finalize_pred(tmp_out_path, out_pred_filepath)
                outputs.append(output)
                images_for_fusion.append(sitk.ReadImage(output, sitk.sitkUInt8))
            except Exception as e:
                self.exception = e
                self.failed = True
                self.traceback = traceback.format_exc()
                return

            #cleanup
            shutil.rmtree(tmp_data_path, ignore_errors=True)
            shutil.rmtree(tmp_out_path, ignore_errors=True)


        fused_outputs = generate_tumor_segmentation_fused_images(images_for_fusion, out_pred_path, subject_id)
        outputs += fused_outputs

        for output in outputs:
            # save the screenshot
            tumor_mask_id = os.path.basename(output).replace(".nii.gz", "")
            save_screenshot(
                input_modalities,
                os.path.join(
                    out_path,
                    f"{tumor_mask_id}_summary.png",
                ),
                output,
            )



    def __update_state(
        self, index: Union[str, int], report: pd.DataFrame
    ) -> Tuple[pd.DataFrame, bool]:
        if self.failed:
            del_paths = self.__get_paths(index, self.out_path, self.subpath)
            report, success = self.__report_failure(index, report)
        else:
            del_paths = self.__get_paths(index, self.prev_path, self.prev_subpath)
            report, success = self.__report_success(index, report)

        for del_path in del_paths:
            shutil.rmtree(del_path, ignore_errors=True)

        return report, success

    def __report_success(
        self, index: Union[str, int], report: pd.DataFrame
    ) -> Tuple[pd.DataFrame, bool]:
        data_path, labels_path = self.__get_paths(index, self.out_path, self.subpath)
        report_data = {
            "status": self.status_code,
            "status_name": "TUMOR_EXTRACT_FINISHED",
            "comment": "",
            "data_path": data_path,
            "labels_path": labels_path,
        }
        update_row_with_dict(report, report_data, index)
        return report, True

    def __report_failure(
        self, index: Union[str, int], report: pd.DataFrame
    ) -> Tuple[pd.DataFrame, bool]:
        prev_data_path, prev_labels_path = self.__get_paths(index, self.prev_path, self.prev_subpath)
        msg = f"{str(self.exception)}: {self.traceback}"

        report_data = {
            "status": -self.status_code,
            "status_name": "TUMOR_EXTRACT_FAILED",
            "comment": msg,
            "data_path": prev_data_path,
            "labels_path": prev_labels_path,
        }
        update_row_with_dict(report, report_data, index)
        return report, False