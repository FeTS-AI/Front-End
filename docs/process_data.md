# Processing the Data for FeTS

**Note** the `${fets_root_dir}` from [Setup](./setup.md#set-up-the-environment).

## Table of Contents
- [Application Path](#application-path)
- [Data Arrangement](#data-arrangement)
- [Running Pre-processing](#pre-processing)
- [Starting from pre-processed data](#starting-from-pre-processed-data)

## Application Path

```bash
cd ${download_location}
${fets_root_dir}/bin/FeTS # launches application
```

Please add the following path to your `LD_LIBRARY_PATH` when using FeTS: `${fets_root_dir}/lib`:
```bash
export LD_LIBRARY_PATH=${fets_root_dir}/lib:$LD_LIBRARY_PATH
```

[Back To Top &uarr;](#table-of-contents)

## Data Arrangement
**Note** the `${fets_root_dir}` from [Setup](./setup.md#set-up-the-environment).

For the first application of FeTS in volumetric brain tumor MRI scans, you should follow the pre-processing pipeline defined in the [International Brain Tumor Segmentation (BraTS) Challenge](http://braintumorsegmentation.org/):
- Download the DICOM files of the structural multi-parametric MRI scans:
  - T1-weighted pre-contrast (T1)
  - T1-weighted post-contrast (T1Gd)
  - T2-weighted (T2)
  - T2 Fluid Attenuated Inversion Recovery (T2-FLAIR)
- Each subject's DICOM (or NIfTI) image is to be put in a separate folder:
```
Input_Data
│
└───Patient_001
│   │
│   └───T1
│   │   │ image_001.dcm
│   │   │ image_002.dcm
│   │   │ ...
│   └───T1GD
│   │   │ image_001.dcm
│   │   │ image_002.dcm
│   │   │ ...
│   └───T2
│   │   │ image_001.dcm
│   │   │ image_002.dcm
│   │   │ ...
│   └───T2FLAIR
│   │   │ image_001.dcm
│   │   │ image_002.dcm
│   │   │ .....
│   
└───Pat_JohnDoe
│   │ ...
│   
│ ...   
│   
└───SmithJoe
│   │ ...
```
- Ensure that unexpected "special" string characters that could hamper the file system's processing are not present (read [this article](https://www.mtu.edu/umc/services/websites/writing/characters-avoid/#:~:text=Illegal%20Filename%20Characters) for more details). Some examples of this include (and are not limited to): `!`, `@`, `#`, `$`, `%`, `^`, `&`, `*`, `'`, `;`, `:`, `,`, `~`
- Construct a CSV (let's call this **raw_data.csv**) containing the first DICOM images from each modality:
```
PatientID,T1,T1GD,T2,T2FLAIR
Patient_001,/path/to/Patient_001/T1/image_001.dcm,/path/to/Patient_001/T1GD/image_001.dcm,/path/to/Patient_001/T2/image_001.dcm,/path/to/Patient_001/T2FLAIR/image_001.dcm
Pat_JohnDoe,/path/to/Pat_JohnDoe/T1/image_001.dcm,/path/to/Pat_JohnDoe/T1GD/image_001.dcm,/path/to/Pat_JohnDoe/T2/image_001.dcm,/path/to/Pat_JohnDoe/T2FLAIR/image_001.dcm
...
SmithJoe,/path/to/SmithJoe/T1/image_001.dcm,/path/to/SmithJoe/T1GD/image_001.dcm,/path/to/SmithJoe/T2/image_001.dcm,/path/to/SmithJoe/T2FLAIR/image_001.dcm
```
  - One way to do this would be via the following command:
  ```bash
  echo "PatientID,T1,T1GD,T2,T2FLAIR" > raw_data.csv
  for d in $Input_data/*; do sub=`basename $d`; t1=`ls -1 $d/T1/* | head -n1`; tlce=`ls -1 $d/T1GD/* | head -n1`; t2=`ls -1 $d/T2/* | head -n1`; flair=`ls -1 $d/T2FLAIR/* | head -n1`; echo $sub,$t1,$tlce,$t2,$flair >> raw_data.csv; done
  ```
  - Please note that subject IDs should be unique strings:
    - **Problematic**: "Patient_1", "Patient_2", ..., "Patient_10", ...,"Patient_20", ...,"Patient_100"
    - **Acceptable**: "Patient_001", "Patient_002", ..., "Patient_010", ..., "Patient_020", ..., "Patient_100", ...

## Pre-processing

- Pass **raw_data.csv** as an input, along with an output directory, to the `PrepareDataset` executable (which internally calls the `BraTSPipeline` executable):
```bash
${fets_root_dir}/bin/PrepareDataset -i /path/to/raw_data.csv -o /path/to/output
```
- Two output directories will be created under the specified output directory with the following structure::
  ```bash
  /path/to/output
  │   │
  │   └───DataForFeTS # this is to be passed for inference/training
  │   │   │
  │   │   └───Patient_001 # this is constructed from the ${PatientID} header of CSV
  │   │   │   │ Patient_001_brain_t1.nii.gz
  │   │   │   │ Patient_001_brain_t1ce.nii.gz
  │   │   │   │ Patient_001_brain_t2.nii.gz
  │   │   │   │ Patient_001_brain_flair.nii.gz
  │   │   │   
  │   │   └───Pat_JohnDoe # this is constructed from the ${PatientID} header of CSV
  │   │   │   │ ...
  │   │
  │   │
  │   └───DataForQC # this is to be used for quality-control
  │   │   │
  │   │   └───Patient_001 # this is constructed from the ${PatientID} header of CSV
  │   │   │   │ raw_${modality}.nii.gz
  │   │   │   │ raw_rai_${modality}.nii.gz
  │   │   │   │ raw_rai_n4_${modality}.nii.gz
  │   │   │   │ ${modality}_to_SRI.nii.gz
  │   │   │   │ brainMask_SRI.nii.gz # generated using BrainMaGe [https://github.com/CBICA/BrainMaGe/] or DeepMedic [https://cbica.github.io/CaPTk/seg_DL.html]
  │   │   │   │ log.txt
  │   │   │   
  │   │   └───Pat_JohnDoe # this is constructed from the ${PatientID} header of CSV
  │   │   │   │ ...
  ```

**NOTE**: For some OS variants, we have seen `PrepareDataset` executable to cause issues, for which we have an alternative with `${fets_root_dir}/bin/PrepareDataset.py`, which has the exact same API and can be invoked in the following way:
```bash
cd ${fets_root_dir}/bin
./OpenFederatedLearning/venv/bin/python \ # virtual environment that was set up in previous section
  ./PrepareDataset.py -i /path/to/raw_data.csv -o /path/to/output
```

[Back To Top &uarr;](#table-of-contents)

## Starting from Pre-processed data

**NOTE**: Skip this step if you have used ```PrepareDataset``` as described in https://fets-ai.github.io/Front-End/process_data#pre-processing

If you have processed data from a prior study that you would like to include in the FeTS federation, please ensure that all the data is co-registered within each patient and the annotations are in the same space. Once that is assured, follow these steps:
1. [Arrange your data](#data-arrangement) by using the processed image files in respective columns
2. Run `PrepareDataset` [as shown above](#pre-processing)
3. In `DataForQC`, under each patient, the transformation matrices will be generated per modality. Use the **T1CE_to_SRI.mat** file (the assumption here is that the data is co-registered within each patient) to transform the annotation (which is in the patient space) in the following manner:
```bash
${fets_root_dir}/bin/Preprocessing \
  -i /path/to/patient_X/annotation.nii.gz \
  -rFI ${fets_root_dir}/data/sri24/atlasImage.nii.gz \
  -o /path/to/output/DataForFeTS/patient_X/annotation_final_seg.nii.gz \
  -rIA /path/to/output/DataForQC/patient_X/T1CE_to_SRI.mat \
  -reg Rigid -rIS 1 -rSg 1
```
4. Load the transformed images and corresponding annotation to ensure they are aligned correctly:
  - `/path/to/output/DataForFeTS/patient_X/brain_*`
  - `/path/to/output/DataForFeTS/patient_X/annotation_final_seg.nii.gz`

[Back To Top &uarr;](#table-of-contents)
---
[Next: Run Application](./runningApplication.md)

---
