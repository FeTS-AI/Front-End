# Processing the Data for FeTS

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
│   └───Patient_001_T1
│       │   image_001.dcm
│       │   image_002.dcm
│       │   ...
│   └───Patient_001_T1GD
│       │   image_001.dcm
│       │   image_002.dcm
│       │   ...
│   └───Patient_001_T2
│       │   image_001.dcm
│       │   image_002.dcm
│       │   ...
│   └───Patient_001_T2FLAIR
│       │   image_001.dcm
│       │   image_002.dcm
│       │   ...
│   
└───Patient_002
│   │ ...
```
- Construct a CSV (let's call this **raw_data.csv**) containing the first DICOM images from each modality:
```
PatientID,T1,T1GD,T2,T2FLAIR
Patient_001,/path/to/T1/image_001.dcm,/path/to/T1GD/image_001.dcm,/path/to/T2/image_001.dcm,/path/to/T2FLAIR/image_001.dcm
Patient_002,/path/to/T1/image_001.dcm,/path/to/T1GD/image_001.dcm,/path/to/T2/image_001.dcm,/path/to/T2FLAIR/image_001.dcm
...
Patient_X,/path/to/T1/image_001.dcm,/path/to/T1GD/image_001.dcm,/path/to/T2/image_001.dcm,/path/to/T2FLAIR/image_001.dcm
```
- Pass **raw_data.csv** as an input, along with an output directory, to the `PrepareDataset` executable (which interally calls the `BraTSPipeline` executable):
```bash
${fets_root_dir}/bin/PrepareDataset -i /path/to/raw_data.csv -o /path/to/output
```
- Two output directories will be created under `/path/to/output`:
  - `DataForFeTS`: this is to be passed for inference/training:
  ```
  DataForFeTS
  │
  └───Patient_001
  │   │ brain_t1.nii.gz
  │   │ brain_t1gd.nii.gz
  │   │ brain_t2.nii.gz
  │   │ brain_t2flair.nii.gz
  │   
  └───Patient_002
  │   │ ...
  ```
  - `DataForQC`: this is to be used for quality-control:
  ```
  DataForQC 
  │
  └───Patient_001
  │   │ raw_${modality}.nii.gz
  │   │ raw_rai_${modality}.nii.gz
  │   │ raw_rai_n4_${modality}.nii.gz
  │   │ ${modality}_to_SRI.nii.gz
  │   │ brainMask_SRI.nii.gz # generated using DeepMedic
  │   │ log.txt
  │   
  └───Patient_002
  │   │ ...
  ```
