# Processing the Data for FeTS

**Note** the `${fets_root_dir}` from [Setup](./setup.md#set-up-the-environment).

For the first application of FeTS in volumetric brain tumor MRI scans, you should follow the pre-processing pipeline defined in the [International Brain Tumor Segmentation (BraTS) Challenge](http://braintumorsegmentation.org/):
- Download the DICOM files of the structural multi-parametric MRI scans:
  - T1-weighted pre-contrast (T1)
  - T1-weighted post-contrast (T1Gd)
  - T2-weighted (T2)
  - T2 Fluid Attenuated Inversion Recovery (T2-FLAIR)
- Convert the DICOM files to NIfTI: feel free to use your preferred tool or via the `Utilities` executable:
```bash
${fets_root_dir}/bin/Utilities -i C:/test/1.dcm -o C:/test.nii.gz -d2n
```
- The complete official BraTS pre-processing pipeline is provided through FeTS in a single executable called `BraTSPipeline`, which performs the following steps:
  
  1. Re-orientation to LPS/RAI
  2. N4 Bias correction (only for facilitating the registration process)
  3. Co-registration to T1Gd
  4. Registration to [SRI-24 anatomical atlas](https://www.nitrc.org/projects/sri24/) (Note that the 2 registration steps are applied as one transformation matrix, thereby avoiding multiple intensity interpolations)
  5. Application of the registration/transformation matrix to the re-oriented image (prior to N4 bias correction) to maximize image fidelity

  ```bash
  ${fets_root_dir}/bin/BraTSPipeline \
    -t1 C:/test/t1.nii.gz \ # you can pass first dicom image series 
    -t1c C:/test/t1gd.nii.gz \ # you can pass first dicom image series 
    -t2 C:/test/t2.nii.gz \ # you can pass first dicom image series 
    -fl C:/test/flair.nii.gz \ # you can pass first dicom image series 
    -o C:/test/outputDir 
  ```
- After finishing the pre-processing, the data needs to be organized in the BraTS format:
```
/data_folder/ 
            /patient_1/
                      /patient_1_t1.nii.gz
                      /patient_1_t2.nii.gz
                      /patient_1_t1ce.nii.gz
                      /patient_1_flair.nii.gz
                      /patient_1_seg.nii.gz
            /patient_2/
                      ...
            ...
            /patient_n/
```
