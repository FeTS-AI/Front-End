# FeTS: Federated Tumor Segmentation 

The Federated Tumor Segmentation (FeTS) platform, describes an on-going under development open-source toolkit, with a user-friendly graphical user interface (GUI), aiming at:

- bringing pre-trained segmentation models of numerous deep learning algorithms and their fusion, closer to clinical experts and researchers, thereby enabling easy quantification of new radiographic scans and comparative evaluation of new algorithms.
- allowing secure multi-institutional collaborations via federated learning to improve these pre-trained models without sharing patient data, thereby overcoming legal, privacy, and data-ownership challenges.

Successful completion of this project will lead to an easy-to-use potentially-translatable tool enabling easy, fast, objective, repeatable and accurate tumor segmentation, without requiring a computational background by the user, and while facilitating further analysis of tumor radio-phenotypes towards accelerating discovery. 

FeTS is developed and maintained by the [Center for Biomedical Image Computing and Analytics (CBICA)](https://www.cbica.upenn.edu/) at the University of Pennsylvania, in collaboration with [Intel Labs](https://www.intel.com/content/www/us/en/research/overview.html), [Intel AI](https://www.intel.com/ai), and Intel Internet of Things Group.

For more details, please visit us at https://www.fets.ai/

## Status & Timeline

- FeTS is currently undergoing a Phase-1 evaluation with a limited number of international collaborating institutions.
- The Phase-1 evaluation is expected to go on until the end of Q4 2020.
- Phase-2 evaluation, including all committed international collaborators is expected to follow the end of Phase-1. 

## Supporting Grant

This work is in part supported by the National Institutes of Health / National Cancer Institute / Informatics Technology for Cancer Research (NIH/NCI/ITCR), under grant award number U01-CA242871.

## Application Setup

### Requirements

- OS: Linux (Ubuntu 18.04) [note that Ubuntu 20.04 does **not** work, and we will support Windows for Phase-2]
- Python 3.6+: we have tested on [Python distributed by Anaconda](https://www.anaconda.com/products/individual) 3.6, 3.7
Example installation on Ubuntu 18.04:
```bash
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt update
sudo apt install python3.6 python3.6-venv python3.6-dev
```
- (OPTIONAL) GPU: for faster training and inference
  - [CUDA 9.2+](https://developer.nvidia.com/cuda-toolkit)
  - [Compatible cuDNN](https://developer.nvidia.com/cudnn)
- Read/write access to the data for processing
- Data requirements for Phase-1 Validation: 
	- Glioblastoma (GBM) patients
	- 4 structural modalities - T1 pre and post contrast, T2 and Flair
	- At least 90 cases
	- Tumor regions (Edema - ED, Enhancing - ET and Non-Enhancing- NET)
  - No instrumentation, pre-resection

### Set up the binaries and environment

- Download the pre-built binaries from [this link](https://www.nitrc.org/frs/downloadlink.php/11892)
- Run the following commands:
```bash
cd ${download_location}
chmod +x ./FeTS_${version}.bin # optional addition of execution permission
./FeTS_${version}.bin --target ${install_path} # change ${install_path} to appropriate location
# accept license
cd ${install_path}/squashfs-root/usr/ # this is the ${fets_root_dir}
cd bin/OpenFederatedLearning
make install_openfl 
make install_fets
./venv/bin/pip install torch torchvision # installs latest stable pytorch (we have tested with 1.6.0 with cuda-10.2), change to appropriate cuda version; see https://pytorch.org/get-started/locally/
# for cuda 9.2, this would be './venv/bin/pip install torch==1.6.0+cu92 torchvision==0.7.0+cu92 -f https://download.pytorch.org/whl/torch_stable.html'
# for cuda 10.1, this would be './venv/bin/pip install torch==1.6.0+cu101 torchvision==0.7.0+cu101 -f https://download.pytorch.org/whl/torch_stable.html'
# after this, the federated backend is ready
```

### Setup the collaborator

- Run the following commands to generate the [Certificate Signing Request (CSR)](https://en.wikipedia.org/wiki/Certificate_signing_request):
```bash
cd ${fets_root_dir} # see previous step for where this would be
cd bin/OpenFederatedLearning/bin/federations/pki
bash create-collaborator.sh ${collaborator_common_name} # keep a note of the ${collaborator_common_name}, as it will be used for authentication and to send/receive jobs to/from the aggregator at UPenn
# this command will generate the following items, which needs to be saved for collaborator verification:
## a CSR file at `./client/${collaborator_common_name}.csr`
## a string of alpha-numeric numbers for hash verification 
```

- Send the CSR file (`${fets_root_dir}/bin/OpenFederatedLearning/bin/federations/pki/client/${collaborator_common_name}.csr`), and **not the verificiation hash**, to your UPenn point-of-contact, Sarthak Pati
- Set up a call/meeting with Sarthak and provide the verification hash to him.
  - *On Aggregator*: `bash sign-csr.sh ${collaborator_common_name}.csr ${hashVerification}`
  - Sarthak will send the following file back: `${collaborator_common_name}.crt`
  - Copy this file to the following location: `${fets_root_dir}/bin/OpenFederatedLearning/bin/federations/pki/client/`
  - (**Verification**) This location should contain the following files:
    - `${collaborator_common_name}.crt`
    - `${collaborator_common_name}.csr`
    - `${collaborator_common_name}.key`

### Process the data

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
${fets_root_dir}/bin/BraTSPipeline -t1 C:/test/t1.nii.gz -t1c C:/test/t1gd.nii.gz -t2 C:/test/t2.nii.gz -fl C:/test/flair.nii.gz -o C:/test/outputDir 
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

## Running the Application

### Application Path

```bash
cd ${download_location}
./squashfs-root/usr/bin/FeTS # launches application
```
### Inference

<p align="center">
    <img src="https://github.com/FETS-AI/Front-End/blob/master/docs_sources/images/fets_inference.png?raw=true" />
</p>

- Ensure the data has been pre-processed and is organized in the format described above.
- Click to the **Segmentation** tab in the UI
- Set up the inputs:
  - Select the input directory, which would be `/data_folder` in the above illustration
- From the "Algorithms" group, select the appropritate model (for Phase-1, we are providing pre-trained models using the 3DResUNet architecture):
  - Skull Stripping 
  - Brain Tumor Segmentation
- Select if you want to run the inference on GPU or CPU
- Click "Fusion + Save" to start
- Use the drawing tools provided in the "Drawing" tab for corrections and save the final segmentation.
- **Note**: 
  - For CPU, a run-time requirement of at least 48GB of RAM would be needed
  - For GPU, at least 11GB of dedicated VRAM would be needed

### Training

<p align="center">
    <img src="https://github.com/FETS-AI/Front-End/blob/master/docs_sources/images/fets_training.png?raw=true" />
</p>

- Ensure the data has been pre-processed and is organized in the format described above.
  - Note that there should be **only 1** file that ends in `_seg.nii.gz`, as this will be used for training the Brain Tumor Segmentation model
- Click to the **Segmentation** tab in the UI
- Set up the inputs:
  - Select the input directory, which would be `/data_folder` in the above illustration
  - Type the correct common collaboration name, which would be `${collaborator_common_name}`
- For Phase-1, we currently support training on the 3DResUNet architecture for Brain Tumor Segmentation.
- Select if you want to run the inference on GPU or CPU
- Click "Train + Save" to start

## The Federated Tumor Segmentation (FeTS) and Open Federated Learning (OpenFL)

To enable distributed testing and training machine learning models without direct access to collaborator data, FeTS uses a technique called Federated Learning. The [Open Federated Learning (OpenFL) framework](https://github.com/IntelLabs/OpenFederatedLearning) is developed as part of a collaboration between Intel and the University of Pennsylvania (UPenn), as a part of Intel’s commitment in supporting the grant awarded to the [Center for Biomedical Image Computing and Analytics](https://www.cbica.upenn.edu/) at UPenn (PI: S.Bakas) from the [Informatics Technology for Cancer Research (ITCR)](https://itcr.cancer.gov/) program of the National Cancer Institute (NCI) of the National Institutes of Health (NIH), for the development of the [Federated Tumor Segmentation (FeTS, www.fets.ai) platform](https://www.fets.ai/) (grant award number: U01-CA242871).

## Disclaimer

- The software has been designed for research purposes only and has neither been reviewed nor approved for clinical use by the Food and Drug Administration (FDA) or by any other federal/state agency.
- Certain part of the code for this user interface (excluding dependent libraries) is governed by the license provided in https://www.med.upenn.edu/sbia/software-agreement.html unless otherwise specified.


For more details, please visit us at https://www.fets.ai/

For issues, please visit https://github.com/FETS-AI/Front-End/issues 

## Downloads

By downloading FeTS, you agree to our [License](./LICENSE). 

## Contact

For more information, please contact <a href="mailto:software@cbica.upenn.edu">CBICA Software</a>.

## GitHub Distribution

We currently provide only our tagged versions of the code via GitHub. Check the "tags" using your favorite Git client after cloning our repository. The analogous commands are as follows:

```bash
git clone https://github.com/FETS-AI/Front-End.git
latesttag=$(git describe --tags)
echo checking out ${latesttag}
git checkout ${latesttag}
```
