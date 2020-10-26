# Application Setup

## Requirements

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
  - 11GB dedicated VRAM
  - 40GB RAM
- 80GB RAM for CPU-only tasks
- Read/write access to the data for processing
- Data requirements for Phase-1 Validation: 
  - Glioblastoma (GBM) patients
  - 4 structural modalities - T1 pre and post contrast, T2 and Flair
  - Scanned along axial axis
  - At least 90 cases
  - Tumor regions (as defined in BraTS challenge):
    - Edema - ED
    - Enhancing - ET
    - Necrotic + Non-Enhancing-Tumor core - NET
  - No instrumentation, pre-resection

## Set up the Environment

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

## Set up the Collaborator

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
