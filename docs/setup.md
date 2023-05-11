# Application Setup

## Table of Contents
- [Application Setup](#application-setup)
  - [Table of Contents](#table-of-contents)
  - [Requirements](#requirements)
    - [Windows](#windows)
  - [Set up the Environment](#set-up-the-environment)
    - [Troubleshooting](#troubleshooting)
    - [Note for Ubuntu 20.04 users](#note-for-ubuntu-2004-users)
    - [Optional instructions for Federation backend](#optional-instructions-for-federation-backend)
  - [Set up the Collaborator](#set-up-the-collaborator)

## Requirements

- OS: Linux (Ubuntu 18.04) [note that Ubuntu 20.04 does **not** work, and we will support Windows for Phase-2]
- Python 3.6+: we have tested on [Python distributed by Anaconda](https://www.anaconda.com/products/individual) 3.6, 3.7, 3.8.
Example installation of non-Anaconda Python on Ubuntu 18.04:
```bash
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt update
sudo apt install python3.6 python3.6-venv python3.6-dev python3-setuptools
```
  - When using Anaconda, simply activate an environment with either of the following versions of python **before** starting the installation (ensure that the command `python3 --version` returns the correct information):
    - 3.6.5
    - 3.7
    - 3.8
  - **NOTE**: Python 3.9 does **not** work with underlying dependencies
- GPU: for faster training and inference
  - [CUDA 9.2 - 10.2](https://developer.nvidia.com/cuda-toolkit)
    - **Note**: [NVIDIA Ampere cards](https://www.nvidia.com/en-us/data-center/ampere-architecture/) (including all 30X consumer series cards) require [PyTorch 1.9+](https://pytorch.org/get-started/locally/) installed with CUDA 11.
  - 11GB dedicated VRAM
  - 40GB RAM (**Note**: 120G if you want to run [DeepScan](https://doi.org/10.1007/978-3-030-11726-9_40) inference)
- Read/write access to the data for processing
- Data requirements: 
  - Glioblastoma (GBM) patients
  - 4 structural modalities - T1 pre and post contrast, T2 and T2-Flair
  - No instrumentation, pre-resection
  - Scanned along axial axis
  - At least 90 cases
  - Tumor regions (as defined in BraTS challenge):
    -	Edematous/Invaded tissue + Non-Enhancing-Tumor - ED
    - Enhancing - ET
    - Necrotic tumor core - NET

### Windows 

Since FeTS is Linux-only at the moment, Windows users will need to enable [Windows subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/install-win10). Here are some detailed instructions to get WSL to recognize the NVIDIA GPU:

- Check the build information for your OS by searching for "about PC" in start menu or running `systeminfo | findstr /B /C:"OS Version"` in powershell.
- If your build is `< 20145`:
	- Enroll in [Windows Insider Program](https://insider.windows.com/en-us/) by searching for `Windows Insider Program` in settings.
    - Once enrollment is finished, activate the **Dev** channel.
	- This will take some time, so once you see the updates (in Windows Update), ensure you finish the download and start installation, following which you can happily go for breakfast/lunch/dinner.
	- Verify that build number (search for "about PC" in start menu or run `systeminfo | findstr /B /C:"OS Version"` in powershell) is `>= 20145`.
- Install [Ubuntu WSL 18.04](https://www.microsoft.com/en-us/p/ubuntu-1804-lts/9n9tngvndl3q).
- Ensure WSL can see the GPU:
	- Download appropriate drivers from https://developer.nvidia.com/cuda/wsl/download (you will need a free NVIDIA developer account).
	- Install the driver.
	- Test that WSL version is greater than `4.19.121` using the command `wsl cat /proc/version` in the powershell terminal.
  - Run the following commands:
  ```bash
  sudo apt-key adv --fetch-keys http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/7fa2af80.pub
  sudo sh -c 'echo "deb http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64 /" > /etc/apt/sources.list.d/cuda.list'
  sudo apt-get update
  sudo apt-get install -y cuda-toolkit-11-0
  ```
	- This will also take some time (enough for a coffee/tea)
- Verify that WSL can see the GPU:
  ```bash
	cd /usr/local/cuda/samples/4_Finance/BlackScholes
	sudo make
	./BlackScholes
  ```
	- If you see "test passed", that means all is well with the world and we can now move on to do some real work.
	- **NOTE**: For DL training, ensure that you install PyTorch with CUDA 11


[Back To Top &uarr;](#table-of-contents)

## Set up the Environment

- Download the pre-built binaries from [the latest releases page](https://github.com/FETS-AI/Front-End/releases)
- Run the following commands:
```bash
sudo apt-get install wget zip unzip # required to download and unzip model weights
cd ${download_location}
chmod +x ./FeTS_${version}_Installer.bin # optional addition of execution permission
./FeTS_${version}_Installer.bin --target ${install_path} # change ${install_path} to appropriate location
# accept license
cd ${install_path}/squashfs-root/usr/ # this is the ${fets_root_dir}
cd bin/OpenFederatedLearning
nvidia-smi # note the ${cuda_version}
# install pytorch using appropriate ${cuda_version}; see https://pytorch.org/get-started/locally/
./venv/bin/pip install torch torchvision # installs latest stable pytorch using CUDA 10.2 (we have tested with 1.7.1 with cuda-10.2)
# for cuda 9.2, this would be './venv/bin/pip install torch==1.7.1+cu92 torchvision==0.8.2+cu92 torchaudio==0.7.2 -f https://download.pytorch.org/whl/torch_stable.html'
# for cuda 10.1, this would be './venv/bin/pip install torch==1.7.1+cu101 torchvision==0.8.2+cu101 torchaudio==0.7.2 -f https://download.pytorch.org/whl/torch_stable.html'
```
### Troubleshooting

If you run into the following error (or something similiar, related to the cryptography package):
```bash
Command "python setup.py egg_info" failed with error code 1 in /tmp/pip-build-m0u0eez3/cryptography/
```
Run the following commands for the solution:
```bash
cd ${fets_root_dir}/bin/OpenFederatedLearning
./venv/bin/python3 -m pip install -U pip
make install_openfl
make install_fets
cd ../LabelFusion
./venv/bin/python3 -m pip install -U pip
./venv/bin/pip install .
```
### Note for Ubuntu 20.04 users

We have **not** tested with Ubuntu 20.04 and there might be unforseen stability issues with dependencies. That being said, if there is no other way around it, there are some pointers that can be followed to get FeTS up and running on this platform:

```bash
echo "deb http://archive.ubuntu.com/ubuntu xenial main" | sudo tee /etc/apt/sources.list.d/xenial.list
sudo apt-get update
sudo update-alternatives --remove-all gcc
sudo update-alternatives --remove-all g++
sudo apt-get install gcc-5 g++-5
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 10
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-5 10
```

### Optional instructions for Federation backend

These commands are run along with the installer, but in case you receive an error during the python environment setup, please follow these instructions:
```bash
cd ${fets_root_dir}
cd bin/OpenFederatedLearning
make install_openfl 
./venv/bin/pip install opencv-python==4.2.0.34 # https://stackoverflow.com/a/63669919/1228757
make install_fets
# after this, the federated backend is ready
./venv/bin/pip install -e ./submodules/fets_ai/Algorithms/GANDLF # gandlf
cd ../LabelFusion
python3 -m venv venv
./venv/bin/pip install -e . # label fusion
```

[Back To Top &uarr;](#table-of-contents)

## Set up the Collaborator

**NOTE:** We are currently working on a follow-up study, and we will be releasing new information on how to perform certificate signing soon!

[Back To Top &uarr;](#table-of-contents)

---
[Next: Process Data](./process_data.md)

---
