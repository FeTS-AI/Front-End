# Hausdorff95

A simple CLI that calculates Hausdorff-95 distance for 2 binary images (to address https://github.com/CBICA/CaPTk/issues/1130)

## Installation

```powershell
conda create -p ./venv python=3.6.5 -y
conda activate ./venv
pip install -e .
```

## Usage

```powershell
python ./Hausdorff95.py -gt ./data/gt.nii.gz -m ./data/mask.nii.gz
```

## Compiling into executable

Follow [instructions](##Installation) to setup the environment and activate it.

```powershell
python runPyInstaller.py
```

The full-contained executable should be generated under `./dist/`.