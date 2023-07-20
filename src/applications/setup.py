#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

__version__ = "2.0.0"

requirements = [
    "black",
    "FigureGenerator==0.0.4",
    "gandlf==0.0.16",
    "labelfusion==1.0.14",
    "numpy==1.22.0",
    "SimpleITK!=2.0.*",
    "SimpleITK!=2.2.1",  # https://github.com/mlcommons/GaNDLF/issues/536
    "scikit-learn>=0.23.2",
    "scikit-image>=0.19.1",
    "tqdm",
    "setuptools",
    "pandas<2.0.0",
    "pyyaml",
    "pytest",
    "pytest-cov",
]

if __name__ == "__main__":
    setup(
        name="FeTS_Tool_Helper",
        version=__version__,
        author="FeTS-AI",
        author_email="admin@fets.ai",
        python_requires=">=3.8",
        packages=find_packages(),
        py_modules=[],
        entry_points={
            "console_scripts": [
                "sanitycheck=SanityCheck:main",
                "preparedataset=PrepareDataset:main",
                "createcsvfordicoms=CreateCSVForDICOMs:main",
            ],
        },
        classifiers=[
            "Development Status :: 5 - Production/Stable",
            "Intended Audience :: Science/Research",
            "Natural Language :: English",
            "Operating System :: OS Independent",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
            "Topic :: Scientific/Engineering :: Medical Science Apps.",
        ],
        description=("Helper scripts for the FeTS Tool."),
        install_requires=requirements,
        license="https://raw.githubusercontent.com/FeTS-AI/Front-End/master/LICENSE",
        long_description="Helper scripts for the FeTS Tool",
        long_description_content_type="text",
        include_package_data=True,
        keywords="brain, mri, neuroimaging, machine learning, federated learning",
        zip_safe=False,
    )


## download various models
import shutil, posixpath
from pathlib import Path
from zipfile import ZipFile
from urllib.request import urlretrieve

models_dir = posixpath.join(Path(__file__).parent.resolve(), "data_prep_models")
Path(models_dir).mkdir(parents=True, exist_ok=True)

brain_extraction_models_zip = posixpath.join(models_dir, "brain_extraction_models.zip")
brain_extraction_models_dir = posixpath.join(models_dir, "brain_extraction_models")

if not Path(brain_extraction_models_zip).exists():
    print("Downloading brain extraction models")
    url = "https://upenn.box.com/shared/static/cp5xz726mtb6gwwym8ydcxmw52zfngun"
    urlretrieve(url, brain_extraction_models_zip)
    z = ZipFile(brain_extraction_models_zip)
    z.extractall(models_dir)
    z.close()
    shutil.move(
        posixpath.join(models_dir, "brain_extraction_models", "brain_extraction_models"),
        brain_extraction_models_dir,
    )
