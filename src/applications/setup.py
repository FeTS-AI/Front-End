#!/usr/bin/env python

"""The setup script."""


import sys, re
from setuptools import setup, find_packages

__version__ = "2.0.0"

requirements = [
    "black",
    "numpy==1.24.0",
    "SimpleITK!=2.0.*",
    "SimpleITK!=2.2.1",  # https://github.com/mlcommons/GaNDLF/issues/536
    "scikit-learn>=0.23.2",
    "scikit-image>=0.19.1",
    "tqdm",
    "setuptools",
    "pyyaml",
    "pytest",
    "pytest-cov",
    "pyradiomics==3.1",
]

if __name__ == "__main__":
    setup(
        name="FeTS_Tool_Helper",
        version=__version__,
        author="FeTS.AI",
        author_email="admin@fets.ai",
        python_requires=">=3.8",
        packages=find_packages(),
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
