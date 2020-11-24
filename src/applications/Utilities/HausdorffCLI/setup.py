#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup

setup(name='Hausdorff95',
      version='1.0.0.Alpha',
      description='Get Hausdorff 95 from 2 annotations',
      url='https://github.com/sarthakpati/test_hausdorff',
      python_requires='>=3.6',
      author='Sarthak Pati',
      author_email='software@cbica.upenn.edu',
      license='BSD-3-Clause',
      zip_safe=False,
      install_requires=[
      'MedPy==0.4.0',
      'setuptools>=47',
      'PyInstaller==3.6'
      ],
      scripts=['Hausdorff95.py'],
      classifiers=[
          'Intended Audience :: Science/Research',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering',
          'Operating System :: Unix'
      ]
      )