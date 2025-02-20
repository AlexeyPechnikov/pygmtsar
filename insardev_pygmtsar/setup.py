#!/usr/bin/env python
# ----------------------------------------------------------------------------
# InSAR.dev
# 
# This file is part of the InSAR.dev project: https://InSAR.dev
# 
# Copyright (c) 2025, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from setuptools import setup
import urllib.request

def get_version():
    with open("insardev_pygmtsar/__init__.py", "r") as f:
        for line in f:
            if line.startswith("__version__"):
                version = line.split('=')[1]
                version = version.replace("'", "").replace('"', "").strip()
                return version

# read the contents of local README file
#from pathlib import Path
#this_directory = Path(__file__).parent
#long_description = (this_directory / "README.md").read_text()

upstream_url = 'https://raw.githubusercontent.com/AlexeyPechnikov/pygmtsar/pygmtsar2/README.md'
response = urllib.request.urlopen(upstream_url)
long_description = response.read().decode('utf-8')

setup(
    name='insardev_pygmtsar',
    version=get_version(),
    description='InSAR.dev (Python InSAR): PyGMTSAR backend',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/AlexeyPechnikov/pygmtsar',
    author='Alexey Pechnikov',
    author_email='alexey@pechnikov.dev',
    license='BSD-3-Clause',
    packages=['insardev_pygmtsar'],
    include_package_data=True,
    package_data={
        'insardev_pygmtsar': ['data/geoid_egm96_icgem.grd','data/google_colab.sh'],
    },
    install_requires=['insardev_toolkit',
                      'xarray',
                      'numpy',
                      'pandas>=2.2',
                      'geopandas',
                      'distributed>=2024.1.0',
                      'dask[complete]>=2024.4.1',
                      'opencv-python',
                      'scipy',
                      'shapely>=2.0.2',
                      'xmltodict',
                      'rioxarray',
                      'statsmodels>=0.14.0',
                      ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'License :: OSI Approved :: BSD License',
        'Operating System :: POSIX',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11'
    ],
    python_requires='>=3.10',
    keywords='satellite interferometry, InSAR, remote sensing, geospatial analysis, Sentinel-1, SBAS, PSI'
)
