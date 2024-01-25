#!/usr/bin/env python
# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2021, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------

from setuptools import setup
import urllib.request

# read the contents of your README file
#from pathlib import Path
#this_directory = Path(__file__).parent
#long_description = (this_directory / "README.md").read_text()

upstream_url = 'https://raw.githubusercontent.com/mobigroup/gmtsar/pygmtsar2/README.md'
response = urllib.request.urlopen(upstream_url)
long_description = response.read().decode('utf-8')

setup(
    name='pygmtsar',
    version='2024.1.21.post3',
    description='PyGMTSAR (Python GMTSAR) - Easy and Fast Satellite Interferometry For Everyone',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/mobigroup/gmtsar',
    author='Alexey Pechnikov',
    author_email='pechnikov@mobigroup.ru',
    license='BSD-3-Clause',
    packages=['pygmtsar'],
    install_requires=['xarray>=0.19.0',
                      'numpy>=1.22.4',
                      'numba',
                      'pandas>=1.5',
                      'geopandas',
                      'distributed>=2022.11.1',
                      'dask[complete]',
                      'dask_image',
                      'joblib',
                      'tqdm',
                      'scipy>=1.9.1',
                      'shapely>=2.0.1',
                      'scikit-learn',
                      'xmltodict',
                      'rioxarray',
                      'ipywidgets',
                      'h5netcdf>=1.2.0',
                      'h5py',
                      'nc-time-axis',
                      'statsmodels>=0.13.5',
                      'pygmt',
                      'remotezip',
                      'asf_search',
                      'imageio',
                      'matplotlib',
                      'adjustText',
                      'seaborn'
                      ],
    extras_require={
                      'vtk_support': ['vtk', 'panel']
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: POSIX',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.8'
    ],
    python_requires='>=3.8'
)
