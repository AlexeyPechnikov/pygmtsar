#!/usr/bin/env python
# TODO: xarray 2023.3.0 requires pandas<2,>=1.4

from setuptools import setup

# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name='pygmtsar',
    version='2023.4.20',
    description='PyGMTSAR (Python GMTSAR) - Easy and Fast Satellite Interferometry For Everyone',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/mobigroup/gmtsar',
    author='Alexey Pechnikov',
    author_email='pechnikov@mobigroup.ru',
    license='BSD-3-Clause',
    packages=['pygmtsar'],
    install_requires=['xarray>=0.19.0',
                      'importlib-metadata',
                      'numpy>=1.22.4',
                      'pandas>=1.4',
                      'geopandas',
                      'distributed>=2022.11.1',
                      'dask[complete]',
                      'dask_image',
                      'joblib',
                      'tqdm',
                      'sentineleof',
                      'scipy>=1.9.1',
                      'shapely>=2.0.1',
                      'scikit-learn',
                      'elevation',
                      'xmltodict',
                      'rioxarray',
                      'ipywidgets',
                      'h5netcdf==1.0.2',
                      'h5py',
                      'nc-time-axis'
                      ],

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
