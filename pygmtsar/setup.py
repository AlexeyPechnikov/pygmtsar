#!/usr/bin/env python

from setuptools import setup

# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name='pygmtsar',
    version='2022.10.9',
    description='PyGMTSAR (Python GMTSAR) - Easy and Fast Satellite Interferometry For Everyone',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/mobigroup/gmtsar',
    author='Alexey Pechnikov',
    author_email='pechnikov@mobigroup.ru',
    license='GPL-3.0',
    packages=['pygmtsar'],
    install_requires=['xarray>=0.19.0',
                      'importlib-metadata<=4.12.0',
                      'numpy',
                      'pandas',
                      'geopandas',
                      'dask',
                      'dask_image',
                      'joblib',
                      'tqdm',
                      'sentineleof',
                      'scipy',
                      'shapely',
                      'sklearn',
                      'elevation',
                      'xmltodict',
                      'rioxarray',
                      'ipywidgets',
                      'h5netcdf',
                      'nc-time-axis'
                      ],

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: POSIX',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.7'
    ],
    python_requires='>=3.7'
)
