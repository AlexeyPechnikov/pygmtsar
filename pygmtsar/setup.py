#!/usr/bin/env python

from setuptools import setup

setup(
    name='pygmtsar',
    version='2022.09.14',
    description='PyGMTSAR (Python GMTSAR) - Easy and Fast Satellite Interferometry For Everyone',
    url='https://github.com/mobigroup/gmtsar',
    author='Alexey Pechnikov',
    author_email='pechnikov@mobigroup.ru',
    license='GPL-3.0',
    packages=['pygmtsar'],
    install_requires=['xarray>=0.19.0',
                      'numpy',
                      'pandas',
                      'geopandas',
                      'joblib',
                      'tqdm',
                      'sentineleof',
                      'scipy',
                      'shapely',
                      'sklearn',
                      'elevation',
                      'xmltodict'
                      ],

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)'
        'Operating System :: POSIX',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.7'
    ],
    python_requires='>=3.7'
)
