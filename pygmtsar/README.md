[![GMTSAR tests](https://github.com/mobigroup/gmtsar/actions/workflows/gmtsar.yml/badge.svg)](https://github.com/mobigroup/gmtsar/actions/workflows/gmtsar.yml)
[![MacOS tests](https://github.com/mobigroup/gmtsar/actions/workflows/macos.yml/badge.svg)](https://github.com/mobigroup/gmtsar/actions/workflows/macos.yml)
[![Ubuntu tests](https://github.com/mobigroup/gmtsar/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/mobigroup/gmtsar/actions/workflows/ubuntu.yml)
[![PyPI tests](https://github.com/mobigroup/gmtsar/actions/workflows/pypi.yml/badge.svg)](https://github.com/mobigroup/gmtsar/actions/workflows/pypi.yml)
[![Available on pypi](https://img.shields.io/pypi/v/pygmtsar.svg)](https://pypi.python.org/pypi/pygmtsar/)

**PyGMTSAR** (Python GMTSAR) is an open source project and Python package that provides easy and fast Sentinel-1 Satellite Interferometry for everyone! 

The goal of the project is easy and fast satellite interferometry (InSAR) processing for Sentinel-1 radar scenes everywhere as on localhost as on cloud environments like to Google Cloud VM and AI Notebooks and Amazon EC2 and on free of charge cloud environment Google Colab. GMTSAR binary command line tools are used under the hood but all GMTSAR scripts and GMT command replaced by Python code using modern and robust algorithms.

## Why PyGMTSAR?

PyGMTSAR itself combines powerful Python instrumentary for sophisticated multidementional processing (xarray library) and lazy calculations (dask library) plus parallel computing (joblib library) to perform fast and interactive processing on huge datasets. And the best algorithms and numerical computation approaches applied for all the processing steps. There are progressbars and preview plots for the every long operation and that's easy to save intermediate results and continue work later on the same or other host. For an example, using the dump/restore features some work like to initial raw Sentinel-1 scenes downloading and preprocessing can be performed on a cloud instance and continued on a much smaller subset locally. And (thanks to joblib library) that's safe to interrupt the execution at any time without memory leaks (common for dask and dask-based libraries).

PyGMTSAR is really fast and that's possible to complete SBAS analysis for 5 years on 800 interferograms in just one day even on Apple Air or Apple iMac (8 cores and 16 GB RAM) using 2 TB raw Sentinel-1 scenes. And see the live Google Colab notebooks to find how dramatically PyGMTSAR enhaces the results in comparision to GMTSAR.

PyGMTSAR uses modified GMTSAR command line tools and all the required patches and enhancements pulled into the original GMTSAR project. In case when you have the recent GMTSAR installation PyGMTSAR will work with it. Otherwise, install by the same way PyGMTSAR or GMTSAR  binaries from the GitHub repositories. PyGMTSAR live example Google Colab Notebooks install it automatically and you would just copy the instructions for Ubuntu 18.04 LTS. Also, see Debian 10 installation script in the repository. 

## Learn more

- Documentation: https://github.com/mobigroup/gmtsar

- Issue tracker: https://github.com/mobigroup/gmtsar/issues

- Source code: https://github.com/mobigroup/gmtsar

