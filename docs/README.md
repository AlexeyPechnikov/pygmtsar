[![MacOS tests](https://github.com/mobigroup/gmtsar/actions/workflows/macos.yml/badge.svg)](https://github.com/mobigroup/gmtsar/actions/workflows/macos.yml)
[![Ubuntu tests](https://github.com/mobigroup/gmtsar/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/mobigroup/gmtsar/actions/workflows/ubuntu.yml)
[![PyPI tests](https://github.com/mobigroup/gmtsar/actions/workflows/pypi.yml/badge.svg)](https://github.com/mobigroup/gmtsar/actions/workflows/pypi.yml)
[![Available on pypi](https://img.shields.io/pypi/v/pygmtsar.svg)](https://pypi.python.org/pypi/pygmtsar/)
[![Docker](https://badgen.net/badge/icon/docker?icon=docker&label)](https://hub.docker.com/r/mobigroup/pygmtsar)
[![DOI](https://zenodo.org/badge/398018212.svg)](https://zenodo.org/badge/latestdoi/398018212)

## About

PyGMTSAR (Python GMTSAR) is an open-source project and Python package that aims to provide accessible and efficient Sentinel-1 Satellite Interferometry for users. While it is built on a pure Python package, it utilizes the GMTSAR binary tools, which need to be installed.

The main objective of PyGMTSAR is to enable easy and fast satellite interferometry (InSAR) processing using Python scripts and Jupyter Notebooks. It supports Sentinel-1 SLC scenes and can be used on local machines, as well as cloud environments such as Google Cloud VM, AI Notebooks, Amazon EC2, and even the free cloud environment Google Colab. This means that PyGMTSAR-based interferometry processing is readily available in various computing environments, including Google Colab notebooks and Docker images (see below for more information).

Initially, PyGMTSAR was forked from the GMTSAR GitHub repository and underwent significant changes to seamlessly integrate the binary tools within a Python API. Currently, all the modifications developed for the PyGMTSAR project have been merged into GMTSAR. However, PyGMTSAR maintains a feature-rich Python API for interactive and batch computations, while GMTSAR primarily focuses on providing shell scripts for batch processing.

PyGMTSAR is equipped with numerous state-of-the-art features for InSAR data processing, including detrending, flexible weighted and unweighted least-squares processing for SBAS (Small Baseline Subset) time series analysis, Seasonal-Trend decomposition using LOESS (STL), VTK export, and more. These tools are efficiently parallelized and designed to be memory-effective, allowing for effective processing on a wide range of hardware setups, from standard laptops to powerful workstations or servers.

While there are many new features in the roadmap, it's important to note that PyGMTSAR is developed by a solo developer with limited free time. As a result, there is no specific timeline for the development of these features. However, the developer is passionate about exploring new theoretical mathematics and physics approaches and hopes to implement many more functions in the future. 

You can sponsor PyGMTSAR software development and get access to lots of real world use cases on [Patreon](https://www.patreon.com/pechnikov).

## Why PyGMTSAR?

PyGMTSAR offers several compelling reasons for its usage. Firstly, it leverages powerful Python libraries such as xarray for multidimensional processing, dask for lazy calculations and parallel computing, and joblib for efficient parallelization. These libraries enable fast and interactive processing on large datasets, while applying the best algorithms and numerical computation approaches for each step of the processing pipeline.

PyGMTSAR also provides features such as progress bars and preview plots for visualizing the processing steps, allowing users to save intermediate results and resume their work later on the same or different host. Furthermore, the use of joblib ensures that the execution can be safely interrupted at any time without memory leaks, which is a common issue with dask-based solutions.

The combination of powerful Python libraries, optimized algorithms, and user-friendly functions makes PyGMTSAR fast and efficient. With its human-readable and concise code syntax and powerful computing capabilities, PyGMTSAR is a versatile tool that can be used in various domains, from education to research and beyond.

## Documentation

With PyGMTSAR, you can leverage the ready-to-use interactive examples in Live Google Colab notebooks. These examples serve as a starting point for your own analysis and can be easily modified to suit your specific needs. The advantage of using PyGMTSAR in a Live Jupyter notebook environment is that you can immediately access the functionality without the hassle of software installation and configuration.

PyGMTSAR provides self-documented functions using Python docstrings. To access the complete documentation for a specific function, you can use the `help()` function in Jupyter notebook cells or in a Python editor. Simply pass the function name as an argument to `help()` and it will display the docstring, which contains detailed information about the function's usage, parameters, and return values.

By using `help()`, you can easily access the comprehensive documentation for each function in PyGMTSAR and gain a better understanding of its functionality and how to use it effectively in your code.

In addition to the `help()` function, you can also explore the PyGMTSAR sources and track any reported issues on the [PyGMTSAR GitHub](https://github.com/mobigroup/gmtsar) repository. Documentation is available on the [PyGMTSAR GitHub Pages](https://mobigroup.github.io/gmtsar/), providing further guidance on using PyGMTSAR and installing GMTSAR.

### Tutorials: Live Examples in Docker images

You can download the PyGMTSAR Docker image from DockerHub by clicking on the Docker badge [![Docker](https://badgen.net/badge/icon/docker?icon=docker&label)](https://hub.docker.com/r/mobigroup/pygmtsar) Alternatively, you can build the Docker image yourself using the Dockerfiles provided in the [docker](https://github.com/mobigroup/gmtsar/tree/pygmtsar/docker) directory of the PyGMTSAR repository.

For detailed instructions on how to use the Docker images, including pulling the image from DockerHub and running it, you can refer to the DockerHub page linked above. It provides step-by-step guidance on setting up and running the PyGMTSAR Docker image for accessing the live examples.

## Tutorials: Live Examples on Google Colab

By visiting the GitHub page at https://github.com/mobigroup/gmtsar, you'll have access to various notebooks that demonstrate different processing examples. You can click on the notebooks to open them directly in Google Colab and interact with the code. It's a convenient and interactive way to explore the capabilities of PyGMTSAR without the need for local software installation. Enjoy the live examples and delve into the world of satellite interferometry with PyGMTSAR!

## Reference

PyGMTSAR defines a high-level Python `SBAS` class for interferometry processing. The following functions are grouped by common tasks in the order of the full pipeline processing steps. For simple pipelines, some steps can be skipped.

The detailed documentation is available interactively in docstring format. Use `help(function_name)` to read it in your code editor.

### SBAS class functions

Firstly, import the class from the PyGMTSAR Python package as:

```
from pygmtsar import SBAS
```

#### Initialisation (PyGMTSAR original)

- `SBAS`: Initialise `SBAS` object using scenes directory and filters.
- `set_master`: Set master scene for `SBAS` object (the first scene is used by default).

#### Manage orbits (PyGMTSAR original)

- `download_orbits`: Download missed orbits for all the `SBAS` scenes.

#### List scenes and orbits (PyGMTSAR original)

- `to_dataframe`: Return Pandas DataFrame for all `SBAS` scenes.

#### Manage DEM (PyGMTSAR original)

- `set_dem`: Set existing WGS84 DEM file in geographic coordinates for `SBAS` object.
- `download_dem`: Download missed DEM file in geographic coordinates and convert it for the WGS84 ellipsoid.
- `get_dem`: Return WGS84 DEM in geographic coordinates as a single or list of Xarray DataArrays.

#### Manage topography in radar coordinates (combined GMTSAR wrapper and PyGMTSAR original)

- `topo_ra_parallel`: Build WGS84 DEM in radar coordinates for interferogram processing.
- `get_topo_ra`: Return WGS84 DEM in radar coordinates as a single or list of Xarray DataArrays depending on the subswaths count.

#### Framing (combined GMTSAR wrapper and PyGMTSAR original)

- `set_pins`: Set pins for scene framing. Pins are lon/lat points to crop extra bursts.
- `get_pins`: Return pins list. Pins are lon/lat points.
- `reframe_parallel`: Reorder bursts from sequential scenes to cover the full orbit area between pins only.

#### Build Interferogram (GMTSAR wrapper)

- `stack_parallel`: Stack and align scenes.
- `intf_parallel`: Build interferograms for all the subswaths separately.
- `pixel_decimator`: Return function for pixel decimation to the specified output resolution.

#### Merging (GMTSAR wrapper)

- `merge_parallel`: Merge all the separate subswath interferograms into one large interferogram.

#### Manage landmask (PyGMTSAR original)

- `set_landmask`: Set existing land mask file name for `SBAS` object.
- `download_landmask`: Download missed land mask file in geographic coordinates.
- `get_landmask`: Return WGS84 land mask in geographic coordinates as an Xarray DataArray.

#### Unwrapping (SNAPHU wrapper)

- `snaphu_config`: Return SNAPHU text configuration for unwrapping.
- `unwrap_parallel`: Perform parallel phase unwrapping using SNAPHU.

#### Detrending (PyGMTSAR original)

- `detrend_parallel`: Detrend and save to files a set of unwrapped interferograms combining topography and linear components removal and Gaussian filtering.
- `detrend`: Detrend and return output for a single unwrapped interferogram combining topography and linear components removal and Gaussian filtering.

#### SBAS (GMTSAR wrapper)

- `baseline_pairs`: Generate SBAS baseline pairs.
- `sbas`: Perform SBAS unwrapped interferogram timeseries analysis. Refer to GMTSAR.

## Installation

PyGMTSAR project includes GMTSAR binary tools plus Python library and the installation requires as the binaries as the library. Miss the binary installation in case when you already have the recent GMTSAR installed.

### Install PyGMTSAR Python library on MacOS (Apple Silicon and Intel)

[Homebrew](https://brew.sh) package manager is required and it should be installed first, follow the link to do it. As the installation step Apple command line developer tools needt to be installed:

```
xcode-select --install
```

The commands below tested on Intel BigSur and Apple Silicon Monterey systems and perhaps work on older MacOS versions too. Python 3.10+ required for parallel computations on Apple Silicon and it's recommended on Intel Macs too while that's possible to use Python 3.8+ as on Ubuntu 18.04 (bionic). 

```
# Python 3.10 itself
brew install python@3.10
# NetCDF supporting packages
brew install hdf5 netcdf
# install system GIS packages
brew install proj gdal
```

NumPy library can be significantly accelerated (~4x) using [Apple Accelerate Framework](https://developer.apple.com/documentation/accelerate) for using custom compilation as

```
python3 -m pip uninstall -y numpy
python3 -m pip install cython pybind11
# numba requires recent numpy 1.22.4
python3 -m pip install --no-binary :all: --no-use-pep517 numpy==1.22.4
python3 -m pip install numba
# check the output for string '-Wl,Accelerate' to be sure the acceleration works well
python3 -c "import numpy; numpy.show_config()"
```

And install some exact library versions to prevent dependencies conflicts:

```
# this version requires proj (8+) including proj.h
python3 -m pip install cartopy==0.20.2
python3 -m pip install xarray==0.19.0
python3 -m pip install scipy==1.8.1
```

All other packages installation is straightforward:

```
python3 -m pip install \
    h5py netcdf4 h5netcdf \
    rasterio scikit-image \
    distributed zarr nc-time-axis \
    matplotlib seaborn geoviews hvplot datashader bokeh
```

Note: SciPy does not support [Apple Accelerate Framework](https://developer.apple.com/documentation/accelerate) on MacOS and SciPy custom compilation is useless.

Install the library using PIP packages manager as

```
python3 -m pip install pygmtsar
```

Optionally check the installation from you Python shell or Jupyter notebook:

```
from pygmtsar import SBAS
```

### Install PyGMTSAR Python library on Debian 10

Install PyGMTSAR dependencies including system packages and Python libraries:

```
apt install -y python3-pip
python3 -m pip install sentineleof elevation
eio selfcheck

# for python tools and sbas command line arguments calculation
apt install -y gdal-bin python-gdal python3-netcdf4 python3-scipy bc
# for additional Python-coded utilities
python3 -m pip install xarray numpy scipy pytest --upgrade
# for notebooks
apt -y install libgmt5 netcdf-bin python3-pip python3-netcdf4 python3-scipy python3-dev
apt -y install gdal-bin python3-gdal libcharls2 libgdal-dev libproj-dev proj-data proj-bin libgeos-dev
```

Install the library using PIP packages manager as

```
python3 -m pip install pygmtsar
```

Optionally check the installation from you Python shell or Jupyter notebook:

```
from pygmtsar import SBAS
```

### Install PyGMTSAR Python library on Ubuntu 18.04 (bionic)

Install PyGMTSAR dependencies including system packages and Python libraries:

```
# use the old versions vbecause the modern ones does not work properly on the outdated Ubuntu
python3 -m pip install cartopy==0.19.0.post1  1>/dev/null 2>/dev/null
python3 -m pip install install xarray==0.19.0 1>/dev/null 2>/dev/null
python3 -m pip install install scipy==1.7.1   1>/dev/null 2>/dev/null

python3 -m pip install \
    h5py netcdf4 h5netcdf \
    rasterio scikit-image \
    distributed zarr nc-time-axis \
    matplotlib seaborn geoviews hvplot datashader bokeh
```

Install the library using PIP packages manager as

```
pip3 install pygmtsar
```

Optionally check the installation from you Python shell or Jupyter notebook:

```
from pygmtsar import SBAS
```

### Install PyGMTSAR Python library on some other operation system

Install the library using PIP packages manager as

```
python3 -m pip install pygmtsar
```

This command will try to install all the required Python libraries but some system packages might be required to install manually. Probably you know how to do it when you use a rare operation system.

Optionally check the installation from you Python shell or Jupyter notebook:

```
from pygmtsar import SBAS
```

### Install GMTSAR binaries on MacOS (Apple Silicon and Intel)

[Homebrew](https://brew.sh) package manager is required and it should be installed first, follow the link to do it. As the installation step Apple command line developer tools needt to be installed:

```
xcode-select --install
```

And some system dependencies required:

```
brew install wget libtiff hdf5 gmt ghostscript autoconf
```

After that the installation is straightforward, run the same commands for BigSur on Intel chips and Monterey on Apple Silicon chips:

```
# create installation directory
sudo mkdir /usr/local/GMTSAR
sudo chown $(whoami) /usr/local/GMTSAR
# install recent GMTSAR
cd /usr/local
git clone --branch master https://github.com/gmtsar/gmtsar GMTSAR
cd /usr/local/GMTSAR
autoconf
./configure --with-orbits-dir=/tmp
make
make install
```

Note: that's possible to install GMTSAR to /opt directory instead on Apple Silicon Macs.

### Install GMTSAR binaries on Debian 10

```
GMTSAR=/usr/local/GMTSAR
GIT=https://github.com/gmtsar/gmtsar
BRANCH=master

# prepare system
apt update

# https://github.com/gmtsar/gmtsar/wiki/GMTSAR-Wiki-Page
apt install -y csh subversion autoconf libtiff5-dev libhdf5-dev wget
apt install -y liblapack-dev
apt install -y gfortran
apt install -y g++
apt install -y libgmt-dev
apt install -y gmt-dcw gmt-gshhg
# gmt-gshhg-full should be installed automatically (it is required to use GMTSAR landmask)
apt install -y gmt
# fix for missed git and make tools
apt install -y git make

cd $(dirname "$GMTSAR")
git clone --branch "$BRANCH" "$GIT" GMTSAR
cd GMTSAR
autoconf
./configure --with-orbits-dir=/tmp
make
make install
```

### Install GMTSAR binaries on Ubuntu 18.04, 20.04, 22.04

```
apt install -y git csh autoconf make gfortran \
    libtiff5-dev libhdf5-dev liblapack-dev libgmt-dev gmt-dcw gmt-gshhg gmt
cd /usr/local && git clone --branch master https://github.com/gmtsar/gmtsar GMTSAR
cd /usr/local/GMTSAR \
&& autoconf \
&& ./configure --with-orbits-dir=/tmp CFLAGS='-z muldefs' LDFLAGS='-z muldefs' \
&& make \
&& make install
# add the binaries search path
export PATH=/usr/local/GMTSAR/bin:$PATH
# to install PyGMTSAR
apt install -y python3 python3-pip gdal-bin libgdal-dev
```

Note: only Ubuntu 22.04 requires additional flags "-z muldefs" while these are safe to add on any Ubuntu version.

Note: Ubuntu 18.04 (bionic) operation system installed on Google Colab platform and the Live Google Colab notebooks include the installation commands.

### Install GMTSAR binaries on some other operation system

GMTSAR compilation is possible on any POSIX-compatible Unix/Linux/BSD system while some system packages should be installed. See for more recipes [GMTSAR Wiki Page](https://github.com/gmtsar/gmtsar/wiki/GMTSAR-Wiki-Page)



### @ Alexey Pechnikov, May, 2023

<!-- Setup title and logo -->

<script>
document.querySelector("header h1").textContent = 'PyGMTSAR'
document.querySelector("header p").textContent = 'Easy and Fast Satellite Interferometry For Everyone'
//document.querySelector("header p a").textContent = 'View on GitHub mobigroup/gmtsar'
this.img = document.createElement("img");
this.img.style.cssText = 'width: 200px';
this.img.src = "https://user-images.githubusercontent.com/7342379/190416030-5388bf04-e322-4616-9d33-adf84247d976.png";
src = document.querySelector("p.view");
src.appendChild(this.img);
</script>

<style>
div.wrapper {width: 960px;}
section {width: 720px;}
header {width: 230px;}
footer {width: 230px;}
</style>
