[![MacOS tests](https://github.com/mobigroup/gmtsar/actions/workflows/macos.yml/badge.svg)](https://github.com/mobigroup/gmtsar/actions/workflows/macos.yml)
[![Ubuntu tests](https://github.com/mobigroup/gmtsar/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/mobigroup/gmtsar/actions/workflows/ubuntu.yml)
[![PyPI tests](https://github.com/mobigroup/gmtsar/actions/workflows/pypi.yml/badge.svg)](https://github.com/mobigroup/gmtsar/actions/workflows/pypi.yml)
[![Available on pypi](https://img.shields.io/pypi/v/pygmtsar.svg)](https://pypi.python.org/pypi/pygmtsar/)
[![Docker](https://badgen.net/badge/icon/docker?icon=docker&label)](https://hub.docker.com/r/mobigroup/pygmtsar)
[![DOI](https://zenodo.org/badge/398018212.svg)](https://zenodo.org/badge/latestdoi/398018212)

## About

PyGMTSAR (Python GMTSAR) is an open source project and Python package that provides Sentinel-1 Satellite Interferometry for everyone! While it's pure Python package under the hood it uses GMTSAR binary tools which should be installed.

Initially, PyGMTSAR is forked from GMTSAR GitHub repository and lot's of changes are made to seamlessly call all the binary tools from Python API. For now, all the changes developed for PyGMTSAR project merged into GMTSAR and the both projects use the same binary tools  although PyGMTSAR maintains rich Python API for interactive and batch computations and GMTSAR provides a set of shell scripts for the batch processing only. To prevent confusing, below PyGMTSAR means the Python package only and GMTSAR means the binary core tools (while these are available in PyGMTSAR GitHub repository too). 

The goal of the project is easy and fast satellite interferometry (InSAR) interactive and batch processing in Python scripts and Jupyter Notebooks for Sentinel-1 SLC scenes everywhere as on localhost as on cloud environments like to Google Cloud VM and AI Notebooks and Amazon EC2 and even on free of charge cloud environment Google Colab. By this way, PyGMTSAR-based interferometry processing is available even in Google Colab notebooks and in Docker image (see below).

Hint: You can sponsor PyGMTSAR software development on [Patreon](https://www.patreon.com/pechnikov).

## Why PyGMTSAR?

PyGMTSAR itself combines powerful Python instrumentary for sophisticated multidementional processing (xarray library) and lazy calculations (dask library) plus parallel computing (dask and joblib libraries) to perform fast and interactive processing on huge datasets. And the best algorithms and numerical computation approaches applied for all the processing steps. There are progressbars and preview plots for the every step and that's easy to save intermediate results and continue work later on the same or other host. And (thanks to joblib library) that's safe to interrupt the execution at any time without memory leaks (common for dask-based solutions).

Thanks to all the powerful Python libraries and the best used algorithms PyGMTSAR is fast and its possible to complete SBAS analysis for 5 years on 800 interferograms in just one day on Apple Air or Apple iMac (8 cores and 16 GB RAM) using 2 TB Sentinel-1 SLC scenes. And PyGMTSAR is user-friendly providing functions to download the required satellite orbit files and DEM and so on. This combination of the human-readable and short code and powerful computing is the key to use PyGMTSAR everywhere from education and to research and more. 

## Documentation

Satellite interferometry is a complex field of knowledge and usually that's requires a lot of time to understood it enough to build the meaningful results. Really, there is a big gap between GMTSAR github cloning and, for an example, valid SBAS + PSI results computation. PyGMTSAR resolves the software-related complexity providing high-level functions running in Live Jupyter notebooks on Google Colab so you have immediate access to the ready to use and working interactive examples which can be easily modified to produce your own results without hassle with software installation and configuration. You'd start from the self-documented Live Google Colab examples and view PyGMTSAR functions documentation and GMTSAR installation instructions  later.

See the project sources and bug tracker on [PyGMTSAR GitHub](https://github.com/mobigroup/gmtsar) and documentation on [PyGMTSAR GitHub Pages](https://mobigroup.github.io/gmtsar/)

### Tutorials: Live Examples in Docker image

Configure your Docker runtime (Preferences -> Resources tab for Docker Desktop) to use 2 CPU cores and 8 GB RAM or 4 CPU cores and 16 GB RAM and so on. Download the Docker image (or build it yourself using the Dockerfile in the repository) and run the container forwarding port 8888 to JupyterLab using this commands inside your command line terminal window:

```
docker pull mobigroup/pygmtsar

docker run -dp 8888:8888 --name pygmtsar docker.io/mobigroup/pygmtsar

docker logs pygmtsar
```

See the output for the JupyterLab link and copy and past it into your web browser address line. Also, the donwloaded Docker image can be started in Docker Desktop app - press "RUN" button and define the container name and the port in the opened dialog window (see "Optional settings" for the port number input field) and click on the newly created container to launch it and see the output log with the clickable link.

## Tutorials: Live Examples on Google Colab

Click on the examples below to run the processing in your own browser without any software installation. That's like to magic and it works.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/12LJqlZNBUmvLlRl98rRFCbKveVPg9Ami?usp=sharing) **ASF Downloading 2017 Iran–Iraq Earthquake vs GMTSAR GAMMA SNAP Co-Seismic Interferogram** The notebook **downloads Sentinel-1 Scenes from Alaska Satellite Facility (ASF)** and **compares the results to GMTSAR, SNAP and GAMMA Software**. Note: replace the scene names to produce an **interferogram** and **LOS displacement** for your area of interest.

<img src="https://user-images.githubusercontent.com/7342379/177748605-788889e5-9afd-44d8-bc3c-dc6efe920ea0.png" width="50%">

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1PyYcxvuyzhh-g4NQEbKjcfTDQhREZInn?usp=sharing) **Live Example S1A_2016_Kumamoto Earthquake_Co-Seismic Interferogram vs ESA Sentinel 1 Toolbox on Alaska Satellite Facility**. This is a single subswath processing with **landmask** applied to **interferogram**, **unwapped phase**, and **LOS, east-west, vertical displacement** results.

<img src="https://user-images.githubusercontent.com/7342379/183805898-d7c1ad76-822e-428e-9259-f19cc9e7540e.jpg" width="50%">

<img src="https://user-images.githubusercontent.com/7342379/183816622-1dacce7e-6a2f-46b9-8e67-d701f55bdd30.png" width="50%">

<img src="https://user-images.githubusercontent.com/7342379/183649417-7fcb7f3f-8c8d-45e8-a2c9-9293498ebada.png" width="50%">

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1ZTPV4HY-UoLvDYVx0UGh_Z3B12scSh9E?usp=sharing) **Live Example S1AB 2021 Crete Earthquake Co-Seismic Interferogram vs Centre of EO Research & Satellite Remote Sensing, Greece Report** This is a single **cropped subswath** processing with **landmask** applied to **interferogram**, **unwapped phase**, and **LOS, east-west, vertical displacement** results.

<img src="https://user-images.githubusercontent.com/7342379/177004287-cdd4351c-0834-42ae-8e46-9da5e8b124bf.jpg" width="50%">

<img src="https://user-images.githubusercontent.com/7342379/183645260-f8529ff3-b014-499e-ba2f-ebea4937b2c2.png" width="50%">

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1sljxm2jAMGXynq4EYam6Siz8OLcPLN0h?usp=sharing) **GMTSAR example dataset S1A_Stack_CPGF_T173** This example illustrates **SBAS** and **PSI** analyses and **detrending** approach to remove **atmospheric noise** to produce much better results.

<img src="https://user-images.githubusercontent.com/7342379/135814732-aa0eb142-ae54-4a57-b271-c33b5174a28e.png" width="50%">

<img src="https://user-images.githubusercontent.com/7342379/189961167-bf3901e5-417c-41ce-a5ca-d1c74c239a04.png" width="50%">

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1ZBVwlkiXMhSDS96oojpWrzTyRFIxv8Rp?usp=sharing) **ASF Downloading 2020 Ardabil, Iran Earthquake Co-Seismic Interferogram and LOS Displacement** The notebook **downloads Sentinel-1 Scenes from Alaska Satellite Facility (ASF)** to **crop the area** and **merge subswaths** and **detrend** results. Note: replace the scene names to produce an interferogram for your area of interest.

<img src="https://user-images.githubusercontent.com/7342379/194813466-fc4734a3-770d-4d6e-8012-91a4e5d781ba.png" width="50%">

<img src="https://user-images.githubusercontent.com/7342379/190451656-386d6cb8-f536-447c-8274-71d4f0435408.png" width="50%">

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/17D53uZu3XcEoWz5T73D__t9Ampzn5l3J?usp=sharing) **ASF Downloading 2023-02-06 Türkiye Earthquakes Co-Seismic Interferogram and LOS Displacement Projections** The notebook **downloads Sentinel-1 Scenes from Alaska Satellite Facility (ASF)** to **stitch 2 scenes** and **merge subswaths** and **detrend** results. Note: replace the scene names to produce an interferogram for your area of interest.

<img width="50%" src="https://user-images.githubusercontent.com/7342379/223332314-afd00f9d-0691-4d21-8be5-4c2ac96f8f3c.png">

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1H652deK0W7nujEky9j9K20729vywntuD?usp=sharing) **ASF Downloading 2023-02-06 Türkiye Earthquakes Co-Seismic Interferogram and LOS Displacement Projections**  The notebook **downloads Sentinel-1 Scenes from Alaska Satellite Facility (ASF)** to **stitch 3 scenes** and **merge subswaths** and **detrend** results. Here are some tricks used to process the large amount of data on Google Colab. Note: replace the scene names to produce an interferogram for your area of interest.

<img width="50%" src="https://user-images.githubusercontent.com/7342379/223332721-2dab4ef5-713a-4bc9-8f6b-1a968e481561.png">

<img width="50%" src="https://user-images.githubusercontent.com/7342379/223333303-91f81035-8ba9-4637-b257-ccbeb3502e95.png">

## Reference

PyGMTSAR defines a high-level Python SBAS class for an interferometry processing. Below listed the SBAS user-friendly functions grouped by common tasks in order of the full pipeline processing steps. For simple pipelines some of the steps can be missed.

The detailed documentation available interactively in docstring format, use help(function_name) to read it in your code editor.

### SBAS class functions

Firstly, import the class from PyGMTSAR Python package as

```
from pygmtsar import SBAS
```

#### Initialisation (PyGMTSAR original)

```
SBAS: Initialise SBAS object using scenes directory and some filters.
set_master: Set master scene for SBAS object (first scene is used by default).
```

#### Manage orbits (PyGMTSAR original)

```
download_orbits: Download missed orbits for all the SBAS scenes.
```

#### List scenes and orbits (PyGMTSAR original)

```
to_dataframe: Return Pandas Dataframe for all SBAS scenes.
```

#### Manage DEM (PyGMTSAR original)

```
set_dem: Set existing WGS84 DEM file in geographic coordinates for SBAS object.
download_dem: Download missed DEM file in geographic coordinates and convert it for WGS84 ellipsoid.
get_dem: Return WGS84 DEM in geographic coordinates as just one or list of Xarray Dataarray.
```

#### Manage topography in radar coordinates (combined GMTSAR wrapper and PyGMTSAR original)

```
topo_ra_parallel: Build WGS84 DEM in radar coordinates for interferogram processing.
get_topo_ra: Return WGS84 DEM in radar coordinates as one or list of Xarray Dataarray depending of the subswaths count.

```

#### Framing (combined GMTSAR wrapper and PyGMTSAR original)

```
set_pins: Set pins for scene framing. Pins are lon/lat points to crop extra bursts.
get_pins: Return pins list. Pins are lon/lat points.
reframe_parallel: Reorder bursts from sequential scenes to cover the full orbit area between pins only.
```

#### Build Interferogram (GMTSAR wrapper)

```
stack_parallel: Stack and align scenes.
intf_parallel: Build interferograms for all the subswaths separately.
pixel_decimator: Return function for pixel decimation to the specified output resolution.
```

#### Merging (GMTSAR wrapper)

```
merge_parallel: Merge all the separate subswath interferograms into one large interferogram.
```

#### Manage landmask (PyGMTSAR original)

```
set_landmask: Set existing land mask file name for SBAS object.
download_landmask: Download missed land mask file in geographic coordinates.
get_landmask: Return WGS84 land mask in geographic coordinates as Xarray Dataarray.
```

#### Unwrapping (SNAPHU wrapper)

```
snaphu_config: Return SNAPHU text configuration to use it for unwrapping.
unwrap_parallel: parallel phase unwrapping using SNAPHU.
```

#### Detrending (PyGMTSAR original)

```
detrend_parallel: Detrend and save to files a set of unwrapped interferograms combining topography and linear components removal and Gaussian filtering.
detrend: Detrend and return output for a single unwrapped interferogram combining topography and linear components removal and Gaussian filtering.

```

#### SBAS (GMTSAR wrapper)

```
baseline_pairs: Generate SBAS baseline pairs.
sbas: SBAS unwrapped interferograms timeseries analysis, see GMTSAR documentation for the details.
```

#### SBAS (PyGMTSAR original)

```
sbas_parallel: SBAS unwrapped and detrended interferograms timeseries analysis using pixel-wise correlation-weighted least squares approach.
```

#### Look vectors and Incidence angle (PyGMTSAR original)

```
sat_look_parallel: Build and save to file satellite look vectors in geographic coordinates.
get_sat_look: Return satellite look vectors in geographic coordinates as just one or list of Xarray Datasets.
incidence_angle: Compute incidence angle grid in geographic coordinates. 
```

#### Displacements (PyGMTSAR original)

```
los_displacement_mm: Compute LOS displacement in millimeters.
vertical_displacement_mm: Compute vertical displacement in millimeters in geographic coordinates.
eastwest_displacement_mm: Compute East-West displacement in millimeters. 
```

#### Data output (PyGMTSAR original)

```
open_grids: Lazy open PyGMTSAR produced NetCDF grids as 3D data cube and apply a set of functions on it.
```

#### Geocoding (PyGMTSAR original)

```
geocode_parallel: Build and save to files direct and inverse geocoding matrices.
intf_ra2ll: Geocoding function based on interferogram geocode matrix.
intf_ll2ra: Inverse geocoding function based on interferogram inverse geocode matrix.
```

#### Backup and restore (PyGMTSAR original)
```
dump: Dump SBAS object state to pickle file (SBAS.pickle in the processing directory by default).
restore: Restore SBAS object state from pickle file (SBAS.pickle in the processing directory by default).
backup: Backup framed SBAS scenes, orbits, DEM and landmask files to build a minimal reproducible dataset.
```

#### Helpers (PyGMTSAR original)

```
pixel_size: Compute ground pixel size in meters for the default processing grid or the defined one.
nearest_grid: Nearest neighbor interpolation.
cropna: Crop raster valid extent removing all rows and columns containing NODATA values only.
as_geo: Add geospatial attributes (CRS and spatial dimentions) to allow RioXarray raster operations.
```


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
