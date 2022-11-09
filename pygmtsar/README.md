[![MacOS tests](https://github.com/mobigroup/gmtsar/actions/workflows/macos.yml/badge.svg)](https://github.com/mobigroup/gmtsar/actions/workflows/macos.yml)
[![Ubuntu tests](https://github.com/mobigroup/gmtsar/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/mobigroup/gmtsar/actions/workflows/ubuntu.yml)
[![PyPI tests](https://github.com/mobigroup/gmtsar/actions/workflows/pypi.yml/badge.svg)](https://github.com/mobigroup/gmtsar/actions/workflows/pypi.yml)
[![Available on pypi](https://img.shields.io/pypi/v/pygmtsar.svg)](https://pypi.python.org/pypi/pygmtsar/)
[![Docker](https://badgen.net/badge/icon/docker?icon=docker&label)](https://hub.docker.com/r/mobigroup/pygmtsar)

**PyGMTSAR** (Python GMTSAR) is an open source project and Python package that provides Sentinel-1 Satellite Interferometry for everyone.

The goal of the project is easy and fast satellite interferometry (InSAR) processing for Sentinel-1 radar scenes everywhere as on localhost as on cloud environments like to Google Cloud VM and AI Notebooks and Amazon EC2 and on free of charge cloud environment Google Colab and in Docker images. GMTSAR binary command line tools are used under the hood but all GMTSAR scripts and GMT command replaced by Python code using modern and robust algorithms.

<img src="https://user-images.githubusercontent.com/7342379/195076831-5a1eece1-16d5-4a44-96f8-e34a69432c1b.jpg" width="40%" />

## Why PyGMTSAR?

PyGMTSAR combines powerful Python instrumentary for sophisticated multidementional processing (xarray library) and lazy calculations (dask library) plus parallel computing (dask and joblib libraries) to perform fast and interactive processing on huge datasets. And the best algorithms and numerical computation approaches applied for all the processing steps. There are progressbars and preview plots for the every step and that's easy to save intermediate results and continue work later on the same or other host. And (thanks to joblib library) that's safe to interrupt the execution at any time without memory leaks (common for dask-based solutions).

Thanks to all the powerful Python libraries and the best used algorithms PyGMTSAR is fast and its possible to complete SBAS analysis for 5 years on 800 interferograms in just one day on Apple Air or Apple iMac (8 cores and 16 GB RAM) using 2 TB Sentinel-1 SLC scenes. And PyGMTSAR is user-friendly providing functions to download the required satellite orbit files and DEM and so on. This combination of the human-readable and short code and powerful computing is the key to use PyGMTSAR everywhere from education and to research and more. 

### Live Examples in Docker image

Configure your Docker runtime (Preferences -> Resources tab for Docker Desktop) to use 2 CPU cores and 8 GB RAM or 4 CPU cores and 16 GB RAM and so on. Download the Docker image (or build it yourself using the Dockerfile in the repository) and run the container forwarding port 8888 to JupyterLab using this commands inside your command line terminal window:

```
docker pull mobigroup/pygmtsar

docker run -dp 8888:8888 --name pygmtsar docker.io/mobigroup/pygmtsar

docker logs pygmtsar
```

See the output for the JupyterLab link and copy and past it into your web browser address line. Also, the donwloaded Docker image can be started in Docker Desktop app - press "RUN" button and define the container name and the port in the opened dialog window (see "Optional settings" for the port number input field) and click on the newly created container to launch it and see the output log with the clickable link.

## Live Examples on Google Colab

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

## Learn more

- Documentation: https://github.com/mobigroup/gmtsar

- Issue tracker: https://github.com/mobigroup/gmtsar/issues

- Source code: https://github.com/mobigroup/gmtsar

- Docker Images: https://hub.docker.com/repository/docker/mobigroup/pygmtsar

- PyPI Python library: https://pypi.org/project/pygmtsar/
