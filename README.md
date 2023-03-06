[![MacOS tests](https://github.com/mobigroup/gmtsar/actions/workflows/macos.yml/badge.svg)](https://github.com/mobigroup/gmtsar/actions/workflows/macos.yml)
[![Ubuntu tests](https://github.com/mobigroup/gmtsar/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/mobigroup/gmtsar/actions/workflows/ubuntu.yml)
[![PyPI tests](https://github.com/mobigroup/gmtsar/actions/workflows/pypi.yml/badge.svg)](https://github.com/mobigroup/gmtsar/actions/workflows/pypi.yml)
[![Available on pypi](https://img.shields.io/pypi/v/pygmtsar.svg)](https://pypi.python.org/pypi/pygmtsar/)
[![Docker](https://badgen.net/badge/icon/docker?icon=docker&label)](https://hub.docker.com/r/mobigroup/pygmtsar)

## PyGMTSAR (Python GMTSAR) - Sentinel-1 Satellite Interferometry For Everyone

Hint: You can sponsor PyGMTSAR software development on [Patreon](https://www.patreon.com/pechnikov).

<img src="https://user-images.githubusercontent.com/7342379/194891967-be2b56b5-c30c-4040-8ef8-39b448ce2390.jpg" width="40%" />

This repository based on forked original GMTSAR and extended by my patches to binary tools and Python library PyGMTSAR. I commit my changes to binary tools to GMTSAR upstream and so it's possible to use the original GMTSAR master branch installation plus PyGMTSAR Python package via PIP. The project documentation including installation instructions available by the link: https://mobigroup.github.io/gmtsar/

The goal of the project is easy and fast satellite interferometry (InSAR) processing for Sentinel-1 radar scenes everywhere as on localhost as on cloud environments like to Google Cloud VM and AI Notebooks and Amazon EC2 and on free of charge cloud environment Google Colab and in Docker images. GMTSAR binary command line tools are used under the hood but all GMTSAR scripts and GMT command replaced by Python code using modern and robust algorithms.

## PyGMTSAR Docker Images for Everyone on DockerHub

<img src="https://user-images.githubusercontent.com/7342379/203853391-b0dd50e5-3b07-4655-b5f4-c08109f23ffd.png" width="50%">

See Docker basic image for merged subswaths and cropped scenes and SBAS time series processing on [DockerHub PyGMTSAR for Everyone](https://hub.docker.com/r/mobigroup/pygmtsar) This image is the right choice to start and perform lots of common interferometry tasks. 

## PyGMTSAR Docker Images for Experts on DockerHub

<img src="https://user-images.githubusercontent.com/7342379/205281391-682816e5-8f1e-44bd-b7b6-478df7453bb1.png" width="50%">

See Docker images for multiple stitched scenes and long SBAS time series processing on [DockerHub PyGMTSAR for Experts](https://hub.docker.com/r/mobigroup/pygmtsar-large)

## PyGMTSAR Live Examples on Google Colab

The notebooks are interactive examples available directly in your web browser. All the steps automated including the software installation on Google Colab cloud host and downloading of Sentinel-1 orbit files, SRTM DEM (and its conversion to ellispoidal heights using EGM96 model), a landmask (to mask low-coherence water surfaces), Sentinel-1 SLC scenes from Alaska Satellite Facility (ASF) datastore and, of course, the complete interferometry processing and the results mapping.

### Notebooks to Compare Results to GMTSAR, SNAP and GAMMA Software

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

### More Complex Notebooks Still Available on Google Colab

The notebooks processing more than a single subswath or scene. It's possible on Google Colab limited resources using prepared datasets produced by PyGMTSAR "backup" command described in the notebooks.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1ZBVwlkiXMhSDS96oojpWrzTyRFIxv8Rp?usp=sharing) **ASF Downloading 2020 Ardabil, Iran Earthquake Co-Seismic Interferogram and LOS Displacement** The notebook **downloads Sentinel-1 Scenes from Alaska Satellite Facility (ASF)** to **crop the area** and **merge subswaths** and **detrend** results. Note: replace the scene names to produce an interferogram for your area of interest.

<img src="https://user-images.githubusercontent.com/7342379/194813466-fc4734a3-770d-4d6e-8012-91a4e5d781ba.png" width="50%">

<img src="https://user-images.githubusercontent.com/7342379/190451656-386d6cb8-f536-447c-8274-71d4f0435408.png" width="50%">

### Long Timeseries Analysis is not available on Google Colab 

See a separate GitHub repository for Yamchi Dam area dynamic model [YamchiDam](https://github.com/mobigroup/YamchiDam) Software tools [PyGMTSAR](https://github.com/mobigroup/gmtsar) [N-Cube ParaView plugin for 3D/4D GIS Data Visualization](https://github.com/mobigroup/ParaView-plugins) combined together for the 4D analysis and visualization:

<img src="https://user-images.githubusercontent.com/7342379/144747743-a24d72ec-8875-4272-91f9-ec1f937bb798.gif" width="50%">

## Learn more

- Documentation: https://github.com/mobigroup/gmtsar

- Issue tracker: https://github.com/mobigroup/gmtsar/issues

- Source code: https://github.com/mobigroup/gmtsar

- Docker Images: https://hub.docker.com/repository/docker/mobigroup/pygmtsar

- PyPI Python library: https://pypi.org/project/pygmtsar/

## About me

I have STEM master’s degree in radio physics and in 2004 I have got the first prize of the All-Russian Physics competition for significant results in forward and inverse modeling for non-linear optics and holography, also applicable for modeling of Gravity, Magnetic, and Thermal fields and satellite interferometry processing. And I’m data scientist and software developer with 20 year’s experience in science and industrial development. I had been working on government contracts and universities projects and on projects for LG Corp, Google Inc, etc. By the way, I left Russia many years ago and I work remotely for about 20 years.

[Geological models on YouTube channel](https://www.youtube.com/channel/UCSEeXKAn9f_bDiTjT6l87Lg)

[Augmented Reality (AR) Geological Models](https://mobigroup.github.io/ParaView-Blender-AR/)

[GitHub repositories](https://github.com/mobigroup)

[DockerHub repositories](https://hub.docker.com/u/mobigroup)

[English posts and publications on LinkedIn](https://www.linkedin.com/in/alexey-pechnikov/)

[English publications on Medium](https://medium.com/@pechnikov)

[Russian publications on Habr](https://habr.com/ru/users/N-Cube/posts/)

@ Alexey Pechnikov, 2023
