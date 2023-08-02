## Announce: ebook "PyGMTSAR: Sentinel-1 Python InSAR: An Introduction" is now available on [Amazon, Apple, Kobo and many other bookstores](https://books2read.com/b/PyGMTSAR-introduction).

[![MacOS tests](https://github.com/mobigroup/gmtsar/actions/workflows/macos.yml/badge.svg)](https://github.com/mobigroup/gmtsar/actions/workflows/macos.yml)
[![Ubuntu tests](https://github.com/mobigroup/gmtsar/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/mobigroup/gmtsar/actions/workflows/ubuntu.yml)
[![PyPI tests](https://github.com/mobigroup/gmtsar/actions/workflows/pypi.yml/badge.svg)](https://github.com/mobigroup/gmtsar/actions/workflows/pypi.yml)
[![Available on pypi](https://img.shields.io/pypi/v/pygmtsar.svg)](https://pypi.python.org/pypi/pygmtsar/)
[![Docker](https://badgen.net/badge/icon/docker?icon=docker&label)](https://hub.docker.com/r/mobigroup/pygmtsar)
[![DOI](https://zenodo.org/badge/398018212.svg)](https://zenodo.org/badge/latestdoi/398018212)

## PyGMTSAR (Python GMTSAR) - Sentinel-1 Satellite Interferometry for Everyone

<img src="https://user-images.githubusercontent.com/7342379/194891967-be2b56b5-c30c-4040-8ef8-39b448ce2390.jpg" width="40%" />

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

## PyGMTSAR Docker Images for Everyone on DockerHub

The PyGMTSAR project provides Docker images that are readily available on DockerHub. These Docker images are designed to facilitate the use of PyGMTSAR for various interferometry tasks, including merged subswaths, cropped scenes, and SBAS time series processing.

You can find the PyGMTSAR Docker images on [DockerHub PyGMTSAR for Everyone](https://hub.docker.com/r/mobigroup/pygmtsar). These images provide a convenient and easily accessible platform for running PyGMTSAR without the need for complex software installation or configuration. Simply download the Docker image and start using PyGMTSAR right away.

The PyGMTSAR Docker images are a great choice for getting started with satellite interferometry tasks and performing common processing operations. They offer a seamless and hassle-free experience, allowing you to focus on your data analysis and insights.

![Docker Image](https://user-images.githubusercontent.com/7342379/203853391-b0dd50e5-3b07-4655-b5f4-c08109f23ffd.png)

By utilizing the PyGMTSAR Docker images, you can leverage the power of PyGMTSAR in a convenient and easily deployable manner, making it accessible to everyone interested in satellite interferometry.

## PyGMTSAR Docker Images for Experts on DockerHub

The PyGMTSAR project also provides Docker images tailored for experts and advanced users who require more advanced functionalities and processing capabilities. These Docker images are available on [DockerHub PyGMTSAR for Experts](https://hub.docker.com/r/mobigroup/pygmtsar-large).

The Docker images for experts are designed to handle multiple stitched scenes and support long SBAS time series processing. These images offer enhanced processing power and scalability, allowing users to perform more complex and resource-intensive interferometry tasks.

By utilizing the PyGMTSAR Docker images for experts, users can leverage the full potential of PyGMTSAR for their advanced processing needs. The images provide a comprehensive platform with the necessary tools and libraries to tackle sophisticated satellite interferometry tasks efficiently.

![Docker Image](https://raw.githubusercontent.com/mobigroup/articles/main/205281391-682816e5-8f1e-44bd-b7b6-478df7453bb1.png)

If you are an expert user and require extensive processing capabilities and scalability, the PyGMTSAR Docker images for experts are the ideal choice. They enable you to harness the full power of PyGMTSAR for your advanced satellite interferometry projects.

## PyGMTSAR Live Examples on Google Colab

These notebooks offer an automated environment where all the necessary software installations, data downloads, and processing steps are handled seamlessly. By opening the notebooks in your web browser, you can directly interact with the examples and explore the functionalities of PyGMTSAR.

The automated process includes the installation of the required software on Google Colab's cloud host, such as the downloading of Sentinel-1 orbit files, SRTM DEM data (which is also converted to ellipsoidal heights using the EGM96 model), and a landmask to mask low-coherence water surfaces. Additionally, the notebooks provide access to Sentinel-1 SLC scenes from the Alaska Satellite Facility (ASF) datastore, enabling you to perform complete interferometry processing and visualize the results.

With these Live Examples, you can gain hands-on experience with satellite interferometry using PyGMTSAR without the need for local installations or configurations. The step-by-step instructions and interactive nature of the notebooks make it easy to follow along and understand the various processing steps involved.

### Notebooks to Compare Results with GMTSAR, SNAP, and GAMMA Software

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/12LJqlZNBUmvLlRl98rRFCbKveVPg9Ami?usp=sharing) **ASF Downloading 2017 Iran–Iraq Earthquake vs GMTSAR GAMMA SNAP Co-Seismic Interferogram** This notebook **downloads Sentinel-1 scenes from the Alaska Satellite Facility (ASF)** to produce **interferogram** and **coherence** map for a **single cropped subswath** and **export NetCDF rasters** (compatible to QGIS, GDAL, and other GIS software) and **compare the results to GMTSAR, SNAP, and GAMMA software**. Note: Replace the scene names to generate an **interferogram** and **LOS displacement** for your area of interest.

<img src="https://user-images.githubusercontent.com/7342379/177748605-788889e5-9afd-44d8-bc3c-dc6efe920ea0.png" width="50%">

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1PyYcxvuyzhh-g4NQEbKjcfTDQhREZInn?usp=sharing) **Live Example: S1A_2016_Kumamoto Earthquake_Co-Seismic Interferogram vs ESA Sentinel 1 Toolbox on Alaska Satellite Facility**. This example demonstrates **single subswath** processing with a **landmask** applied to the **interferogram**, **unwrapped phase**, and **LOS, east-west, and vertical displacement** results.

<img src="https://user-images.githubusercontent.com/7342379/183805898-d7c1ad76-822e-428e-9259-f19cc9e7540e.jpg" width="50%">

<img src="https://user-images.githubusercontent.com/7342379/183816622-1dacce7e-6a2f-46b9-8e67-d701f55bdd30.png" width="50%">

<img src="https://user-images.githubusercontent.com/7342379/183649417-7fcb7f3f-8c8d-45e8-a2c9-9293498ebada.png" width="50%">

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1ZTPV4HY-UoLvDYVx0UGh_Z3B12scSh9E?usp=sharing) Live Example: S1AB 2021 Crete Earthquake Co-Seismic Interferogram vs Centre of EO Research & Satellite Remote Sensing, Greece Report. This example features the processing of a **single cropped subswath** with a **landmask** applied to the **interferogram**, **unwrapped phase**, and **LOS**, **east-west**, and **vertical displacement** results. The output compared to the Centre of EO Research & Satellite Remote Sensing, Greece Report.

<img src="https://user-images.githubusercontent.com/7342379/177004287-cdd4351c-0834-42ae-8e46-9da5e8b124bf.jpg" width="50%">

<img src="https://user-images.githubusercontent.com/7342379/183645260-f8529ff3-b014-499e-ba2f-ebea4937b2c2.png" width="50%">

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1sljxm2jAMGXynq4EYam6Siz8OLcPLN0h?usp=sharing) This example demonstrates the **SBAS analysis** and **detrending** approach used to remove atmospheric noise, resulting in significantly improved outcomes.

<img src="https://user-images.githubusercontent.com/7342379/135814732-aa0eb142-ae54-4a57-b271-c33b5174a28e.png" width="50%">

<img src="https://user-images.githubusercontent.com/7342379/189961167-bf3901e5-417c-41ce-a5ca-d1c74c239a04.png" width="50%">

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1xhVedrIvNS66jGKgS30Dgqy0S31uJ8gm?usp=sharing) This example demonstrates the **SBAS analysis** and **detrending** approach to remove atmospheric noise for improved results. Additionally, an **OpenStreetMap road** mask is used to **unwrap** and analyze only the roads.

<img src="https://user-images.githubusercontent.com/7342379/233364597-db66cf85-a748-4188-8cd1-f9fa739da228.png" width="50%">

<img src="https://user-images.githubusercontent.com/7342379/233364568-eb42f7a6-6685-46b6-8f8b-bae9de1f6a0b.png" width="50%">

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1eJzRtpMp7_7QsJ5jyzdb-fiUgo8F09Xm?usp=sharing) Pico do Fogo Volcano Eruption on Cape Verde's Fogo Island, 2014. This example features the processing of a **single cropped subswath** with a **landmask** applied to the **interferogram**, **unwrapped phase**, and **LOS**, **east-west**, and **vertical displacement** maps.

<img src="https://github.com/mobigroup/gmtsar/assets/7342379/dcd412ac-ae8b-4dc8-b9bb-13c9877bbe38" width="50%">

### More Complex Notebooks Still Available on Google Colab

The notebooks processing more than a single subswath or scene. It's possible on Google Colab limited resources using prepared datasets produced by PyGMTSAR "backup" command described in the notebooks.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1ZBVwlkiXMhSDS96oojpWrzTyRFIxv8Rp?usp=sharing) **ASF Downloading 2020 Ardabil, Iran Earthquake Co-Seismic Interferogram and LOS Displacement** The notebook **downloads Sentinel-1 Scenes from Alaska Satellite Facility (ASF)** to **crop the area** and **merge subswaths** and **detrend** phase and **export NetCDF rasters** (compatible to QGIS, GDAL, and other GIS software). Note: replace the scene names to produce an interferogram for your area of interest.

<img src="https://user-images.githubusercontent.com/7342379/194813466-fc4734a3-770d-4d6e-8012-91a4e5d781ba.png" width="50%">

<img src="https://user-images.githubusercontent.com/7342379/190451656-386d6cb8-f536-447c-8274-71d4f0435408.png" width="50%">

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/17D53uZu3XcEoWz5T73D__t9Ampzn5l3J?usp=sharing) **ASF Downloading 2023-02-06 Türkiye Earthquakes Co-Seismic Interferogram and LOS Displacement Projections** The notebook **downloads Sentinel-1 Scenes from Alaska Satellite Facility (ASF)** to **stitch 2 scenes** and **merge subswaths** and **detrend** phase and **lazy export NetCDF rasters** (compatible to QGIS, GDAL, and other GIS software). Note: replace the scene names to produce an interferogram for your area of interest.

<img width="50%" src="https://user-images.githubusercontent.com/7342379/223332314-afd00f9d-0691-4d21-8be5-4c2ac96f8f3c.png">

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1H652deK0W7nujEky9j9K20729vywntuD?usp=sharing) **ASF Downloading 2023-02-06 Türkiye Earthquakes Co-Seismic Interferogram and LOS Displacement Projections**  The notebook **downloads Sentinel-1 Scenes from Alaska Satellite Facility (ASF)** to **stitch 3 scenes** and **merge subswaths** and **detrend** phase and **lazy export NetCDF rasters** (compatible to QGIS, GDAL, and other GIS software). Here are some tricks used to process the large amount of data on Google Colab. Note: replace the scene names to produce an interferogram for your area of interest.

<img width="50%" src="https://user-images.githubusercontent.com/7342379/223332721-2dab4ef5-713a-4bc9-8f6b-1a968e481561.png">

<img width="50%" src="https://user-images.githubusercontent.com/7342379/223333303-91f81035-8ba9-4637-b257-ccbeb3502e95.png">

### Preview Notebooks for Development PyGMTSAR Version

One of the most awaited features in PyGMTSAR (Python InSAR) is the combined analysis of Persistent Scatterers (PS or PSI) and the Small Baseline Subset (SBAS).

Each of the PS and SBAS techniques has unique advantages and drawbacks — with SBAS performing better in rural areas and PS in urban ones. My vision involves merging the benefits of both methods and mitigating their shortcomings in a unified PS-SBAS process. In the development version, PyGMTSAR offers persistent scatterer analysis and weighted interferogram processing. This emphasizes stable pixels, enhancing phase and coherence. This not only improves the accuracy of results but also simplifies SBAS analysis by maintaining high coherence, even in rural areas.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/18KGGO9jy_EE8RejgYhKgqyxryhrQdmbS?usp=sharing) The updated SBAS notebook provides an illustration of the PS-SBAS analysis.

<img width="50%" src="https://raw.githubusercontent.com/mobigroup/articles/main/avg_corr_15m_adi.40%25.png">

### Long Timeseries Analysis is not available on Google Colab 

For long timeseries analysis, you can explore the separate GitHub repository for the Yamchi Dam area dynamic model, [YamchiDam](https://github.com/mobigroup/YamchiDam). This repository combines the software tools [PyGMTSAR](https://github.com/mobigroup/gmtsar) and the [N-Cube ParaView plugin for 3D/4D GIS Data Visualization](https://github.com/mobigroup/ParaView-plugins) to provide comprehensive 4D analysis and visualization capabilities. By exploring these repositories, you can learn more about the tools and methods used in the Yamchi Dam area dynamic model.

<img src="https://user-images.githubusercontent.com/7342379/144747743-a24d72ec-8875-4272-91f9-ec1f937bb798.gif" width="50%">

## Documentation

With PyGMTSAR, you can leverage the ready-to-use interactive examples in Live Google Colab notebooks. These examples serve as a starting point for your own analysis and can be easily modified to suit your specific needs. The advantage of using PyGMTSAR in a Live Jupyter notebook environment is that you can immediately access the functionality without the hassle of software installation and configuration.

PyGMTSAR provides self-documented functions using Python docstrings. To access the complete documentation for a specific function, you can use the `help()` function in Jupyter notebook cells or in a Python editor. Simply pass the function name as an argument to `help()` and it will display the docstring, which contains detailed information about the function's usage, parameters, and return values.

By using `help()`, you can easily access the comprehensive documentation for each function in PyGMTSAR and gain a better understanding of its functionality and how to use it effectively in your code.

In addition to the `help()` function, you can also explore the PyGMTSAR sources and track any reported issues on the [PyGMTSAR GitHub](https://github.com/mobigroup/gmtsar) repository. Documentation is available on the [PyGMTSAR GitHub Pages](https://mobigroup.github.io/gmtsar/), providing further guidance on using PyGMTSAR and installing GMTSAR.

## Learn more

Here are some additional resources to learn more about PyGMTSAR:

- Documentation: [PyGMTSAR Documentation](https://github.com/mobigroup/gmtsar)
- Issue tracker: [PyGMTSAR Issue Tracker](https://github.com/mobigroup/gmtsar/issues)
- Source code: [PyGMTSAR GitHub Repository](https://github.com/mobigroup/gmtsar)
- Docker Images: [PyGMTSAR Docker Images](https://hub.docker.com/repository/docker/mobigroup/pygmtsar)
- PyPI Python library: [PyGMTSAR on PyPI](https://pypi.org/project/pygmtsar/)

These resources will provide you with detailed documentation, issue tracking, source code, Docker images, and the PyPI package for PyGMTSAR. They will help you explore, troubleshoot, and utilize the functionalities of PyGMTSAR effectively.

## About me

I hold a master's degree in radio physics with a focus on forward and inverse modeling for non-linear optics, holography, gravity, magnetic, and thermal fields, and satellite interferometry processing. In 2004, I was awarded the first prize in the All-Russian Physics competition for my significant contributions in these areas. With over 20 years of experience as a data scientist and software developer, I have worked on various projects for government contracts, universities, and companies like LG Corp and Google Inc.

I also have a YouTube channel where I share geological models and an Augmented Reality (AR) Geological Models project. You can find my GitHub repositories and DockerHub repositories for my open-source projects. I have published articles and posts in English on LinkedIn and Medium, as well as in Russian on Habr.

Feel free to explore my work and connect with me on these platforms.

[Geological models on YouTube channel](https://www.youtube.com/channel/UCSEeXKAn9f_bDiTjT6l87Lg)

[Augmented Reality (AR) Geological Models](https://mobigroup.github.io/ParaView-Blender-AR/)

[GitHub repositories](https://github.com/mobigroup)

[DockerHub repositories](https://hub.docker.com/u/mobigroup)

[English posts and publications on LinkedIn](https://www.linkedin.com/in/alexey-pechnikov/)

[English publications on Medium](https://medium.com/@pechnikov)

[Russian publications on Habr](https://habr.com/ru/users/N-Cube/posts/)

@ Alexey Pechnikov, 2023
