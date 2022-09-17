## PyGMTSAR (Python GMTSAR) - Easy and Fast Satellite Interferometry For Everyone

**Note**: this repository based on forked original GMTSAR plus my patches, shell scripts for SBAS pipeline on cloud hosts and Python codes. I commit my changes to GMTSAR and so it's possible to install original GMTSAR master branch following the upstream installation instructions and PyGMTSAR Python package via PIP. To simplify the process the Live example notebooks for Google Colab below automatically install all the software. MacOS and Linux Debian (Ubuntu) are my prefered OSes and I support both of them, see Debian cloud hosts init scripts in [gmtsar/sh](https://github.com/mobigroup/gmtsar/tree/master/gmtsar/sh) directory.

The goal of the project is easy and fast satellite interferometry (InSAR) processing everywhere as on localhost as on cloud environments like to Google Cloud VM and AI Notebooks and Amazon EC2 and on free of charge cloud environment Google Colab. GMTSAR binary command line tools are used under the hood but all GMTSAR scripts and GMT command replaced by Python code using modern and robust algorithms. By my opinion, GMTSAR is great project and has only one big problem - it is not user friendly. I build rich Python API instead of GMTSAR CSH scripts and GMT toolkit calls and by this way we have 3-10 times faster processing and more accurate results. How is it possible? The key is parallelization and algorithms - all the scenes and subswaths processing using all available processor cores and SNAPHU unwrapping configurations allow to use all the cores as for a single interferogram unwrapping as for multiple ones. Also, slow GMTSAR direct and inverse geocoding is replaced by fast matrix geocoding transformation and DEM processing using smoth cubic spline with smooth derivative instead of inaccurate and slow GMT tension surfaces (which derivative is not continuos on any scale so any numeric operations on it produces big artifacts). By this way, PyGMTSAR = GMTSAR - GMT + Python. And a couple of words about the really popular question - would be Python fast enough to use it for intensive numerical computations? The answer is exactly "yes" because basic scientific Python libraries coded using compiled C and Fortran languages and use as code parallelization as vectorization. Actually, GMTSAR binaries don't use any of these and Python processing can be significantly faster. Be careful Jupyter notebooks use lot of RAM for optional interactive maps (by default enabled on MacOS only where 16 GB RAM is enough). To limit RAM consuming use PyGMTSAR in command-line scripts still having progress indicators and other advantages.

PyGMTSAR automatically downloads Sentinel-1 orbit files and SRTM DEM (and converts it to ellispoidal heights using EGM96 model) and even Sentinel-1 raw scenes using Alaska Satellite Facility (ASF) datastore.

The project documentation including installation instructions available by the link: https://mobigroup.github.io/gmtsar/

### Live Examples on Google Colab - just click on the links to run the processing in your own browser without any software installation

The notebooks include as common static plots as additional interactive plots which can't be processing on Google Colab. Run the notebooks on MacOS or Linux Debian/Ubuntu to visualize all the content. Note: I work locally on Apple Silicon and Intel MacOS and use Debian and Ubuntu Linux on cloud servers so PyGMTSAR is well tested on all of them. 

#### Simple Notebooks on Google Colab to Compare Results to GMTSAR, SNAP and GAMMA Software

* [ASF Downloading 2017 Iran–Iraq Earthquake vs GMTSAR GAMMA SNAP Co-Seismic Interferogram](https://colab.research.google.com/drive/12LJqlZNBUmvLlRl98rRFCbKveVPg9Ami?usp=sharing). The notebook **downloads Sentinel-1 Scenes from Alaska Satellite Facility (ASF)** and **compares the results to GMTSAR, SNAP and GAMMA Software**. Note: replace the scene names to produce an **interferogram** and **LOS displacement** for your area of interest.

<img src="https://user-images.githubusercontent.com/7342379/177748605-788889e5-9afd-44d8-bc3c-dc6efe920ea0.png" width="50%">

* [Live Example S1A_2016_Kumamoto Earthquake_Co-Seismic Interferogram vs ESA Sentinel 1 Toolbox on Alaska Satellite Facility](https://colab.research.google.com/drive/1PyYcxvuyzhh-g4NQEbKjcfTDQhREZInn?usp=sharing). This is a single subswath processing with **landmask** applied to **interferogram**, **unwapped phase**, and **LOS, east-west, vertical displacement** results.

<img src="https://user-images.githubusercontent.com/7342379/183805898-d7c1ad76-822e-428e-9259-f19cc9e7540e.jpg" width="50%">

<img src="https://user-images.githubusercontent.com/7342379/183816622-1dacce7e-6a2f-46b9-8e67-d701f55bdd30.png" width="50%">

<img src="https://user-images.githubusercontent.com/7342379/183649417-7fcb7f3f-8c8d-45e8-a2c9-9293498ebada.png" width="50%">

* [Live Example S1AB 2021 Crete Earthquake Co-Seismic Interferogram vs Centre of EO Research & Satellite Remote Sensing, Greece Report](https://colab.research.google.com/drive/1ZTPV4HY-UoLvDYVx0UGh_Z3B12scSh9E?usp=sharing) This is a single **cropped subswath** processing with **landmask** applied to **interferogram**, **unwapped phase**, and **LOS, east-west, vertical displacement** results.

<img src="https://user-images.githubusercontent.com/7342379/177004287-cdd4351c-0834-42ae-8e46-9da5e8b124bf.jpg" width="50%">

<img src="https://user-images.githubusercontent.com/7342379/183645260-f8529ff3-b014-499e-ba2f-ebea4937b2c2.png" width="50%">

* [GMTSAR example dataset S1A_Stack_CPGF_T173](https://colab.research.google.com/drive/1sljxm2jAMGXynq4EYam6Siz8OLcPLN0h?usp=sharing) This example illustrates **SBAS** and **PSI** analyses and **detrending** approach to remove **atmospheric noise** to produce much better results.

<img src="https://user-images.githubusercontent.com/7342379/135814732-aa0eb142-ae54-4a57-b271-c33b5174a28e.png" width="50%">

<img src="https://user-images.githubusercontent.com/7342379/189961167-bf3901e5-417c-41ce-a5ca-d1c74c239a04.png" width="50%">

#### More Complex Notebooks Still Available on Google Colab

The notebooks processing more than a single subswath or scene. It's possible on Google Colab limited resources using prepared datasets produced by PyGMTSAR "backup" command described in the notebooks.

* [ASF Downloading 2020 Ardabil, Iran Earthquake Co-Seismic Interferogram and LOS Displacement](https://colab.research.google.com/drive/1ZBVwlkiXMhSDS96oojpWrzTyRFIxv8Rp?usp=sharing). The notebook **downloads Sentinel-1 Scenes from Alaska Satellite Facility (ASF)** to **crop the area** and **merge subswaths** and **detrend** results. Note: replace the scene names to produce an interferogram for your area of interest.

#### Long Timeseries Analysis is not available on Google Colab 

See a separate GitHub repository for Yamchi Dam area dynamic model [YamchiDam](https://github.com/mobigroup/YamchiDam) Here two of my software tools [PyGMTSAR](https://github.com/mobigroup/gmtsar) [N-Cube ParaView plugin for 3D/4D GIS Data Visualization](https://github.com/mobigroup/ParaView-plugins) are combined together for 4D analysis and visualization:

<img src="https://user-images.githubusercontent.com/7342379/144747743-a24d72ec-8875-4272-91f9-ec1f937bb798.gif" width="50%">

### About me

I have STEM master's degree in radio physics and in 2004 I was awarded first prize of the All-Russian Physics competition for significant results in Inverse modeling for non-linear optics and holography, also applicable for Inverse Modeling of Gravity, Magnetic, and Thermal fields. To create laser-induced holograms in non-linear optical composites I worked on interferograms numerical modeling and development of satellite interferometry processing software is very close task and so I build PyGMTSAR. Also, that's the related to inverse modeling of potensial fields like to gravity and I build Geomed3D geophisical modeling software too. In addition to my fundamental science knowledge, I’m world class data scientist and software developer with 20 years experience in science and industrial development. I have worked on government contracts and universities projects and on projects for LG Corp, Google Inc, etc. You are able to find some of my software and results on LinkedIn and GitHub and Upwork, see the links below. By the way, I left Russia many years ago and I work remotely for about 20 years.

To order some research, development and support see my profile on freelance platform [Upwork](https://www.upwork.com/freelancers/~01e65e8e7221758623) And of cource you are able to use my Open Source software for you scientific research and geological exploration projects and beyond.

 [Geological models on YouTube channel](https://www.youtube.com/channel/UCSEeXKAn9f_bDiTjT6l87Lg)

 [Augmented Reality (AR) Geological Models](https://mobigroup.github.io/ParaView-Blender-AR/)

 [GitHub repositories](https://github.com/mobigroup)

 [English posts and articles on LinkedIn](https://www.linkedin.com/in/alexey-pechnikov/)

[Russian articles on Habr](https://habr.com/ru/users/N-Cube/posts/)

@ Alexey Pechnikov, 2022
