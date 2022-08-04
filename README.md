## PyGMTSAR (Python GMTSAR) - Easy and Fast Satellite Interferometry For Everyone

The goal of the project is easy and fast satellite interferometry (InSAR) processing everywhere as on localhost as on cloud environments like to Google Cloud VM and AI Notebooks and Amazon EC2 and on free of charge cloud environment Google Colab. GMTSAR binary command line tools are used under the hood but all GMTSAR scripts and GMT command replaced by Python code using modern and robust algorithms. By my opinion, GMTSAR is great project and has only one big problem - it is not user friendly. I build rich Python API instead of GMTSAR CSH scripts and GMT toolkit calls and by this way we have 3-10 times faster processing and more accurate results. How is it possible? The key is parallelization and algorithms - all the scenes and subswaths processing using all available processor cores and SNAPHU unwrapping configurations allow to use all the cores as for a single interferogram unwrapping as for multiple ones. Also, slow GMTSAR geocoding is replaced by fast matrix geocoding transformation and DEM processing using smoth cubic spline with smooth derivative instead of inaccurate and very slow GMT tension surfaces (which derivative is not continuos on any scale so any numeric operations on it produces big artifacts). By this way, PyGMTSAR = GMTSAR - GMT + Python. And a couple of words about the really popular question - would be Python fast enough to use it for intensive calculation? The answer is exactly yes because basic scientific Python libraries based on compiled C and Fortran codes and using the both code parallelization and vectorization. Actually, GMTSAR binaries don't use any of these and Python processing can be significantly faster. I even don't speak about GMT binaries which are extremely slow due to inedffective iterative algorithms usage instead of hardware accelerated linear algebra operations. At the same time, PyGMTSAR uses more RAM sometimes although the most of RAM consumption requires for optional interactive maps and Jupyter Notebooks user interface. To limit RAM usage that's possible to use PyGMTSAR in command-line scripts still having progress indicators and other advantages. 

MacOS and Linux Debian are my prefered OSes and I support both of them.

PyGMTSAR automatically downloads Sentinel-1 orbit files and SRTM DEM (and converts it to ellispoidal heights using EGM96 model) and even Sentinel-1 raw scenes using Alaska Satellite Facility (ASF) datastore.

### Live Examples on Google Colab - just click on the links to run the processing in your own browser without any software installation

The notebooks include additional interactive plots which can't be processing on Google Colab. Run the notebooks on MacOS or Linux Debian/Ubuntu to visualize all the content.

#### Notebooks to Compare Results

* [ASF Downloading 2017 Iran–Iraq Earthquake vs GMTSAR GAMMA SNAP Co-Seismic Interferogram and LOS Displacement](https://colab.research.google.com/drive/184_x1hFyn9vCg03cH_q11POyb0bcjyZC?usp=sharing). The notebook downloads Sentinel-1 Scenes from Alaska Satellite Facility (ASF) and compares the results to SNAP, GMTSAR and GAMMA Software. Note: replace the scene names to produce an interferogram for your area of interest.

<img src="https://user-images.githubusercontent.com/7342379/177748605-788889e5-9afd-44d8-bc3c-dc6efe920ea0.png" width="40%">

* [Live Example S1A_2016_Kumamoto Earthquake_Co-Seismic Interferogram vs ESA Sentinel 1 Toolbox on Alaska Satellite Facility](https://colab.research.google.com/drive/1NjWXrMkFRB6E0-gViBLmmxprTqqRbDl5?usp=sharing)

<img src="https://user-images.githubusercontent.com/7342379/182755940-1a48f116-3673-4f88-939d-8b9521ccc949.jpg" width="50%">

* [Live Example S1AB 2021 Crete Earthquake Co-Seismic Interferogram vs Centre of EO Research & Satellite Remote Sensing, Greece Report](https://colab.research.google.com/drive/1oYHzNeH5mLuYFXcUxCGDrVBQYm9PUEfK?usp=sharing)

<img src="https://user-images.githubusercontent.com/7342379/177004287-cdd4351c-0834-42ae-8e46-9da5e8b124bf.jpg" width="40%">

* [GMTSRAR example dataset S1A_Stack_CPGF_T173](https://colab.research.google.com/drive/1S7qJzn-1CFv02yU9NE8ohyDli1JM11cl?usp=sharing)

<img src="https://user-images.githubusercontent.com/7342379/135814732-aa0eb142-ae54-4a57-b271-c33b5174a28e.png" width="50%">

<img src="https://user-images.githubusercontent.com/7342379/136558388-cffbcea5-e5a7-44d5-ba36-b51a7e0f10e9.png" width="50%">

#### Complex Processing

The notebooks processing more than a single subswath or scene. It's possible on Google Colab limited resources using prepared datasets produced by PyGMTSAR "backup" command described in the notebooks.

* [ASF Downloading 2020 Ardabil, Iran Earthquake Co-Seismic Interferogram and LOS Displacement](https://colab.research.google.com/drive/1gr6pHFu1_gRltjfSuGLGoYqgwDQHNl5M?usp=sharing). The notebook downloads Sentinel-1 Scenes from Alaska Satellite Facility (ASF) to crop the area and merge subswaths. Note: replace the scene names to produce an interferogram for your area of interest.
  
* [Yamchi DAM Interferograms Timeseries Analysis](https://colab.research.google.com/drive/16DYYnOcfC5PJk4kdEwZD8WECqfUZy_NS?usp=sharing)
  
  <img src="https://user-images.githubusercontent.com/7342379/182756524-c4c5aa0e-c6a7-4bda-b736-6ac872071433.jpg" width="50%">

See also a separate GitHub repository for the processing results: [YamchiDam](https://github.com/mobigroup/YamchiDam) That's combine two my software tools [PyGMTSAR](https://github.com/mobigroup/gmtsar) [N-Cube ParaView plugin for 3D/4D GIS Data Visualization](https://github.com/mobigroup/ParaView-plugins) for 4D analysis and visualization:

<img src="https://user-images.githubusercontent.com/7342379/144747743-a24d72ec-8875-4272-91f9-ec1f937bb798.gif" width="50%">

### Installation on MacOS

Use the same commands for BigSur on Intel chips and Monterey on Apple Silicon chips:

```
# create installation directory
sudo mkdir /usr/local/GMTSAR
sudo chown $(whoami) /usr/local/GMTSAR
# prepare system dependencies
brew install libtiff hdf5 gmt ghostscript
# install recent PyGMTSAR
cd /usr/local
git clone --branch master https://github.com/mobigroup/gmtsar GMTSAR
cd /usr/local/GMTSAR
autoconf
./configure --with-orbits-dir=/tmp
make
make install
```

Use your preferred way to install Python libraries (via PIP, Conda, HomeBrew packages and so on). Note: that's possible to install PyGMTSAR to /opt directory instead. For this case we need to change "GMTSAR" path in the notebooks. For compatibility reasons I save it the same for all the supported platforms.

### Installation on Debian and Ubuntu Linux

See the notebooks above where the installation commands included. Also, there is the cloud initialization scripit [GMTSAR.install.debian10.sh](https://github.com/mobigroup/gmtsar/blob/master/gmtsar/sh/GMTSAR.install.debian10.sh) to install and configure all the dependencies and GMTSAR on cloud and local Debian 10 hosts.

### About me

I have STEM master's degree in radio physics and in 2004 I was awarded first prize of the All-Russian Physics competition for significant results in Inverse modeling for non-linear optics and holography, also applicable for Inverse Modeling of Gravity, Magnetic, and Thermal fields. In addition to my fundamental science knowledge, I’m world class data scientist and software developer with 20 years experience in science and industrial development. I have worked on government contracts and universities projects and on projects for LG Corp, Google Inc, etc. You are able to find some of my software and results on LinkedIn and GitHub and Upwork, see the links below. By the way, I left Russia many years ago and I work remotely for about 20 years.

To order some research, development and support see my profile on freelance platform [Upwork](https://www.upwork.com/freelancers/~01e65e8e7221758623)

 [Geological models on YouTube channel](https://www.youtube.com/channel/UCSEeXKAn9f_bDiTjT6l87Lg)

 [Augmented Reality (AR) Geological Models](https://mobigroup.github.io/ParaView-Blender-AR/)

 [GitHub repositories](https://github.com/mobigroup)

 [English posts and articles on LinkedIn](https://www.linkedin.com/in/alexey-pechnikov/)

[Russian articles on Habr](https://habr.com/ru/users/N-Cube/posts/)

@ Alexey Pechnikov, 2022
