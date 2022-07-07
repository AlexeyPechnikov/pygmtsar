## PyGMTSAR (Python GMTSAR) - Easy and Fast Satellite Interferometry For Everyone

The goal of the project is easy and fast satellite interferometry (InSAR) processing everywhere
from local host to cloud environments. GMTSAR binary command line tools are used under the hood. By my opinion, GMTSAR is great project and has only one big problem - it is not user friendly. I've built rich Python API instead of GMTSAR CSH scripts and GMT toolkit calls and by this way we have 3x faster processing and more accurate results. Really, now we are able to use even free of charge services like to Google Colab which are powerful enough for processing about 10 interferograms in 30 minutes.
MacOS and Linux Debian are my prefered OSes and I support both of them. Please don't ask me about Windows support while you don't ready to pay for it.

PyGMTSAR automatically downloads Sentinel-1 orbit files and SRTM DEM (and converts it to ellispoidal heights using EGM96 model) and so on. You need just 2+ raw Sentinel scenes for the processing. See below **Live Example S1A_Stack_CPGF_T173 on Google Colab comparision to report from Centre of EO Research & Satellite Remote Sensing, Greece** for the fully automated processing and other examples for faster processing using pre-downloaded DEM and orbits.

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

Note: that's possible to install PyGMTSAR to /opt directory instead. For this case we need to change "GMTSAR" path in the notebooks.

### Installation on Debian and Ubuntu Linux

See the notebooks below where the installation commands included.

### Live Examples

* [Live Example Kumamoto Earthquake on Apr, 2016 Co-Seismic Interferogram on Google Colab](https://colab.research.google.com/drive/1yCeYYfh-HPehpld1lWJ860MBCkq0JDE4?usp=sharing)

<img src="https://user-images.githubusercontent.com/7342379/177028662-c24463dd-f58d-4be0-80e7-f9d0b38e8eb4.jpg" width="50%">

* [Live Example S1A_Crete_Earthquake on Google Colab](https://colab.research.google.com/drive/1reRd-BJxa3Vxz_hmCMrn_Jv1Dbpazwa-?usp=sharing)

<img src="https://user-images.githubusercontent.com/7342379/177004243-123f7127-2c42-48d4-996b-b4200929a011.jpg" width="50%">

* [Live Example S1A_Stack_CPGF_T173 on Google Colab comparision to report from Centre of EO Research & Satellite Remote Sensing, Greece](https://colab.research.google.com/drive/1asddx-b3f7jS6TK-fkTm_WuytqBZzovc?usp=sharing)

<img src="https://user-images.githubusercontent.com/7342379/177004287-cdd4351c-0834-42ae-8e46-9da5e8b124bf.jpg" width="50%">

* [Live Example S1A_Stack_CPGF_T173 on Google Colab](https://colab.research.google.com/drive/1Ab2I9A2kmZIHuNwIortaga1zZ90hBjMa?usp=sharing)

<img src="https://user-images.githubusercontent.com/7342379/135814732-aa0eb142-ae54-4a57-b271-c33b5174a28e.png" width="50%">

<img src="https://user-images.githubusercontent.com/7342379/136558388-cffbcea5-e5a7-44d5-ba36-b51a7e0f10e9.png" width="50%">

* [Yamchi DAM Interferograms Persistent Scatterer Interferometry (PSI) Analysis](https://colab.research.google.com/drive/1ant72nEGxARIqxkXfVvwoMg1yxEkImrr?usp=sharing)
See also a separate GitHub repository for the processing results: [YamchiDam](https://github.com/mobigroup/YamchiDam)

<img src="https://user-images.githubusercontent.com/7342379/144747743-a24d72ec-8875-4272-91f9-ec1f937bb798.gif" width="50%">

<img src="https://user-images.githubusercontent.com/7342379/143826953-987c115e-21d7-4396-a1b2-3ecff85dad7e.jpg" width="50%">

> I'm a freelancer and that's my free time Open Source project. In 2005 my master's thesis was awarded first prize of the
> All-Russian Physics competition for significant results in Inverse modeling for non-linear optics and holography and so
> I know a lot about interferometry modeling. Sure, I have varios ideas about new features but my free time is very limited
> and I can't promise anything. You are able to sponsor my projects on [Patreon: Become a Patron!](https://www.patreon.com/bePatron?u=54500608) and order research, development and support on [Upwork](https://www.upwork.com/freelancers/~01e65e8e7221758623)
>
> @ Alexey Pechnikov, 2021
>
> [Geological models on YouTube channel](https://www.youtube.com/channel/UCSEeXKAn9f_bDiTjT6l87Lg)
>
> [Augmented Reality (AR) Geological Models](https://mobigroup.github.io/ParaView-Blender-AR/)
>
> [GitHub repositories](https://github.com/mobigroup)
>
> [English posts and articles on LinkedIn](https://www.linkedin.com/in/alexey-pechnikov/)
>
> [Russian articles on Habr](https://habr.com/ru/users/N-Cube/posts/)
