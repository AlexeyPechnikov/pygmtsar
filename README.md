## Announcements

The **e-book, titled 'PyGMTSAR: Sentinel-1 Python InSAR: An Introduction'** is now available for the stable PyGMTSAR release on various platforms, including [Amazon, Apple, Kobo, and many other bookstores](https://books2read.com/b/PyGMTSAR-introduction). If you'd like a preview of the content, you can check out the [PyGMTSAR Introduction Preview](https://github.com/mobigroup/gmtsar/blob/pygmtsar2/book/PyGMTSAR_preview.pdf) uploaded in the repository.

You have the option to **support the development of PyGMTSAR software on [Patreon](https://www.patreon.com/pechnikov) and [Buy Me a Coffee](https://www.buymeacoffee.com/pechnikov) platforms**. These platforms also offer additional documentation and use cases.

PyGMTSAR (Python InSAR) **Video Lessons and Educational Notebooks** available on [Patreon](https://www.patreon.com/collection/12458) and [YouTube](https://www.youtube.com/channel/UCSEeXKAn9f_bDiTjT6l87Lg).

## PyGMTSAR (Python InSAR) - Sentinel-1 Satellite Interferometry for Everyone

In its current development phase, PyGMTSAR (Python InSAR) aims to cater to the needs of both occasional users and experts in Sentinel-1 Satellite Interferometry. It offers a range of features, including SBAS, PSI, PSI-SBAS, and more. It is available in [pygmtsar2 branch](https://github.com/mobigroup/gmtsar/tree/pygmtsar2). I  share Jupyter notebook examples on [Patreon](https://www.patreon.com/pechnikov) and updates on its progress through my [LinkedIn](https://www.linkedin.com/in/alexey-pechnikov/).

<img src="https://user-images.githubusercontent.com/7342379/194891967-be2b56b5-c30c-4040-8ef8-39b448ce2390.jpg" width="40%" />

## About Development PyGMTSAR

PyGMTSAR offers accessible, reproducible, and powerful Sentinel-1 SBAS interferometry for everyone, regardless of their location. It encompasses various interferometry approaches, including SBAS, PSI, PSI-SBAS, and time series and trend analysis, all wrapped into a single Python package. Whether you're using Google Colab, DockerHub, or any other platform, PyGMTSAR is readily available for your needs.

The latest version of PyGMTSAR is currently in development, and I regularly share Jupyter notebook examples on [Patreon](https://www.patreon.com/pechnikov). You can also stay updated on its progress through my [LinkedIn profile](https://www.linkedin.com/in/alexey-pechnikov/). The new version is capable of performing powerful and fast analysis on hundreds of Sentinel-1 scenes and thousands of interferograms, although some functions and their options may be subject to renaming or changes.

One of the most awaited features in PyGMTSAR (Python InSAR) is the combined analysis of Persistent Scatterers (PS or PSI) and the Small Baseline Subset (SBAS). Each of the PS and SBAS techniques has unique advantages and drawbacks — with SBAS performing better in rural areas and PS in urban ones. My vision involves merging the benefits of both methods and mitigating their shortcomings in a unified PS-SBAS process. In the development version, PyGMTSAR offers persistent scatterer analysis and weighted interferogram processing. This emphasizes stable pixels, enhancing phase and coherence. This not only improves the accuracy of results but also simplifies SBAS analysis by maintaining high coherence, even in rural areas.

## PyGMTSAR Live Examples on Google Colab

These notebooks provide interactive examples directly in your web browser. All steps are automated, including software installation on Google Colab's cloud host, downloading of Sentinel-1 orbit files, SRTM DEM (and its conversion to ellipsoidal heights using the EGM96 model), a landmask (to mask low-coherence water surfaces), Sentinel-1 SLC scenes from the Alaska Satellite Facility (ASF) datastore, and of course, the complete interferometry processing and result mapping.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1shNGvUlUiXeyV7IcTmDbWaEM6XrB0014?usp=sharing) **ASF Downloading 2017 Iran–Iraq Earthquake Co-Seismic Interferogram**. Note: you can replace the scene names to generate an **interferogram** and **LOS displacement** for your area of interest. To compare the results with GMTSAR, SNAP, and GAMMA Software see the same notebook for the stable PyGMTSAR.

<img src="https://user-images.githubusercontent.com/7342379/177748605-788889e5-9afd-44d8-bc3c-dc6efe920ea0.png" width="50%">

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1h4XxJZwFfm7EC8NUzl34cCkOVUG2uJr4?usp=sharing) **GMTSAR example dataset S1A_Stack_CPGF_T173** This example demonstrates the SBAS analysis and detrending approach used to remove atmospheric noise, resulting in significantly improved outcomes.

<img src="https://user-images.githubusercontent.com/7342379/135814732-aa0eb142-ae54-4a57-b271-c33b5174a28e.png" width="50%">

<img src="https://user-images.githubusercontent.com/7342379/189961167-bf3901e5-417c-41ce-a5ca-d1c74c239a04.png" width="50%">

## See Stable PyGMTSAR

The stable PyGMTSAR is available on GitHub, PyPI, DockerHub and Google Colab, see the project home page [PyGMTSAR GitHub Repository](https://github.com/mobigroup/gmtsar)

@ Alexey Pechnikov, 2023
