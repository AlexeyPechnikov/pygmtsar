## Announcements

The **e-book, titled 'PyGMTSAR: Sentinel-1 Python InSAR: An Introduction'** is now available for the stable PyGMTSAR release on various platforms, including [Amazon, Apple, Kobo, and many other bookstores](https://books2read.com/b/PyGMTSAR-introduction). If you'd like a preview of the content, you can check out the [PyGMTSAR Introduction Preview](https://github.com/mobigroup/gmtsar/blob/pygmtsar2/book/PyGMTSAR_preview.pdf) uploaded in the repository.

You have the option to **support the development of PyGMTSAR software on [Patreon](https://www.patreon.com/pechnikov) and [Buy Me a Coffee](https://www.buymeacoffee.com/pechnikov) platforms**. These platforms also offer additional documentation and use cases.

PyGMTSAR (Python InSAR) **Video Lessons and Educational Notebooks** available on [Patreon](https://www.patreon.com/collection/12458) and [YouTube](https://www.youtube.com/channel/UCSEeXKAn9f_bDiTjT6l87Lg).

## PyGMTSAR (Python InSAR) - Sentinel-1 Satellite Interferometry for Everyone

PyGMTSAR (Python InSAR) aims to cater to the needs of both occasional users and experts in Sentinel-1 Satellite Interferometry. It offers a range of features, including SBAS, PSI, PSI-SBAS, and more. It is available in [pygmtsar2 branch](https://github.com/mobigroup/gmtsar/tree/pygmtsar2). I  share Jupyter notebook examples on [Patreon](https://www.patreon.com/pechnikov) and updates on its progress through my [LinkedIn](https://www.linkedin.com/in/alexey-pechnikov/).

<img src="https://github.com/mobigroup/gmtsar/assets/7342379/3a7d8fda-a3e1-4282-b5ae-2b1c362b891d" width="40%" />

## About Development PyGMTSAR

PyGMTSAR offers accessible, reproducible, and powerful Sentinel-1 SBAS interferometry for everyone, regardless of their location. It encompasses various interferometry approaches, including SBAS, PSI, PSI-SBAS, and time series and trend analysis, all wrapped into a single Python package. Whether you're using Google Colab, DockerHub, or any other platform, PyGMTSAR is readily available for your needs.

The latest version of PyGMTSAR is currently in development, and I regularly share Jupyter notebook examples on [Patreon](https://www.patreon.com/pechnikov). You can also stay updated on its progress through my [LinkedIn profile](https://www.linkedin.com/in/alexey-pechnikov/). The new version is capable of performing powerful and fast analysis on hundreds of Sentinel-1 scenes and thousands of interferograms, although some functions and their options may be subject to renaming or changes.

One of the most awaited features in PyGMTSAR (Python InSAR) is the combined analysis of Persistent Scatterers (PS or PSI) and the Small Baseline Subset (SBAS). Each of the PS and SBAS techniques has unique advantages and drawbacks — with SBAS performing better in rural areas and PS in urban ones. My vision involves merging the benefits of both methods and mitigating their shortcomings in a unified PS-SBAS process. In the development version, PyGMTSAR offers persistent scatterer analysis and weighted interferogram processing. This emphasizes stable pixels, enhancing phase and coherence. This not only improves the accuracy of results but also simplifies SBAS analysis by maintaining high coherence, even in rural areas.

## PyGMTSAR Live Examples on Google Colab

Google Colab is a free service, and these notebooks offer interactive examples that are accessible directly in your web browser, available to everyone. You don't need a powerful computer, extensive disk space, a fast internet connection, or to install any required software. Almost any internet-connected device, such as a desktop, laptop, smartphone, or even a smart TV, is sufficient for InSAR processing using PyGMTSAR. Moreover, you can save the results and the processing Jupyter notebook on your local computer or server to run it locally or in the cloud.

All steps are automated, which includes software installation on Google Colab's cloud host (Linux Ubuntu 22, Python 3.10), downloading of Sentinel-1 SLC and orbit files from the Alaska Satellite Facility (ASF) datastore, obtaining SRTM DEM data and converting it to ellipsoidal heights using the EGM96 model, downloading a landmask for masking low-coherence water surfaces, and, of course, performing complete interferometry processing and result mapping. You can replace the scene names with your own to obtain similar results for your specific area. All the notebooks are accompanied by interactive 3D maps that are available instantly. 

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1TARVTB7z8goZyEVDRWyTAKJpyuqZxzW2?usp=sharing) CENTRAL Türkiye Mw 7.8 & 7.5 Earthquakes Co-Seismic Interferogram, 2023.

<img src="https://github.com/mobigroup/gmtsar/assets/7342379/cce39fa5-0115-467e-836d-8361a37da935" width="50%"><img src="https://github.com/mobigroup/gmtsar/assets/7342379/47543745-e7b1-41cb-b9f3-6f73cb1f9fb3" width="50%">

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1dDFG8BoF4WfB6tOF5sAi5mjdBKRbhxHo?usp=sharing) Pico do Fogo Volcano Eruption on Cape Verde's Fogo Island, 2014.

<img src="https://github.com/mobigroup/gmtsar/assets/7342379/d2eda089-0730-4699-82db-9410712d55ff" width="50%"><img src="https://github.com/mobigroup/gmtsar/assets/7342379/694d9670-36c9-4e56-bfb8-056e0d038d58" width="50%">

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1d9RcqBmWIKQDEwJYo8Dh6M4tMjJtvseC?usp=sharing) La Cumbre Volcano Eruption Interferogram, 2020.

<img src="https://github.com/mobigroup/gmtsar/assets/7342379/93cc9c5c-a654-4cc6-a310-2f3337c95ce2" width="50%"><img src="https://github.com/mobigroup/gmtsar/assets/7342379/fe085c2b-5bd5-4385-a1fe-04144568e1cb" width="50%">

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1shNGvUlUiXeyV7IcTmDbWaEM6XrB0014?usp=sharing) Iran–Iraq Earthquake Co-Seismic Interferogram, 2017.

<img src="https://github.com/mobigroup/gmtsar/assets/7342379/7357a56a-d69f-451b-91ab-367cbf2af410" width="50%"><img src="https://github.com/mobigroup/gmtsar/assets/7342379/bcd807f9-5d48-4bb4-ac13-803305f3b6da" width="50%">

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1h4XxJZwFfm7EC8NUzl34cCkOVUG2uJr4?usp=sharing) Imperial Valley SBAS analysis, 2015.

<img src="https://github.com/mobigroup/gmtsar/assets/7342379/bbe0f043-af09-4724-9e50-5549d3f24adc" width="50%"><img src="https://github.com/mobigroup/gmtsar/assets/7342379/cd1d8c33-3488-41af-aece-985b4d4202ae" width="50%">

## PyGMTSAR Live Examples on Google Colab Pro

Additionally, I share more complex SBAS and PSI use cases on my [Patreon](https://www.patreon.com/pechnikov) for subscribers. Please note that Google Colab Pro is a paid service, and accessing these examples requires a paid membership.

* InSAR analysis on Gastein Valley, Austria, 2021–2023. SBAS and PSI example featuring 58 Sentinel-1 SLC and between 200 to 1400 interferograms.

<img src="https://github.com/mobigroup/gmtsar/assets/7342379/27bd0436-db6f-45cd-88c8-99332ec337d5" width="50%"><img src="https://github.com/mobigroup/gmtsar/assets/7342379/6d1181fb-5bdb-4923-ae64-74f59b48f9b4" width="50%">

* InSAR analysis on Imperial Valley, California, USA, 2015. SBAS and PSI example featuring 58 Sentinel-1 SLC and between 200 to 1400 interferograms.
<img src="https://github.com/mobigroup/gmtsar/assets/7342379/c016c07e-06e9-4aef-ab59-4caa83d10541" width="100%">

## See Stable PyGMTSAR (previous version)

The stable PyGMTSAR is available on GitHub, PyPI, DockerHub and Google Colab, see the project home page [PyGMTSAR GitHub Repository](https://github.com/mobigroup/gmtsar)

@ Alexey Pechnikov, 2023
