[![View on GitHub](https://img.shields.io/badge/GitHub-View%20on%20GitHub-blue)](https://github.com/AlexeyPechnikov/pygmtsar)
[![Available on pypi](https://img.shields.io/pypi/v/pygmtsar.svg)](https://pypi.python.org/pypi/pygmtsar/)
[![Docker](https://badgen.net/badge/icon/docker?icon=docker&label)](https://hub.docker.com/r/pechnikov/pygmtsar)
[![DOI](https://zenodo.org/badge/398018212.svg)](https://zenodo.org/badge/latestdoi/398018212)
[![Support on Patreon](https://img.shields.io/badge/Patreon-Support-orange.svg)](https://www.patreon.com/pechnikov)
[![ChatGPT Assistant](https://img.shields.io/badge/ChatGPT-Assistant-green?logo=openai)](https://insar.dev/ai)

## PyGMTSAR (Python InSAR): Powerful and Accessible Satellite Interferometry

PyGMTSAR (Python InSAR) is designed to meet the needs of both occasional users and experts in Sentinel-1 Satellite Interferometry. It offers a wide range of features, including SBAS, PSI, PSI-SBAS, and more. In addition to the examples provided below, I also share Jupyter notebook examples on [Patreon](https://www.patreon.com/pechnikov) and provide updates on its progress through my [LinkedIn profile](https://www.linkedin.com/in/alexey-pechnikov/).

<img src="https://github.com/AlexeyPechnikov/pygmtsar/assets/7342379/c157c3a6-ed06-4b6d-82ae-c0aefb286d47" width="40%" />

## About PyGMTSAR

PyGMTSAR provides accessible, reproducible, and powerful Sentinel-1 interferometry that is available to everyone, regardless of their location. It encompasses a variety of interferometry approaches, including SBAS, PSI, PSI-SBAS, and time series and trend analysis, all integrated into a single Python package. Whether you're utilizing Google Colab, DockerHub, or any other platform, PyGMTSAR is ready to meet your needs.

One of the most in-demand features in PyGMTSAR (Python InSAR) is the combined analysis of Persistent Scatterers (PS or PSI) and the Small Baseline Subset (SBAS). Each method has its own unique advantages and drawbacks — SBAS typically performs better in rural areas, while PS is more suited to urban environments. My vision is to merge the benefits of both methods while mitigating their shortcomings through a unified PS-SBAS process. Additionally, PyGMTSAR offers weighted interferogram processing using an amplitude stability matrix, which emphasizes stable pixels. This approach enhances phase and coherence, improving the accuracy of results by maintaining high coherence, even in rural areas.

## PyGMTSAR Live Examples on Google Colab

Google Colab is a free service offering interactive notebooks that are accessible directly in your web browser and available to everyone. These notebooks provide live examples of InSAR processing using PyGMTSAR. You don't need a powerful computer, extensive disk space, a fast internet connection, or any special software installations. Almost any internet-connected device, including desktops, laptops, smartphones, or even smart TVs, can effectively handle InSAR processing with PyGMTSAR. Furthermore, you can save the results and the processing Jupyter notebook on your local computer or server to run it locally or in the cloud.

All steps in these notebooks are automated. This includes the software installation on Google Colab's cloud host (Linux Ubuntu 22, Python 3.10), downloading Sentinel-1 SLC and orbit files from the Alaska Satellite Facility (ASF) datastore, obtaining SRTM DEM data and converting it to ellipsoidal heights using the EGM96 model, downloading a land mask for masking low-coherence water surfaces, and of course, carrying out complete interferometry processing and result mapping. You can also customize the notebooks by replacing the scene names to process specific areas of your interest. Additionally, all notebooks are accompanied by interactive 3D maps that are available instantly.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1TARVTB7z8goZyEVDRWyTAKJpyuqZxzW2?usp=sharing) CENTRAL Türkiye Mw 7.8 & 7.5 Earthquakes Co-Seismic Interferogram, 2023.

<img src="https://github.com/AlexeyPechnikov/pygmtsar/assets/7342379/cce39fa5-0115-467e-836d-8361a37da935" width="50%"><img src="https://github.com/AlexeyPechnikov/pygmtsar/assets/7342379/47543745-e7b1-41cb-b9f3-6f73cb1f9fb3" width="50%">

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1dDFG8BoF4WfB6tOF5sAi5mjdBKRbhxHo?usp=sharing) Pico do Fogo Volcano Eruption on Cape Verde's Fogo Island, 2014.

<img src="https://github.com/AlexeyPechnikov/pygmtsar/assets/7342379/d2eda089-0730-4699-82db-9410712d55ff" width="50%"><img src="https://github.com/AlexeyPechnikov/pygmtsar/assets/7342379/694d9670-36c9-4e56-bfb8-056e0d038d58" width="50%">

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1d9RcqBmWIKQDEwJYo8Dh6M4tMjJtvseC?usp=sharing) La Cumbre Volcano Eruption Interferogram, 2020.

<img src="https://github.com/AlexeyPechnikov/pygmtsar/assets/7342379/93cc9c5c-a654-4cc6-a310-2f3337c95ce2" width="50%"><img src="https://github.com/AlexeyPechnikov/pygmtsar/assets/7342379/fe085c2b-5bd5-4385-a1fe-04144568e1cb" width="50%">

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1shNGvUlUiXeyV7IcTmDbWaEM6XrB0014?usp=sharing) Iran–Iraq Earthquake Co-Seismic Interferogram, 2017.

<img src="https://github.com/AlexeyPechnikov/pygmtsar/assets/7342379/7357a56a-d69f-451b-91ab-367cbf2af410" width="50%"><img src="https://github.com/AlexeyPechnikov/pygmtsar/assets/7342379/bcd807f9-5d48-4bb4-ac13-803305f3b6da" width="50%">

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1h4XxJZwFfm7EC8NUzl34cCkOVUG2uJr4?usp=sharing) Imperial Valley SBAS analysis, 2015.

<img src="https://github.com/AlexeyPechnikov/pygmtsar/assets/7342379/bbe0f043-af09-4724-9e50-5549d3f24adc" width="50%"><img src="https://github.com/AlexeyPechnikov/pygmtsar/assets/7342379/cd1d8c33-3488-41af-aece-985b4d4202ae" width="50%">

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1aqAr9KWKzGx9XpVie1M000C3vUxzNDxu?usp=sharing) Flooding [Correlation] Map: Kalkarindji, NT Australia, 2024.

<img src="https://github.com/AlexeyPechnikov/pygmtsar/assets/7342379/368e5fc2-1966-4f98-a03d-e82b50103c05" width="100%">

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1ipiQGbvUF8duzjZER8v-_R48DSpSmgvQ?usp=sharing) PyGMTSAR SBAS and PSI Analyses: Golden Valley, CA.
<img src="https://github.com/AlexeyPechnikov/pygmtsar/assets/7342379/8b416787-4b81-44f8-8956-3a5d596af51b" width="100%">

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1O3aZtZsTrQIldvCqlVRel13wJRLhmTJt?usp=sharing) PyGMTSAR SBAS and PSI Analyses: Lake Sarez Landslides, Tajikistan.
<img src="https://github.com/AlexeyPechnikov/pygmtsar/assets/7342379/b507cad0-db7a-47e6-a679-f74631c5e840" width="100%">

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/19PLuebOZ4gaYX5ym1H7SwUbJKfl23qPr?usp=sharing) PyGMTSAR Elevation Map:  Erzincan, Türkiye.
<img src="https://github.com/AlexeyPechnikov/pygmtsar/assets/7342379/066532d5-7b07-49d2-9478-7b8f966a3752" width="100%">

## PyGMTSAR Live Examples on Google Colab Pro

For subscribers, I share more complex SBAS and PSI use cases on Google Colab Pro through my [Patreon page](https://www.patreon.com/pechnikov). These use cases are suitable for InSAR learners, researchers, and industry specialists working on their challenging projects. Large areas and big stacks for thousands of interferograms, low-coherence territories, and extensive atmospheric phase delays - all these tasks can be addressed with PyGMTSAR. These examples can still be run online on the Google Colab Pro platform, which is cost-effective ($10/month) and provides a good balance between very fast data transfer speeds for downloading dozens of Sentinel-1 SLC scenes, available disk space to store the datasets and process them (approximately 220GB vs. 110GB for the free version of Google Colab), processing speed (8 vCPUs vs. 2 for the free version of Google Colab), and accessible memory (54GB vs. 12GB for the free version of Google Colab). I frequently utilize Google Colab Pro myself to manage up to five parallel InSAR projects, without concerns about disk space, memory, or processing performance limitations. Moreover, all the examples can be executed locally as well as on cloud hosts and remote servers.

## Projects and Publications Using PyGMTSAR

Explore the diverse applications of PyGMTSAR in projects and academic research on the dedicated [Projects and Publications](/pubs/README.md) page.

## Announcements

**E-Book Release: 'PyGMTSAR: Sentinel-1 Python InSAR: An Introduction'**
The e-book is now available for the stable PyGMTSAR release across various platforms, including [Amazon, Apple, Kobo, and many other bookstores](https://books2read.com/b/PyGMTSAR-introduction). For a glimpse of the content, check out the [PyGMTSAR Introduction Preview](https://github.com/AlexeyPechnikov/pygmtsar/blob/pygmtsar2/book/PyGMTSAR_preview.pdf) in the GitHub repository.

<img src="https://github.com/AlexeyPechnikov/pygmtsar/assets/7342379/93859fc8-f867-41d0-a03a-0fd89d854e82" width="40%">

**Educational Resources: Video Lessons and Notebooks**
Find PyGMTSAR (Python InSAR) video lessons and educational notebooks on [Patreon](https://www.patreon.com/collection/12458) and [YouTube](https://www.youtube.com/channel/UCSEeXKAn9f_bDiTjT6l87Lg).

**PyGMTSAR AI Assistant**
The [PyGMTSAR AI Assistant](https://insar.dev/ai), powered by OpenAI GPT-4, is knowledgeable in InSAR processing using PyGMTSAR. It can assist in understanding the theory, finding and explaining InSAR examples, creating an InSAR processing pipeline, and troubleshooting issues in your processing.

<img width="40%" alt="PyGMTSAR AI Assistant" src="https://github.com/AlexeyPechnikov/pygmtsar/assets/7342379/f1b7780d-9a93-4277-b9c3-9e54d9ff3d8b">

The assistant can answer many of your questions, such as:

* How to start with InSAR?

* Where can I find interactive InSAR example?

* Please provide interferogram creation code.

* Show me online InSAR examples on Google Colab.

* Explain to me content https:// [colab.research.google.com/drive/1673p-BhRwsh8g3VBYhqBYLrL5Lso81mj?usp=sharing](http://colab.research.google.com/drive/1673p-BhRwsh8g3VBYhqBYLrL5Lso81mj?usp=sharing)

* Show me open tickets.

* Find the recent ticket about Docker images and display last message.

* Create my AOI as GeoJSON text for a line between the points (-24.42, 14.8) and (-24.54, 14.88).

* Could you explain the global plotting parameters used in https://colab.research.google.com/drive/1dpDWbp3BO-xVWnTcJN4NXTdfZ47oxrM4?usp=sharing

* What specific lines of code need to be modified to compute the interferogram without multilooking in https://colab.research.google.com/drive/1dpDWbp3BO-xVWnTcJN4NXTdfZ47oxrM4?usp=sharing

Furthermore, you have the option to upload a document or a screenshot for discussion, and you can request explanations, such as 'explain the code to me,' among many other possibilities.

## PyGMTSAR Previous Version

The 2023 releases of PyGMTSAR are still available on GitHub, PyPI, DockerHub, and Google Colab. For more information and access to these releases, visit the project's home page at the [PyGMTSAR 2023 GitHub Repository](https://github.com/AlexeyPechnikov/pygmtsar/tree/pygmtsar). Included is a collection of examples that facilitate the comparison of PyGMTSAR's InSAR processing capabilities with those of other InSAR software.

@ Alexey Pechnikov, 2024
