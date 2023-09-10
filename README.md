## Currently in development, PyGMTSAR2 (Python GMTSAR) is aimed at experts in Sentinel-1 Satellite Interferometry, offering features such as SBAS, PSI, PSI-SBAS, and more.

<img src="https://user-images.githubusercontent.com/7342379/194891967-be2b56b5-c30c-4040-8ef8-39b448ce2390.jpg" width="40%" />

## About PyGMTSAR2

PyGMTSAR allows accessible, reproducible, and powerful Sentinel-1 SBAS interferometry for everyone. But how about all the interferometry approaches like SBAS, PSI, PSI-SBAS, and time series and trend analysis, along with visualization, in a single Python package available on Google Colab, DockerHub, and beyond? This is the goal for PyGMTSAR2!

PyGMTSAR2 is under development. I share Jupyter notebook examples on [Patreon](https://www.patreon.com/pechnikov) and updates on its progress through my [LinkedIn](https://www.linkedin.com/in/alexey-pechnikov/).

### Preview Notebooks for Development PyGMTSAR2

One of the most awaited features in PyGMTSAR (Python InSAR) is the combined analysis of Persistent Scatterers (PS or PSI) and the Small Baseline Subset (SBAS).

Each of the PS and SBAS techniques has unique advantages and drawbacks — with SBAS performing better in rural areas and PS in urban ones. My vision involves merging the benefits of both methods and mitigating their shortcomings in a unified PS-SBAS process. In the development version, PyGMTSAR offers persistent scatterer analysis and weighted interferogram processing. This emphasizes stable pixels, enhancing phase and coherence. This not only improves the accuracy of results but also simplifies SBAS analysis by maintaining high coherence, even in rural areas.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/18KGGO9jy_EE8RejgYhKgqyxryhrQdmbS?usp=sharing) The updated SBAS notebook provides an illustration of the PS-SBAS analysis.

<img width="50%" src="https://raw.githubusercontent.com/mobigroup/articles/main/avg_corr_15m_adi.40%25.png">

## See Stable PyGMTSAR

The stable PyGMTSAR is available on GitHub, PyPI, DockerHub and Google Colab, see the project home page [PyGMTSAR GitHub Repository](https://github.com/mobigroup/gmtsar)

@ Alexey Pechnikov, 2023
