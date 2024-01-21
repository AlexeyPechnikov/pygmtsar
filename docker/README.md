[![MacOS tests](https://github.com/mobigroup/gmtsar/actions/workflows/macos.yml/badge.svg)](https://github.com/mobigroup/gmtsar/actions/workflows/macos.yml)
[![Ubuntu tests](https://github.com/mobigroup/gmtsar/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/mobigroup/gmtsar/actions/workflows/ubuntu.yml)
[![Available on pypi](https://img.shields.io/pypi/v/pygmtsar.svg)](https://pypi.python.org/pypi/pygmtsar/)
[![DOI](https://zenodo.org/badge/398018212.svg)](https://zenodo.org/badge/latestdoi/398018212)

## PyGMTSAR (Python InSAR) - Sentinel-1 Satellite Interferometry for Everyone

PyGMTSAR (Python InSAR) is designed to meet the needs of both occasional users and experts in Sentinel-1 Satellite Interferometry. It offers a wide range of features, including SBAS, PSI, PSI-SBAS, and more. In addition to the examples provided below, I also share Jupyter notebook examples on [Patreon](https://www.patreon.com/pechnikov) and provide updates on its progress through my [LinkedIn profile](https://www.linkedin.com/in/alexey-pechnikov/).

<img src="https://github.com/mobigroup/gmtsar/assets/7342379/3a7d8fda-a3e1-4282-b5ae-2b1c362b891d" width="40%" />

## About PyGMTSAR

PyGMTSAR provides accessible, reproducible, and powerful Sentinel-1 interferometry that is available to everyone, regardless of their location. It encompasses a variety of interferometry approaches, including SBAS, PSI, PSI-SBAS, and time series and trend analysis, all integrated into a single Python package. Whether you're utilizing Google Colab, DockerHub, or any other platform, PyGMTSAR is ready to meet your needs.

One of the most in-demand features in PyGMTSAR (Python InSAR) is the combined analysis of Persistent Scatterers (PS or PSI) and the Small Baseline Subset (SBAS). Each method has its own unique advantages and drawbacks — SBAS typically performs better in rural areas, while PS is more suited to urban environments. My vision is to merge the benefits of both methods while mitigating their shortcomings through a unified PS-SBAS process. Additionally, PyGMTSAR offers weighted interferogram processing using an amplitude stability matrix, which emphasizes stable pixels. This approach enhances phase and coherence, improving the accuracy of results by maintaining high coherence, even in rural areas.

## PyGMTSAR Live Examples in Docker Image

Configure your Docker runtime (Preferences -> Resources tab for Docker Desktop) to use 2 CPU cores and 4 GB of RAM (or 4 CPU cores and 8 GB of RAM, and so on), plus 1 GB of swap and at least 16 GB of disk image size. Download the Docker image (or build it yourself using the Dockerfile in the repository), and run the container while forwarding port 8888 to JupyterLab using these commands inside your command line terminal window:

```
docker pull mobigroup/pygmtsar
docker run -dp 8888:8888 --name pygmtsar docker.io/mobigroup/pygmtsar
docker logs pygmtsar
```

See the output for the JupyterLab link and copy and past it into your web browser address line. Also, the donwloaded Docker image can be started in Docker Desktop app - press "RUN" button and define the container name and the port in the opened dialog window (see "Optional settings" for the port number input field) and click on the newly created container to launch it and see the output log with the clickable link.

System operations are available using password-less “sudo” command. To upgrade PyGMTSAR Python library, execute the command below in a notebook cell and restart the notebook:

```
import sys
!sudo {sys.executable} -m pip install -U pygmtsar
```

Alternatively, use the following terminal command in the Docker container's Terminal:

```
sudo --preserve-env=PATH sh -c "pip3 install -U pygmtsar"
```

## Build Docker images

The commands below build multi-arch images using the [Dockerfile](https://github.com/mobigroup/gmtsar/blob/pygmtsar2/docker/pygmtsar.Dockerfile) from the project's GitHub repository and share them to DockerHub in the 'mobigroup' repository with tags 'latest' and the build date (like '2024-01-21'):

```
docker buildx create --name mobigroup
docker buildx use mobigroup
docker buildx inspect --bootstrap
docker buildx build . -f pygmtsar.Dockerfile \
    --platform linux/amd64,linux/arm64 \
    --tag mobigroup/pygmtsar:$(date "+%Y-%m-%d") \
    --tag mobigroup/pygmtsar:latest \
    --pull --push --no-cache
```

And this simpler command builds a Docker image for local use:

```
docker build . -f pygmtsar.Dockerfile -t mobigroup/pygmtsar:latest --no-cache
```

The only requirement for the builds is the [Dockerfile](https://github.com/mobigroup/gmtsar/blob/pygmtsar2/docker/pygmtsar.Dockerfile); there is no need to download the complete PyGMTSAR GitHub repository.

## Learn more

- Google Colab Examples on GitHub: [https://github.com/mobigroup/gmtsar](https://github.com/mobigroup/gmtsar)

- PyPI Python library: [https://pypi.org/project/pygmtsar/](https://pypi.org/project/pygmtsar/)

- Google Colab Pro Examples on Patreon: [https://www.patreon.com/pechnikov](https://www.patreon.com/pechnikov)

- Issue tracker: [https://github.com/mobigroup/gmtsar/issues](https://github.com/mobigroup/gmtsar/issues)
