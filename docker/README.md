[![View on GitHub](https://img.shields.io/badge/GitHub-View%20on%20GitHub-blue)](https://github.com/AlexeyPechnikov/pygmtsar)
[![Available on pypi](https://img.shields.io/pypi/v/pygmtsar.svg)](https://pypi.python.org/pypi/pygmtsar/)
[![DOI](https://zenodo.org/badge/398018212.svg)](https://zenodo.org/badge/latestdoi/398018212)
[![Support on Patreon](https://img.shields.io/badge/Patreon-Support-orange.svg)](https://www.patreon.com/pechnikov)
[![ChatGPT Assistant](https://img.shields.io/badge/ChatGPT-Assistant-green?logo=openai)](https://insar.dev/ai)

## PyGMTSAR (Python InSAR) - Powerful and Accessible Satellite Interferometry

**Note**: Previous builds available at https://hub.docker.com/r/mobigroup/pygmtsar

PyGMTSAR (Python InSAR) is designed to meet the needs of both occasional users and experts in Sentinel-1 Satellite Interferometry. It offers a wide range of features, including SBAS, PSI, PSI-SBAS, and more. In addition to the examples provided below, I also share Jupyter notebook examples on [Patreon](https://www.patreon.com/pechnikov) and provide updates on its progress through my [LinkedIn profile](https://www.linkedin.com/in/alexey-pechnikov/).

<img src="https://github.com/AlexeyPechnikov/pygmtsar/assets/7342379/c157c3a6-ed06-4b6d-82ae-c0aefb286d47" width="40%" />

## About PyGMTSAR

PyGMTSAR provides accessible, reproducible, and powerful Sentinel-1 interferometry that is available to everyone, regardless of their location. It encompasses a variety of interferometry approaches, including SBAS, PSI, PSI-SBAS, and time series and trend analysis, all integrated into a single Python package. Whether you're utilizing Google Colab, DockerHub, or any other platform, PyGMTSAR is ready to meet your needs.

One of the most in-demand features in PyGMTSAR (Python InSAR) is the combined analysis of Persistent Scatterers (PS or PSI) and the Small Baseline Subset (SBAS). Each method has its own unique advantages and drawbacks — SBAS typically performs better in rural areas, while PS is more suited to urban environments. My vision is to merge the benefits of both methods while mitigating their shortcomings through a unified PS-SBAS process. Additionally, PyGMTSAR offers weighted interferogram processing using an amplitude stability matrix, which emphasizes stable pixels. This approach enhances phase and coherence, improving the accuracy of results by maintaining high coherence, even in rural areas.

## PyGMTSAR Live Examples in Docker Image

Configure your Docker runtime (Preferences -> Resources tab for Docker Desktop) to use 4 CPU cores and 8 GB of RAM, plus 1 GB of swap space and a 300 GB disk image to execute all the examples. Some of the examples also work with just 2 CPU cores, 4 GB RAM, and a disk image size of 30 GB.

It's important to note that the examples are coded with different objectives in mind. While some notebooks are designed to provide easy-to-read educational content, others are optimized for the most efficient execution possible. This diversity in coding approaches explains why an annual SBAS+PSI example can be run on a system with 4GB RAM, while a SBAS analysis of a few interferograms may require 8GB RAM.

**Table: Processing Times for InSAR Analyses Example Notebooks on an iMac 2021 (Apple M1, 16 GB RAM, 2 TB SSD) Using Various Docker Configurations** 

| Analysis                 | Notebook                      | Scenes | Subswaths | Interferograms | Time (Native, no Docker) | Time (4 CPUs, 8 GB RAM) | Time (2 CPUs, 4 GB RAM) |
| ------------------------ | ----------------------------- | ------ | --------- | -------------- | ------------------------ | ----------------------- | ----------------------- |
| SBAS and PSI Analyses    | Lake Sarez Landslides         | 19     | 1         | 76             | 31 min                   | 50 min                  | -                       |
| Co-Seismic Interferogram | CENTRAL Türkiye Earthquake    | 4      | 3         | 1              | 12 min                   | 28 min                  | -                       |
| SBAS and PSI Analyses    | Golden Valley Subsidence      | 30     | 1         | 57             | 11 min                   | 19 min                  | 23 min                  |
| SBAS Analysis            | Imperial Valley Groundwater   | 5      | 1         | 9              | 4 min                    | 8 min                   | -                       |
| Co-Seismic Interferogram | Iran–Iraq Earthquake          | 3      | 1         | 1              | 1 min                    | 6 min                   | 6 min                   |
| Co-Seismic Interferogram | La Cumbre Volcano Eruption    | 2      | 2         | 1              | < 1 min                  | 1 min                   | 1 min                   |
| Co-Seismic Interferogram | Pico do Fogo Volcano Eruption | 2      | 1         | 1              | < 1 min                  | 1 min                   | 1 min                   |
| Flooding Map             | Kalkarindji                   | 3      | 1         | 2              | < 1 min                  | 1 min                   | 1 min                   |

Download the Docker image (or build it yourself using the Dockerfile in the repository), and run the container while forwarding port 8888 to JupyterLab using these commands inside your command line terminal window:

```
docker pull pechnikov/pygmtsar
docker run -dp 8888:8888 --name pygmtsar docker.io/pechnikov/pygmtsar
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

The commands below build multi-arch images using the [Dockerfile](https://github.com/AlexeyPechnikov/pygmtsar/blob/pygmtsar2/docker/pygmtsar.Dockerfile) from the project's GitHub repository and share them to DockerHub in the 'pechnikov' repository with tags 'latest' and the build date (like '2024-01-21'):

```
docker buildx create --name pygmtsar
docker buildx use pygmtsar
docker buildx inspect --bootstrap
docker buildx build . -f pygmtsar.Dockerfile \
    --platform linux/amd64,linux/arm64 \
    --tag pechnikov/pygmtsar:$(date "+%Y-%m-%d") \
    --tag pechnikov/pygmtsar:latest \
    --pull --push --no-cache
docker buildx rm pygmtsar
```

And this simpler command builds a Docker image for local use:

```
docker build . -f pygmtsar.Dockerfile -t pygmtsar:latest --no-cache
```

The only requirement for the builds is the [Dockerfile](https://github.com/mobigroup/gmtsar/blob/pygmtsar2/docker/pygmtsar.Dockerfile); there is no need to download the complete PyGMTSAR GitHub repository.

## Learn more

- Source code and bug tracker on [Github](https://github.com/AlexeyPechnikov/pygmtsar)
- Python library on [PyPI](https://pypi.python.org/pypi/pygmtsar/)
- Interactive examples online on Google Colab at [InSAR.dev](https://insar.dev)
- Superior online examples on Google Colab Pro and articles at [pechnikov.dev](https://pechnikov.dev)
- AI InSAR Assistant powered by ChatGPT4 at [InSAR.dev/ai](https://insar.dev/ai)
