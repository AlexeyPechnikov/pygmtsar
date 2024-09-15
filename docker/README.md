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

Configure your Docker runtime (Preferences -> Resources tab for Docker Desktop) using the example configurations in the table below, or adjust it to fit your needs. The disk usage for the examples is shown in the table. To process all the examples, you will need a 120 GB Docker virtual disk limit (114 GB Docker container size in my tests, though I recommend allowing slightly more to be safe). You can check the Docker container’s disk usage with the following command:

```
docker ps -s
```

The examples are coded with different objectives in mind. Some notebooks are designed to provide easy-to-read educational content, while others are optimized for maximum execution efficiency. This variation in coding approaches is why an annual SBAS+PSI example can run on a system with 2GB RAM, whereas an SBAS analysis of a few interferograms may require 8GB RAM.

Some examples use complete Sentinel-1 SLC scenes, while others utilize Sentinel-1 bursts. PyGMTSAR can effectively process both, allowing you to choose the type of source data that best suits your needs.

**Table: Processing Times for InSAR Analysis Example Notebooks on iMac 2021 (Apple M1, 8 CPU cores, 16 GB RAM, 2 TB SSD) Using Various Docker Configurations and Natively Without Docker**

| Analysis                 | Notebook                      | **Scn** | **Brst** | **Subswth** | **Intf** | Time (Native, no Docker) | Time (4 CPUs, 8 GB RAM) | Time (2 CPUs, 4 GB RAM) | Time (1 CPUs, 2 GB RAM) | Disk Usage, GB |
| ------------------------ | ----------------------------- | -------- | --------- | ----------- | ---------------------- | ------------------------ | ----------------------- | ----------------------- | ----------------------- | -------------- |
| SBAS and PSI Analyses    | Lake Sarez Landslides         | 19       | 38        | 1           | 76                     | 11 min                   | 17 min                  | 24 min                  | 38 min*      | 22.2           |
| Co-Seismic Interferogram | CENTRAL Türkiye Earthquake    | 4        | 112       | 3           | 1                      | 15 min                   | 24 min                  | 33 min                  | 62 min*      | 50.3           |
| SBAS and PSI Analyses    | Golden Valley Subsidence      | 30       | 30        | 1           | 57                     | 5 min                    | 7 min                   | 10 min                  | 18 min                  | 12.8           |
| SBAS Analysis            | Imperial Valley Groundwater   | 5        |           | 1           | 9                      | 6 min                    | 9 min                   | 12 min                  | 21 min                  | 5.5            |
| Co-Seismic Interferogram | Iran–Iraq Earthquake          | 3        |           | 1           | 1                      | 1 min                    | 2 min                   | 2 min                   | 3 min*       | 6.4            |
| Co-Seismic Interferogram | La Cumbre Volcano Eruption    | 2        | 8         | 2           | 1                      | 1 min                    | 1 min                   | 1 min                   | 2 min                   | 2.4            |
| Co-Seismic Interferogram | Pico do Fogo Volcano Eruption | 2        |           | 1           | 1                      | 1 min                    | 1 min                   | 1 min                   | 2 min*       | 2.1            |
| Flooding Map             | Kalkarindji                   | 3        |           | 1           | 2                      | 1 min                    | 1 min                   | 1 min                   | 2 min*       | 5.3            |
| Elevation Map            | Erzincan, Türkiye             | 2        |           | 1           | 1                      | 5 min                    | 7 min                   | 9 min                   | 16 min*      | 5.3            |

**Note**: Download time is excluded. The reported times reflect the complete notebook run time when all data has already been downloaded in a previous run. For stable downloading of scenes and bursts with 2GB or 4GB RAM configurations, set the parameter n_jobs=1 for the ASF.download() function call. The default Docker Desktop 1GB swap is used, except for some examples on a 2GB RAM configuration that require a 2GB swap, indicated by an asterisk (*).

Download the Docker image (or build it yourself using the [Dockerfile](https://github.com/AlexeyPechnikov/pygmtsar/blob/pygmtsar2/docker/pygmtsar.Dockerfile) in the repository), and run the container while forwarding port 8888 to JupyterLab using these commands inside your command line terminal window:

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

The commands below build multi-architecture images using the [Dockerfile](https://github.com/AlexeyPechnikov/pygmtsar/blob/pygmtsar2/docker/pygmtsar.Dockerfile) from the project’s GitHub repository and push them to DockerHub in the ‘pechnikov’ repository with the tags ‘latest’ and the build date (e.g., ‘2024-01-21’):

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

For a local build, use the following command:

```
docker build . -f pygmtsar.Dockerfile -t pygmtsar:latest --no-cache
```

The only requirement for these builds is the [Dockerfile](https://github.com/AlexeyPechnikov/pygmtsar/blob/pygmtsar2/docker/pygmtsar.Dockerfile); there is no need to download the entire PyGMTSAR GitHub repository, as the build process will automatically retrieve all necessary files.

## Learn more

- Source code and bug tracker on [Github](https://github.com/AlexeyPechnikov/pygmtsar)
- Python library on [PyPI](https://pypi.python.org/pypi/pygmtsar/)
- Interactive examples online on Google Colab at [InSAR.dev](https://insar.dev)
- Superior online examples on Google Colab Pro and articles at [pechnikov.dev](https://pechnikov.dev)
- AI InSAR Assistant powered by ChatGPT4 at [InSAR.dev/ai](https://insar.dev/ai)
