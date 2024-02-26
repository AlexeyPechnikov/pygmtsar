# https://jupyter-docker-stacks.readthedocs.io/en/latest/using/selecting.html
# https://github.com/jupyter/docker-stacks/blob/main/base-notebook/Dockerfile
# https://hub.docker.com/repository/docker/mobigroup/pygmtsar
# host platform compilation:
# docker build . -f pygmtsar.Dockerfile -t mobigroup/pygmtsar:latest --no-cache
# cross-compilation:
# docker buildx build . -f pygmtsar.Dockerfile -t mobigroup/pygmtsar:latest --no-cache --platform linux/amd64 --load
FROM jupyter/scipy-notebook:ubuntu-22.04

USER root

# install GMTSAR dependencies and some helpful command-line tools
RUN apt-get -y update \
&&  apt-get -y install git gdal-bin libgdal-dev subversion curl jq \
&&  apt-get -y install csh autoconf make gfortran \
&&  apt-get -y install libtiff5-dev libhdf5-dev liblapack-dev libgmt-dev gmt-dcw gmt-gshhg gmt \
&&  apt-get -y install zip htop mc netcdf-bin \
&&  apt-get clean && rm -rf /var/lib/apt/lists/*

# define installation and binaries search paths
ARG GMTSAR=/usr/local/GMTSAR
ARG ORBITS=/usr/local/orbits
ENV PATH=${GMTSAR}/bin:$PATH

# install GMTSAR from git
RUN cd $(dirname ${GMTSAR}) \
&&  git config --global advice.detachedHead false \
&&  git clone --branch master https://github.com/gmtsar/gmtsar GMTSAR \
&&  cd ${GMTSAR} \
&&  git checkout e98ebc0f4164939a4780b1534bac186924d7c998 \
&&  autoconf \
&&  ./configure --with-orbits-dir=${ORBITS} CFLAGS='-z muldefs' LDFLAGS='-z muldefs' \
&&  make \
&&  make install

# install PyGMTSAR and additional libraries
RUN apt-get -y update \
&&  apt-get -y install xvfb libegl1-mesa \
&&  apt-get clean && rm -rf /var/lib/apt/lists/*
RUN conda install -y -c conda-forge vtk panel xvfbwrapper pyvista
# magick fix for the library
RUN pip3 uninstall -y h5py
# use requirements.sh to build the installation command
RUN pip3 install \
    adjustText==1.0.4 \
    asf_search==7.0.4 \
    dask==2024.1.1 \
    distributed==2024.1.1 \
    geopandas==0.14.3 \
    h5netcdf==1.3.0 \
    h5py==3.10.0 \
    imageio==2.31.5 \
    ipywidgets==8.1.1 \
    joblib==1.3.2 \
    matplotlib==3.8.0 \
    nc-time-axis==1.4.1 \
    numba==0.57.1 \
    numpy==1.24.4 \
    pandas==2.2.1 \
    remotezip==0.12.2 \
    rioxarray==0.15.1 \
    scikit-learn==1.3.1 \
    scipy==1.11.4 \
    seaborn==0.13.0 \
    shapely==2.0.3 \
    statsmodels==0.14.0 \
    tqdm==4.66.1 \
    xarray==2024.2.0 \
    xmltodict==0.13.0 \
    pygmtsar

# modify start-notebook.py to start Xvfb
RUN sed -i '/import sys/a \
# Start Xvfb\n\
import xvfbwrapper\n\
display = xvfbwrapper.Xvfb(width=1280, height=1024)\n\
display.start()' /usr/local/bin/start-notebook.py

# grant passwordless sudo rights
RUN echo "${NB_USER} ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

# switch user
USER    ${NB_UID}
WORKDIR "${HOME}"

# Clone only the pygmtsar2 branch
RUN git clone --branch pygmtsar2 --single-branch https://github.com/AlexeyPechnikov/pygmtsar.git \
&& mv pygmtsar/notebooks ./notebooks \
&& mv pygmtsar/README.md ./ \
&& rm -rf pygmtsar work


