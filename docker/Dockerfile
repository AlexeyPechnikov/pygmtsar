# https://jupyter-docker-stacks.readthedocs.io/en/latest/using/selecting.html
# https://github.com/jupyter/docker-stacks/blob/main/base-notebook/Dockerfile
FROM jupyter/scipy-notebook:ubuntu-22.04

USER root

# install GMTSAR dependencies
RUN apt-get -y update \
&&  apt-get -y install git gdal-bin libgdal-dev subversion curl jq \
&&  apt-get -y install csh autoconf make gfortran \
&&  apt-get -y install libtiff5-dev libhdf5-dev liblapack-dev libgmt-dev gmt-dcw gmt-gshhg gmt \
&&  apt-get clean && rm -rf /var/lib/apt/lists/*

# define installation paths
ARG GMTSAR=/usr/local/GMTSAR
ARG ORBITS=/usr/local/orbits

# install GMTSAR from git
RUN cd /usr/local && git clone --branch master https://github.com/gmtsar/gmtsar GMTSAR
RUN cd ${GMTSAR} \
&&  autoconf \
&&  ./configure --with-orbits-dir=${ORBITS} CFLAGS='-z muldefs' LDFLAGS='-z muldefs' \
&&  make \
&&  make install
# define binaries search path
ENV PATH=${GMTSAR}/bin:$PATH

# install PyGMTSAR
RUN pip3 install pygmtsar

# install interactive plot libraries missed in the base image
RUN pip3 install hvplot datashader geoviews

# switch user
USER    ${NB_UID}
WORKDIR "${HOME}"

# download example notebooks
RUN svn export https://github.com/mobigroup/gmtsar/trunk/notebooks
RUN rm -rf notebooks/README.md work
