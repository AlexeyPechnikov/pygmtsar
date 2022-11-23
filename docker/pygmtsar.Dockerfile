# https://jupyter-docker-stacks.readthedocs.io/en/latest/using/selecting.html
# https://github.com/jupyter/docker-stacks/blob/main/base-notebook/Dockerfile
FROM jupyter/scipy-notebook:ubuntu-22.04

USER root

# install GMTSAR dependencies and some helpful command-line tools
RUN apt-get -y update \
&&  apt-get -y install git gdal-bin libgdal-dev subversion curl jq \
&&  apt-get -y install csh autoconf make gfortran \
&&  apt-get -y install libtiff5-dev libhdf5-dev liblapack-dev libgmt-dev gmt-dcw gmt-gshhg gmt \
&&  apt-get -y install zip htop mc \
&&  apt-get clean && rm -rf /var/lib/apt/lists/*

# define installation and binaries search paths
ARG GMTSAR=/usr/local/GMTSAR
ARG ORBITS=/usr/local/orbits
ENV PATH=${GMTSAR}/bin:$PATH

# install GMTSAR from git
RUN cd $(dirname ${GMTSAR}) \
&&  git clone --branch master https://github.com/gmtsar/gmtsar GMTSAR \
&&  cd ${GMTSAR} \
&&  autoconf \
&&  ./configure --with-orbits-dir=${ORBITS} CFLAGS='-z muldefs' LDFLAGS='-z muldefs' \
&&  make \
&&  make install

# install PyGMTSAR and additional plot libraries
RUN pip3 install pygmtsar \
&&  pip3 install matplotlib seaborn hvplot datashader geoviews

# workaround for libgeos-dev issue
RUN ln -s /usr/lib/aarch64-linux-gnu/libgeos_c.so /opt/conda/lib/libgeos_c.so || true

# switch user
USER    ${NB_UID}
WORKDIR "${HOME}"

# download example notebooks and CI test scripts and cleanup work dir
RUN svn export https://github.com/mobigroup/gmtsar/trunk/notebooks \
&&  svn export https://github.com/mobigroup/gmtsar/trunk/tests \
&&  rm -rf notebooks/README.md work
