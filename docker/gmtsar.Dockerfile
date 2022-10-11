# https://jupyter-docker-stacks.readthedocs.io/en/latest/using/selecting.html
# https://github.com/jupyter/docker-stacks/blob/main/base-notebook/Dockerfile
FROM ubuntu:22.04
# unpack arbitrary files like to mans, etc.
RUN yes | unminimize
# fail on pipe errors
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# install GMTSAR dependencies plus rsync
RUN set -e \
&&  apt-get -y update \
&&  apt-get -y install git subversion curl rsync ghostscript \
&&  apt-get -y install csh autoconf make gfortran \
&&  apt-get -y install libtiff5-dev libhdf5-dev liblapack-dev libgmt-dev gmt-dcw gmt-gshhg gmt \
&&  apt-get clean && rm -rf /var/lib/apt/lists/*

# define installation paths
ARG GMTSAR=/usr/local/GMTSAR
ARG ORBITS=/usr/local/orbits

# install GMTSAR from git
RUN set -e \
&&  cd /usr/local \
&&  git clone --branch master https://github.com/gmtsar/gmtsar GMTSAR
RUN set -e \
&&  cd ${GMTSAR} \
&&  autoconf \
&&  ./configure --with-orbits-dir=${ORBITS} CFLAGS='-z muldefs' LDFLAGS='-z muldefs' \
&&  make \
&&  make install

# add script to download orbits by user
RUN echo "#!/bin/sh"                                  > ${GMTSAR}/bin/download_orbits.sh
RUN echo "# download orbit data for ERS and Envisat" >> ${GMTSAR}/bin/download_orbits.sh
RUN echo "curl -O https://topex.ucsd.edu/gmtsar/tar/ORBITS.tar && sudo mkdir -p "$ORBITS" && sudo tar xf ORBITS.tar -C "$ORBITS" && rm ORBITS.tar" >> ${GMTSAR}/bin/download_orbits.sh
# set execution permissions to the script
RUN chmod a+x ${GMTSAR}/bin/download_orbits.sh

# install packages required to user work in terminal
RUN set -e \
&&  apt-get -y update \
&&  apt-get -y upgrade \
&&  apt-get install -y vim nano perl wget tar zip man sudo adduser w3m chafa mc htop \
&&  apt-get clean && rm -rf /var/lib/apt/lists/*

# allow ImageMagick to process PDF (required for chafa)
RUN sed -i 's/rights="none" pattern="PDF"/rights="read | write" pattern="PDF"/' /etc/ImageMagick-*/policy.xml

# create unprivileged user
RUN set -e \
&&  useradd --create-home --shell /bin/bash ubuntu \
&&  usermod -aG sudo ubuntu \
&&  echo "ubuntu ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/ubuntu \
&&  chmod 044 /etc/sudoers.d/ubuntu

# open terminal for the unprivileged user
USER ubuntu:ubuntu
WORKDIR /home/ubuntu
# define binaries search path for unprivileged user only
ENV PATH=${GMTSAR}/bin:$PATH
CMD ["/bin/bash"]
