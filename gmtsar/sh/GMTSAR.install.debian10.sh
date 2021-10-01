#!/bin/sh
# Alexey Pechnikov, Aug, 2021, https://github.com/mobigroup/gmtsar
# see https://topex.ucsd.edu/gmtsar/tar/sentinel_time_series_3.pdf
# the script install GMTSAR from my repository and add plus replaces some tools
# use the script first to prepare the linux system and GMTSAR.setup.sh next
# check init lag at /var/log/daemon.log (/var/log/syslog on Ubuntu)
set -e

ORBITS=/usr/local/orbits
GMTSAR=/usr/local/GMTSAR
GIT=https://github.com/mobigroup/gmtsar
BRANCH=master
#GIT=https://github.com/gmtsar/gmtsar
#BRANCH=6.1

# prepare system
apt update
apt -y upgrade
apt install -y locales
# Uncomment en_US.UTF-8 for inclusion in generation
sed -i 's/^# *\(en_US.UTF-8\)/\1/' /etc/locale.gen
# Generate locale
#locale-gen
locale-gen en_US.UTF-8
# check generated locales
#locale -a
# Export env vars
echo "export LC_ALL=en_US.UTF-8" >> /root/.bashrc
echo "export LANG=en_US.UTF-8" >> /root/.bashrc
echo "export LANGUAGE=en_US.UTF-8" >> /root/.bashrc

# https://github.com/gmtsar/gmtsar/wiki/GMTSAR-Wiki-Page
apt install -y csh subversion autoconf libtiff5-dev libhdf5-dev wget
apt install -y liblapack-dev
apt install -y gfortran
apt install -y g++
apt install -y libgmt-dev
apt install -y gmt-dcw gmt-gshhg
# gmt-gshhg-full should be installed automatically (it is required to use GMTSAR landmask)
apt install -y gmt
# fix for missed "gs" utility and git and make tools
apt install -y imagemagick git make
# https://topex.ucsd.edu/gmtsar/tar/sentinel_time_series_3.pdf
# see page 16
apt install -y parallel

mkdir -p "$ORBITS"
cd "$ORBITS"
#wget --content-disposition https://hu.berlin/s1-orbits
#tar -xzvf S1-Orbits.tar.gz


cd $(dirname "$GMTSAR")
# drop previous installation if needed
rm -fr GMTSAR
git clone --branch "$BRANCH" "$GIT" GMTSAR
cd GMTSAR
autoconf
./configure --with-orbits-dir="$ORBITS"
make
make install

# replace original binary tool by Pythonic one
mv "${GMTSAR}/bin/nearest_grid"          "${GMTSAR}/bin/nearest_grid.orig"
mv "${GMTSAR}/gmtsar/py/nearest_grid.py" "${GMTSAR}/bin/nearest_grid"
# install additional scripts from the repo
find "${GMTSAR}/gmtsar/sh" -name '*.sh' -print0 | xargs -0 -I {} -n 1 ln -f -s "{}" "${GMTSAR}/bin/"

# configure binaries path
cat << EOF >> /root/.bashrc

export GMTSAR="${GMTSAR}"
export PATH=\$GMTSAR/bin:"\$PATH"

EOF

# to download Sentinel-1 orbit files and SRTM 30m and 90m DEM
apt install -y unzip
python3 -m pip install sentineleof elevation
eio selfcheck

# for GMTSAR CSH scripts
touch /root/.cshrc

# allow ImageMagick to process PS and PDF files (it was already fixed vulnerability in Ghostscript)
sed -i '/policy domain="coder" rights="none" pattern="PS"/d'  /etc/ImageMagick-6/policy.xml
sed -i '/policy domain="coder" rights="none" pattern="PDF"/d' /etc/ImageMagick-6/policy.xml

# PyGMTSAR dependencies
# for python tools and sbas command line arguments calculation
apt install -y python3-pip gdal-bin python-gdal python3-netcdf4 python3-scipy bc
# for additional Python-coded utilities
python3 -m pip install xarray numpy scipy pytest --upgrade
# for notebooks
apt -y install libgmt5 netcdf-bin python3-pip python3-netcdf4 python3-scipy python3-dev
apt -y install gdal-bin python3-gdal libcharls2 libgdal-dev libproj-dev proj-data proj-bin libgeos-dev
