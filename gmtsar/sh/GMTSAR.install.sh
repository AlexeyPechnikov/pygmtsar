#!/bin/sh
# use the script first to prepare the linux system and GMTSAR.setup.sh next
set -e

# prepare system
apt install -y locales
# Uncomment en_US.UTF-8 for inclusion in generation
sed -i 's/^# *\(en_US.UTF-8\)/\1/' /etc/locale.gen
# Generate locale
#locale-gen
locale-gen en_US.UTF-8
# check generated locales
#locale -a
# Export env vars
echo "export LC_ALL=en_US.UTF-8" >> ~/.bashrc
echo "export LANG=en_US.UTF-8" >> ~/.bashrc
echo "export LANGUAGE=en_US.UTF-8" >> ~/.bashrc

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
# for python tools and sbas command line arguments calculation
apt install -y python3-pip gdal-bin python-gdal python3-netcdf4 python3-scipy bc
# for my additional Python-coded utilities
python3 -m pip install xarray numpy scipy --upgrade

mkdir -p /usr/local/orbits
cd /usr/local/orbits
#wget --content-disposition https://hu.berlin/s1-orbits
#tar -xzvf S1-Orbits.tar.gz

cd /usr/local
# drop previous installation if needed
rm -fr GMTSAR
#git clone --branch 6.1 https://github.com/gmtsar/gmtsar GMTSAR
git clone --branch master https://github.com/mobigroup/gmtsar GMTSAR
cd GMTSAR
autoconf
./configure --with-orbits-dir=/usr/local/orbits
make
make install

# replace original binary tool by Pythonic one
mv /usr/local/GMTSAR/bin/nearest_grid /usr/local/GMTSAR/bin/nearest_grid.orig
mv /usr/local/GMTSAR/bin/nearest_grid.py /usr/local/GMTSAR/bin/nearest_grid

# configure binaries path
cat << EOF >> ~/.bashrc

export GMTSAR=/usr/local/GMTSAR
export PATH=\$GMTSAR/bin:"\$PATH"

EOF

# to download Sentinel-1 orbit files and SRTM 30m and 90m DEM
python3 -m pip install sentineleof elevation
eio selfcheck
# https://topex.ucsd.edu/gmtsar/tar/sentinel_time_series_3.pdf
# see page 10
#echo "alias wgetasf='echo NOTICE: Call to wgetasf shell alias for orbit files download. Use command eof --path . instead'" >> ~/.bashrc

# for GMTSAR CSH scripts
touch ~/.cshrc

# load the configuration
. ~/.bashrc

# allow ImageMagick to process PS and PDF files (it was already fixed vulnerability in Ghostscript)
sed -i '/policy domain="coder" rights="none" pattern="PS"/d'  /etc/ImageMagick-6/policy.xml
sed -i '/policy domain="coder" rights="none" pattern="PDF"/d' /etc/ImageMagick-6/policy.xml
