#!/bin/sh
# See "Data Selection and Download for Region of Interest" in https://topex.ucsd.edu/gmtsar/tar/sentinel_time_series_3.pdf
# Supported as zipped Sentinel-1 SLC products (.zip) as unpacked ones (.SAFE)
# use after GMTSAR.setup.sh as
# ./GMTSAR.data_orbit.sh /home/jupyter/GMTSAR asc /home/jupyter/YamchiDamAsc
# or
# ./GMTSAR.data_orbit.sh /home/jupyter/GMTSAR asc /home/jupyter/YamchiDamAsc *.zip
set -e

workdir="$1"
orbit="$2"
datadir="$3"
pattern="$4"

cd "$workdir/$orbit/data/"

# define server core count
CPU_CORES=$(lscpu -p=CORE,ONLINE | grep -c 'Y')

# no pattern defined
if [ "$#" = "3" ]
then
    pattern="S1A_IW_SLC_*.zip"
fi
if echo -n "${pattern}" | grep -q ".zip"
then
    echo "Copy zipped data by pattern ${datadir}/${pattern}"
    find "$datadir" -name "${pattern}" -print0 | xargs -0 -I {} -n 1 -P ${CPU_CORES} unzip {}
    basepattern=$(basename "$pattern" .zip)
    find "$datadir" -name "${basepattern}.EOF" -print0 | xargs -0 -I {} -n 1 -P ${CPU_CORES} ln -f -s {} .
fi

# no pattern defined
if [ "$#" = "3" ]
then
    pattern="S1A_IW_SLC_*.SAFE"
fi
if echo -n "${pattern}" | grep -q ".SAFE"
then
    echo "Link unzipped data by pattern ${datadir}/${pattern}"
    find "$datadir" -name "$pattern" -print0 | xargs -0 -I {} -n 1 -P ${CPU_CORES} ln -f -s {} .
    basepattern=$(basename "$pattern" .SAFE)
    find "$datadir" -name "${basepattern}.EOF" -print0 | xargs -0 -I {} -n 1 -P ${CPU_CORES} ln -f -s {} .
fi

# download orbit files
eof
