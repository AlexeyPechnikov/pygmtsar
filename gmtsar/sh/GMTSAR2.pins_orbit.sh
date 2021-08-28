#!/bin/sh
# See "Adjusting Region of Interest" in https://topex.ucsd.edu/gmtsar/tar/sentinel_time_series_3.pdf
# use after GMTSAR.data_orbit.sh and GMTSAR.baseline_orbit.sh next
# ./GMTSAR.pins_orbit.sh /home/jupyter/GMTSAR asc 48.1 37.9 48.0 38.15
# see https://github.com/gmtsar/gmtsar/issues/180
set -e

workdir="$1"
orbit="$2"
lon1="$3"
lat1="$4"
lon2="$5"
lat2="$6"

cd "$workdir"
cd "$2"
cd reframed

cat << EOF > pins.ll
$lon1 $lat1
$lon2 $lat2
EOF

# return
cd ..

cd data
# list full paths to all of your .SAFE files in order of acquisition date:
ls -d $PWD/*.SAFE > SAFE_filelist
# data directory named Fxxxx_Fxxxx create and populate
# with new stitched/trimmed *SAFE data
organize_files_tops_linux.csh SAFE_filelist ../reframed/pins.ll 2
