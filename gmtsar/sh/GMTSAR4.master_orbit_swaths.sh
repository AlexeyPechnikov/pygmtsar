#!/bin/sh
# Alexey Pechnikov, Aug, 2021, https://github.com/mobigroup/gmtsar
# See "Aligning Secondary Images to Reference Image" in https://topex.ucsd.edu/gmtsar/tar/sentinel_time_series_3.pdf
# use after GMTSAR.baseline_orbit.sh and GMTSAR.pairs_orbit.sh next
# ./GMTSAR.master_orbit.sh /home/jupyter/GMTSAR asc F2 20210402
set -e

workdir="$1"
orbit="$2"
swathlist="$3"
master="$4"

cd "$workdir"
cd "$orbit"

for num in $(echo "$swathlist" | grep -o . | tail -n +2)
do
    swath="F${num}"
    echo "Reset master image in ${orbit}/${swath}/batch_tops.config"
    cd "${swath}"
    sed -i -s "s/^master_image = .*\$/master_image = S1_${master}_ALL_${swath}/g" batch_tops.config
    cd ..
    if [ -f "${swath}/raw/data.in" ]
    then
        echo "Processing ${orbit}/${swath}/data.in"
        cd "${swath}/raw"
        scene=$(grep "$master" -- data.in | head -n1)
        echo "master scene: $scene"
        mv data.in data.in.orig
        echo "$scene" > data.in
        grep -v "$scene" -- data.in.orig >> data.in
        # Re-run preproc_batch_tops.csh
        preproc_batch_tops.csh data.in dem.grd 2 | tee preproc_batch_tops.log
        # return
        cd ../..
    else
        echo "Ignore not existing ${orbit}/${swath}/data.in"
    fi
done
