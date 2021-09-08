#!/bin/sh
# Alexey Pechnikov, Aug, 2021, https://github.com/mobigroup/gmtsar
# See "Selecting a Reference Image from a Baseline Plot" in https://topex.ucsd.edu/gmtsar/tar/sentinel_time_series_3.pdf
# use after GMTSAR.pins_orbit.sh and GMTSAR.master_orbit.sh next
# ./GMTSAR.baseline_orbit.sh /home/jupyter/GMTSAR asc F2
set -e

workdir="$1"
orbit="$2"
swathlist="$3"

cd "$workdir"
cd "$orbit"

for num in $(echo "$swathlist" | grep -o . | tail -n +2)
do
    swath="F${num}"
    echo "Prepare ${orbit}/${swath}/raw directory"
    cd "${swath}/raw"
    # Link the data and orbit files that belong to this subswath
    if [ $(find ../../data -name 'F????_F????' -type d) ]
        echo "Found pinned data files"
        ln -f -s ../../data/F????_F????/*.SAFE/*/*iw${num}*vv*xml .
        ln -f -s ../../data/F????_F????/*.SAFE/*/*iw${num}*vv*tiff .
    else
        echo "Found original data files"
        ln -f -s ../../data/*.SAFE/*/*iw${num}*vv*xml .
        ln -f -s ../../data/*.SAFE/*/*iw${num}*vv*tiff .
    fi
    ln -f -s ../../data/*EOF .
    ln -f -s ../topo/dem.grd .
    # prepare data.in where the data file names (no suffix) and orbit file names organized
    # for macos
    #prep_data.csh
    # for linux
    #prep_data_linux.csh
    prep_data.sh
    # generate a table of baselines
    preproc_batch_tops.csh data.in dem.grd 1
    # Save (move) baseline_table.dat for later use
    mv baseline_table.dat ../
    # check baseline.ps
    # return
    cd ../..
done

# notify user via Telegram
# goto last visited subswath
# dotted baseline plot placed in subswath/raw directory
cd "${swath}/raw"
ps2pdf baseline.ps baseline.pdf
telegram_senddocument.sh "Baseline plot for ${swath}" baseline.pdf
# return
cd ../..
