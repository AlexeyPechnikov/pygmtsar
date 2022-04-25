#!/bin/sh
# Alexey Pechnikov, Aug, 2021, https://github.com/mobigroup/gmtsar
# See page 16 "Run All Interferograms" in https://topex.ucsd.edu/gmtsar/tar/sentinel_time_series_3.pdf
# use after GMTSAR.pairs_orbit.sh and GMTSAR.merge_orbit.sh next
# ./GMTSAR.run_orbit.sh /home/jupyter/GMTSAR asc F2
set -e

workdir="$1"
orbit="$2"
swathlist="$3"

cd "$workdir"
cd "$orbit"

# define server cores count (on MacOS lscpu tool is missed)
# see also getconf _NPROCESSORS_ONLN
if which lscpu
then
    CPU_CORES=$(lscpu -p=CORE,ONLINE | grep -c 'Y')
else
    # Note: see also sysctl hw.logicalcpu
    CPU_CORES=$(sysctl -n hw.physicalcpu)
fi

for num in $(echo "$swathlist" | grep -o . | tail -n +2)
do
    swath="F${num}"
    echo "Run in ${orbit}/${swath}"
    cd "$swath"

    # cleanup from previous runs
    rm -fr intf intf_all intf_in/*
    # Run One Interferogram to Test Settings
    sed -i -s "s/^proc_stage.*\$/proc_stage = 1/g" batch_tops.config
    head -1 intf.in > one.in
    # cleanup from previous run
    rm -fr intf intf_all intf_in/*
    intf_tops.csh one.in batch_tops.config
    # Run All Interferograms
    sed -i -s "s/^proc_stage.*\$/proc_stage = 2/g" batch_tops.config
    intf_tops_parallel.csh intf.in batch_tops.config "$CPU_CORES"

    # montage phasefilt images grid
    cmd="montage"
    for fname in intf_all/???????_???????/phasefilt.pdf
    do
        name=$(dirname "${fname}")
        name=$(basename "${name}")
        cmd="${cmd} -label ${name} ${fname}"
    done
    cmd="${cmd} -pointsize 3 -geometry 300x300>+4+3 -density 600 -quality 90 phasefilt_grid.pdf"
    $cmd

    # notify user via Telegram
    telegram_senddocument.sh "Filtered phase images grid for ${swath}" phasefilt_grid.pdf
    cd ..
done
