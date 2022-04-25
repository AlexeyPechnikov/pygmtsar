#!/bin/sh
# Alexey Pechnikov, Aug, 2021, https://github.com/mobigroup/gmtsar
# See page 18 "Run All Interferograms" in https://topex.ucsd.edu/gmtsar/tar/sentinel_time_series_3.pdf
# use after GMTSAR.merge_orbit.sh and GMTSAR.sbas_orbit.sh next
# ./GMTSAR.unwrap_orbit.sh /mnt/GMTSAR asc
# or
# ./GMTSAR.unwrap_orbit.sh /mnt/GMTSAR asc 0.15
set -e

workdir="$1"
orbit="$2"
usemask="$3"
uselandmask="$4"

cd "$workdir"
cd "$orbit"
cd merge

# cleanup for multiple runs
rm -f */*_patch.grd */*_interp.grd */unwrap* */landmask* */tmp* */conncomp.grd */gmt.history
rm -f landmask* log_???????_???????.txt unwrap* tmp* gmt.history

# define server cores count (on MacOS lscpu tool is missed)
# see also getconf _NPROCESSORS_ONLN
if which lscpu
then
    CPU_CORES=$(lscpu -p=CORE,ONLINE | grep -c 'Y')
else
    # Note: see also sysctl hw.logicalcpu
    CPU_CORES=$(sysctl -n hw.physicalcpu)
fi

# create landmask if needed
if [ -n "${uselandmask}" ]
then
    echo "Create landmask"
    region=$(head -n 1 intflist | xargs -I {} -n 1 gmt grdinfo -C {}/phasefilt.grd | cut -s -f2,3,4,5 --output-delimiter '/')
    # xmin/xmax/ymin/ymax
    echo landmask.csh "region_cut$region"
    landmask.csh "region_cut$region"
    # because these could be linked directories (for a single subswath processing) we just copy files instead of symlinking
    cat intflist | xargs -I {} -n 1 cp -f landmask_ra.grd {}/landmask_ra.grd
fi

if [ -z "${usemask}" ]
then
    # a. Unwrapping in Regions of Good Coherence
    echo "Unwrap without mask for Regions of Good Coherence"
    unwrap_parallel.csh unwrap_intf_unmasked.sh intflist "$CPU_CORES"
else
    # b. Unwrapping in Regions of Poor Coherence
    echo "Unwrap with mask for Regions of Poor Coherence"
    # build mask
    gmt grdmath corr_stack.grd "${GMTSAR_THRESHOLD_SNAPHU}" GE 0 NAN = mask_def.grd
    # copy mask because symlinks in linked directories can not work
    cat intflist | xargs -I {} -P "${CPU_CORES}" -n 1 cp --remove-destination mask_def.grd "{}/"
    # Projecting into Latitude/Longitude
    proj_ra2ll.csh trans.dat mask_def.grd mask_def_ll.grd
    gmt grd2cpt mask_def.grd -T= -Z -Cjet > mask_def_ll.cpt
    grd2kml.csh mask_def_ll mask_def_ll.cpt
    # notify user via Telegram
    telegram_sendmediagroup.sh "Unwrap Mask" mask_def_ll.png mask_def_ll.kml
    unwrap_parallel.csh unwrap_intf_masked.sh intflist "$CPU_CORES"
fi

# montage unwrap images grid
cmd="montage"
for fname in ???????_???????/unwrap.pdf
do
    name=$(dirname "${fname}")
    name=$(basename "${name}")
    cmd="${cmd} -label ${name} ${fname}"
done
cmd="${cmd} -pointsize 3 -geometry 300x300>+4+3 -density 600 -quality 90 unwrap_grid.pdf"
$cmd

# notify user via Telegram
telegram_senddocument.sh "Unwrapped phase images grid" unwrap_grid.pdf
