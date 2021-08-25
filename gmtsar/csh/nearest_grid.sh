#!/bin/sh
# Replacement for /usr/local/GMTSAR/bin/nearest_grid
set -e

if [ "$#" -ne "2" -a "$#" -ne "3" ]
then
    echo "Usage: nearest_grid input.grd output.grd [search_radius_pixels]"
    exit
fi

input="$1"
output="$2"
radius="$3"

# set default search radius, pixels
if [ -z "$radius" ]
then
    radius=300
fi

gdal_fillnodata.py -md "$radius" "${input}" "${output}.tif"
gmt grdreformat "${output}.tif" "${output}"
bounds=$(gmt grdinfo -C "${input}" | awk 'BEGIN{OFS = "/"}{print ($2,$3,$4,$5)}')
gmt grdedit "${output}" -R"$bounds" -T
rm -f "${output}.tif"
