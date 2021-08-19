#!/bin/sh
# intflist contains a list of all date1_date2 directories.
set -e

cd "$1"

if [ -z "${GMTSAR_THRESHOLD_SNAPHU_MASKED}" ]
then
    thr=0.001
else
    thr="${GMTSAR_THRESHOLD_SNAPHU_MASKED}"
fi
if [ -z "${GMTSAR_DEFOMAX}" ]
then
    jumps=0
else
    jumps="${GMTSAR_DEFOMAX}"
fi

ln -s ../mask_def.grd .
snaphu_interp.csh "$thr" "$jumps"
cd ..
