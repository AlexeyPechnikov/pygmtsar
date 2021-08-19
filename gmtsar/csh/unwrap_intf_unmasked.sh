#!/bin/sh
# intflist contains a list of all date1_date2 directories.
set -e

cd "$1"
if [ -z "${GMTSAR_THRESHOLD_SNAPHU_UNMASKED}" ]
then
    thr=0.10
else
    thr="${GMTSAR_THRESHOLD_SNAPHU_UNMASKED}"
fi
if [ -z "${GMTSAR_DEFOMAX}" ]
then
    jumps=0
else
    jumps="${GMTSAR_DEFOMAX}"
fi
snaphu_interp.csh "$thr" "$jumps"
cd ..
