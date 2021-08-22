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

# make symlink to file in directory above if it is not exists
# for linked directories we need to copy the mask to here before
if [ ! -f mask_def.grd ]
then
    ln -f -s ../mask_def.grd .
fi
echo snaphu_interp.csh "$thr" "$jumps"
snaphu_interp.csh "$thr" "$jumps"
cd ..
