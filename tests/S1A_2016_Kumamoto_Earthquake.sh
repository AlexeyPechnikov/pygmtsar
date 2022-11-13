#!/bin/sh
set -e

DATADIR=data_kumamoto
PATTERN="*.SAFE/*/s1?-iw1-slc-vv-*"

echo "Downloading iCloud datasets..."
echo "\
https://www.icloud.com/iclouddrive/04946wsi9N5g0W6beTl_KUvQg#S1A_IW_SLC__1SSV_20160408T091355_20160408T091430_010728_01001F_83EB.zip
https://www.icloud.com/iclouddrive/059EHU1iazg_CDEBWd7CPjXnQ#S1A_IW_SLC__1SSV_20160420T091355_20160420T091423_010903_010569_F9CE.zip" \
| ./icloud_download.sh "$DATADIR" "$PATTERN" "$1"

echo "Checking disk space..."
df -h

echo "Running Python test script..."
rm -fr *.jpg
python3 ./S1A_2016_Kumamoto_Earthquake.py && echo SUCCESS
