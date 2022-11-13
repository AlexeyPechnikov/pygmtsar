#!/bin/sh
set -e

DATADIR=data_crete
PATTERN="*.SAFE/*/s1?-iw3-slc-vv-*"

echo "Downloading iCloud datasets..."
echo "\
https://www.icloud.com/iclouddrive/060EcVsMZF5o5nDoyAxHZ-JgQ#S1A_IW_SLC__1SDV_20210929T162245_20210929T162312_039899_04B89A_651D.zip
https://www.icloud.com/iclouddrive/0bcQQkK2QDSiX76_y1O6zfd_Q#S1B_IW_SLC__1SDV_20210923T162204_20210923T162231_028828_0370C4_681D.zip" \
| ./icloud_download.sh "$DATADIR" "$PATTERN" "$1"

echo "Checking disk space..."
df -h

echo "Running Python test script..."
rm -fr *.jpg
python3 ./S1AB_2021_Crete_Earthquake.py && echo SUCCESS
