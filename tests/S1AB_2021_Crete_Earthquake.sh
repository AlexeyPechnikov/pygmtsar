#!/bin/sh
set -e

echo "Downloading iCloud datasets..."
./icloud_download.sh \
    "https://www.icloud.com/iclouddrive/060EcVsMZF5o5nDoyAxHZ-JgQ#S1A_IW_SLC__1SDV_20210929T162245_20210929T162312_039899_04B89A_651D.zip" \
    "https://www.icloud.com/iclouddrive/0bcQQkK2QDSiX76_y1O6zfd_Q#S1B_IW_SLC__1SDV_20210923T162204_20210923T162231_028828_0370C4_681D.zip"
find . -name '*.zip' -exec unzip -j -n '{}' '*.SAFE/*/s1?-iw3-slc-vv-*' -d data_crete \;
if [ "$1" != "" ]
then
    rm *.zip
fi

echo "Checking disk space..."
df -h

echo "Running Python test script..."
rm -fr *.jpg
python3 ./S1AB_2021_Crete_Earthquake.py && echo SUCCESS
