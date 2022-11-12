#!/bin/sh
set -e

if [ "$#" -eq 0 ]
then
    echo "Default launch saves downloaded dataset for future usage"
    force=0
elif [ "$#" -eq 1 ]
then
    echo "Force launch removes downloaded dataset to save disk space"
    force=1
else
    echo "Incorrect number of arguments"
    exit 1
fi

echo "Downloading iCloud datasets..."
cat > icloud_urls.txt <<- END
https://www.icloud.com/iclouddrive/060EcVsMZF5o5nDoyAxHZ-JgQ#S1A_IW_SLC__1SDV_20210929T162245_20210929T162312_039899_04B89A_651D.zip
https://www.icloud.com/iclouddrive/0bcQQkK2QDSiX76_y1O6zfd_Q#S1B_IW_SLC__1SDV_20210923T162204_20210923T162231_028828_0370C4_681D.zip
https://www.icloud.com/iclouddrive/085_X__5ZszrqcnP--SGCH-tw#S1AB_Crete_Earthquake_vs_geObsevatory.jpg
END
cut -f5 -d '/' icloud_urls.txt | cut -f1 -d '#' | sed -E 's/(.*)/{"shortGUIDs":[{"value":"\1"}]}/' > requests.txt
cut -f5 -d '/' icloud_urls.txt | cut -f2 -d '#' > fnames.txt
paste -d ' ' requests.txt fnames.txt  | while read LINE
do
    request=$(echo "$LINE"| cut -d ' ' -f1)
    fname=$(echo "$LINE"| cut -d ' ' -f2)
    echo "request: $request fname: $fname"
    url=$(curl -s 'https://ckdatabasews.icloud.com/database/1/com.apple.cloudkit/production/public/records/resolve' \
    --data-raw "${request}" --compressed | jq -r '.results[0].rootRecord.fields.fileContent.value.downloadURL')
    # continue downloading or miss downloading when the complete file already exists
    wget -q --show-progress --progress=bar:force:noscroll -c -O "${fname}" "${url}"
done
echo "Unpacking selected files only to data directory..."
mkdir -p data
unzip -j -n S1A_IW_SLC__1SDV_20210929T162245_20210929T162312_039899_04B89A_651D.zip \
    '*.SAFE/*/s1?-iw3-slc-vv-*-006.xml' -d data
unzip -j -n S1A_IW_SLC__1SDV_20210929T162245_20210929T162312_039899_04B89A_651D.zip \
    '*.SAFE/*/s1?-iw3-slc-vv-*-006.tiff' -d data
unzip -j -n S1B_IW_SLC__1SDV_20210923T162204_20210923T162231_028828_0370C4_681D.zip \
    '*.SAFE/*/s1?-iw3-slc-vv-*-006.xml' -d data
unzip -j -n S1B_IW_SLC__1SDV_20210923T162204_20210923T162231_028828_0370C4_681D.zip \
    '*.SAFE/*/s1?-iw3-slc-vv-*-006.tiff' -d data
# cleanup
if [ "$force" -eq 1 ]
then
    echo "Removing downloaded datasets to free disk space following 'force' command line argument'..."
    rm *.zip
fi
echo "Checking disk space and downloaded files"
df -h
ls -lh
ls -lh data
echo "Running Python test script..."
rm -fr *.jpg
python3 ./S1AB_2021_Crete_Earthquake.py && echo SUCCESS
