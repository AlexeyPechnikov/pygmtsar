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
https://www.icloud.com/iclouddrive/04946wsi9N5g0W6beTl_KUvQg#S1A_IW_SLC__1SSV_20160408T091355_20160408T091430_010728_01001F_83EB.zip
https://www.icloud.com/iclouddrive/059EHU1iazg_CDEBWd7CPjXnQ#S1A_IW_SLC__1SSV_20160420T091355_20160420T091423_010903_010569_F9CE.zip
https://www.icloud.com/iclouddrive/042d2NP0tDVFz_N39jCu5yyVA#S1A_Kumamoto_Earthquake_ESA_Sentinel_1_Toolbox.jpg
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
    wget -c -O "${fname}" "${url}"
done
echo "Unpacking selected files only to data directory..."
mkdir -p data
unzip -j -n S1A_IW_SLC__1SSV_20160408T091355_20160408T091430_010728_01001F_83EB.zip \
    '*.SAFE/*/s1?-iw1-slc-vv-*-001.xml' -d data
unzip -j -n S1A_IW_SLC__1SSV_20160408T091355_20160408T091430_010728_01001F_83EB.zip \
    '*.SAFE/*/s1?-iw1-slc-vv-*-001.tiff' -d data
unzip -j -n S1A_IW_SLC__1SSV_20160420T091355_20160420T091423_010903_010569_F9CE.zip \
    '*.SAFE/*/s1?-iw1-slc-vv-*-001.xml' -d data
unzip -j -n S1A_IW_SLC__1SSV_20160420T091355_20160420T091423_010903_010569_F9CE.zip \
    '*.SAFE/*/s1?-iw1-slc-vv-*-001.tiff' -d data
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
python3 ./S1A_2016_Kumamoto_Earthquake.py && echo SUCCESS
