#!/bin/sh
set -e

if [ "$#" -eq 2 -o "$#" -eq 3 ]
then
    DATADIR="$1"
    PATTERN="$2"
else
    echo "Use as: echo \"icloud_links\" | $0 DATADIR PATTERN FORCE_CLEANUP_FLAG"
    exit 1
fi

cat - > "${DATADIR}_icloud_urls.txt"

cut -f5 -d '/' "${DATADIR}_icloud_urls.txt" | cut -f1 -d '#' | sed -E 's/(.*)/{"shortGUIDs":[{"value":"\1"}]}/' > "${DATADIR}_requests.txt"
cut -f5 -d '/' "${DATADIR}_icloud_urls.txt" | cut -f2 -d '#' > "${DATADIR}_fnames.txt"
rm -fr "$DATADIR"
mkdir -p "$DATADIR"
paste -d ' ' "${DATADIR}_requests.txt" "${DATADIR}_fnames.txt"  | while read LINE
do
    request=$(echo "$LINE"| cut -d ' ' -f1)
    fname=$(echo "$LINE"| cut -d ' ' -f2)
    echo "request: $request fname: $fname"
    url=$(curl -s 'https://ckdatabasews.icloud.com/database/1/com.apple.cloudkit/production/public/records/resolve' \
        --data-raw "${request}" --compressed | jq -r '.results[0].rootRecord.fields.fileContent.value.downloadURL')
    # continue downloading or miss downloading when the complete file already exists
    wget -q --show-progress --progress=bar:force:noscroll -c -O "${fname}" "${url}"
    echo "Unpacking selected files only to data directory..."
    unzip -j -n "${fname}" "$PATTERN" -d "$DATADIR"
    # cleanup
    if [ "$3" != "" ]
    then
        echo "Removing downloaded datasets to free disk space following 'force' command line argument'..."
        rm "${fname}"
    fi
done
# cleanup
rm "${DATADIR}_icloud_urls.txt" "${DATADIR}_requests.txt" "${DATADIR}_fnames.txt"

echo "Checking downloaded files"
ls -lh
ls -lh "$DATADIR"
