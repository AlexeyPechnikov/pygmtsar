#!/bin/sh
# Use as: echo \"icloud_links\" | icloud_download.sh
set -e

for URL in "$@"
do
    request=$(echo "$URL" | cut -f5 -d '/' | cut -f1 -d '#' | sed -E 's/(.*)/{"shortGUIDs":[{"value":"\1"}]}/' | cut -d ' ' -f1)
    fname=$(echo "$URL" | cut -f5 -d '/' | cut -f2 -d '#' | cut -d ' ' -f2)
    echo "request: $request fname: $fname"
    url=$(curl -s 'https://ckdatabasews.icloud.com/database/1/com.apple.cloudkit/production/public/records/resolve' \
        --data-raw "${request}" --compressed | jq -r '.results[0].rootRecord.fields.fileContent.value.downloadURL')
    # continue downloading or miss downloading when the complete file already exists
    wget -q --show-progress --progress=bar:force:noscroll -c -O "${fname}" "${url}"
done
