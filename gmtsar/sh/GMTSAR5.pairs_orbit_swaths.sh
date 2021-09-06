#!/bin/sh
# Alexey Pechnikov, Aug, 2021, https://github.com/mobigroup/gmtsar
# See "Aligning Secondary Images to Reference Image" in https://topex.ucsd.edu/gmtsar/tar/sentinel_time_series_3.pdf
# check this command for F1 to find the right parameters first, see pages 14-15
# select_pairs.csh baseline_table.dat 40 120
# use after GMTSAR.master_orbit.sh and GMTSAR.run_orbit.sh next
# ./GMTSAR.pairs_orbit.sh /home/jupyter/GMTSAR asc F2 40 120
set -e

workdir="$1"
orbit="$2"
swathlist="$3"
days="$4"
meters="$5"

cd "$workdir"
cd "$orbit"

masterswath=""
for num in $(echo "$swathlist" | grep -o . | tail -n +2)
do
    swath="F${num}"
    if [ -z "$masterswath" ]
    then
        # really configure the subswath
        masterswath="$swath"
        echo "Configure pairs in ${orbit}/${swath}"
        cd "${swath}"
        select_pairs.csh baseline_table.dat "$days" "$meters"

        # notify user via Telegram
        if [ -n "$TELEGRAM_TOKEN" ]
        then
            # paired baseline plot placed in subswath directory
            ps2pdf baseline.ps baseline.pdf
            curl \
                -F "chat_id=${TELEGRAM_CHAT_ID}" \
                -F caption="GMTSAR5 on ${TELEGRAM_SENDER}: Baseline pairs for ${days} days and ${meters} meters for ${swath}" \
                -F document='@baseline.pdf' \
                https://api.telegram.org/bot${TELEGRAM_TOKEN}/sendDocument
        fi

        cd ..
    else
        # copy and edit the intf.in file to ensure it refers to the F2 and F3 directories respectively (instead of F1)
        echo "Copy configured pairs from ${orbit}/${masterswath} to ${orbit}/${swath}"
        cd "${swath}"
        # Copy the intf.in file from the F1 directory to the F2 and F3 directories
        cp "../${masterswath}/intf.in" .
        sed -i -s "s/${masterswath}/${swath}/g" intf.in
        cd ..
    fi
done
