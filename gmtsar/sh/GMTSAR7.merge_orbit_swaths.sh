#!/bin/sh
# Alexey Pechnikov, Aug, 2021, https://github.com/mobigroup/gmtsar
# See page 16 "Run All Interferograms" in https://topex.ucsd.edu/gmtsar/tar/sentinel_time_series_3.pdf
# use after GMTSAR.run_orbit.sh and GMTSAR.unwrap_orbit.sh next
# ./GMTSAR.merge_orbit.sh /home/jupyter/GMTSAR asc F2 1
# TODO: this file is specific for F2 and merges F2 + F3 where F3 is empty
set -e

workdir="$1"
orbit="$2"
swathlist="$3"

cd "$workdir"
cd "$orbit"
# do it with a directory name to prevent a big damage on script fails
rm -fr merge/*
cd merge

# number of subswaths plus 1 (symbols in string F...)
swaths=$(echo -n "${swathlist}" | wc -c | tr -d ' ')

# process one subswath
for num in $(echo "$swathlist" | grep -o . | tail -n +2 | head -n 1)
do
    swath="F${num}"
    echo "Merging master ${orbit}/${swath}"

    # link the required files
    cp "../${swath}/batch_tops.config" .
    ln -f -s ../topo/dem.grd .

    # for sbas
    ln -f -s "../${swath}/intf.in" .
    ln -f -s "../${swath}/baseline_table.dat" .

    # missed step to build a list of interferogram directories in the intf_all directory of F2 (intflist)
    find "../${swath}/intf_all/" -mindepth 1 -maxdepth 1 -name '*' -type d | sed -s "s|^../${swath}/intf_all/||g" > intflist

    # define master image to use for merging below
    master=$(grep master_image "../${swath}/batch_tops.config" | sed -E 's/(.*)(S1_.*_F[123])(.*)/\2/')
    echo "Master image: $master"

done

# Linking vs merging for a single subswath only
if [ "${swaths}" = "2" ]
then
    echo "Linking single subswath ${swaths} to merging directory"

    # we do not need to perform merging here because only one subswath just linked as is
    cat intflist | xargs -I {} -n 1 ln -s "../${swath}/intf_all/{}" .

    # for latlon coordinates calculation
    # add the used filter to the merge directory for proj_ra2ll.csh to define output grid resolution (1/4 of filter size)
    head -n 1 intflist | xargs -I {} -n 1 sh -c 'ln -f -s {}/gauss_* .'
    ln -f -s ../${swath}/topo/trans.dat .

    # for SBAS
    masterprm=$(find -L ../${swath}/intf_all -name "${master}.PRM" | head -n 1)
    echo "Link master PRM: $masterprm"
    ln -f -s "${masterprm}" .
    masterled=$(find -L ../${swath}/intf_all -name "${master}.LED" | head -n 1)
    echo "Link master LED: $masterled"
    ln -f -s "${masterled}" .
fi

if [ "$swaths" -gt 2 ]
then
    echo "Merge subswaths $swathlist"
    # mode 0 is merging all 3 subswaths
    # mode 1 is for F1/F2
    # mode 2 is for F2/F3
    # note: F1/F3 are not connected and so unsopported for the merging
    if echo -n "$swathlist" | grep -q 123
    then
        mode=0
    else
        if echo -n "$swathlist" | grep -q 12
        then
            # mode 1 is for F1/F2
            mode=1
        else
            # mode 2 is for F2/F3
            mode=2
        fi
    fi

    create_merge_input.csh intflist .. "${mode}" > merge_list

    masterline=$(grep "${master}.PRM" merge_list | head -n 1)
    echo "Master line: $masterline"
    mv merge_list merge_list.orig
    echo "$masterline" > merge_list
    grep -v "$masterline" -- merge_list.orig >> merge_list
    rm -f merge_list.orig

    # Run merge_batch.csh
    # switch_land = 1  in batch_tops.config but landmask is not generating,so it is  added below
    merge_batch.csh merge_list batch_tops.config

    # for SBAS
    # here is just one LED file
    masterled=$(find . -mindepth 2 -name '*LED')
    masterdir=$(dirname "$masterled")
    find . -mindepth 2 -name '*LED' | xargs -n 1 -I {} ln -f -s {} .
    # here are multiple supermaster.PRM
    find ./2021019_2021043 -name 'supermaster.PRM' | xargs -n 1 -I {} ln -f -s {} .
fi

# Stacking Coherence Grids for SBAS preprocessing later
ls ???????_???????/corr.grd > corr.grd_list
stack.csh corr.grd_list 1 corr_stack.grd std.grd
# see corr_stack.pdf std.pdf

# notify user via Telegram
if [ -n "$TELEGRAM_TOKEN" ]
then
    curl \
        -F "chat_id=${TELEGRAM_CHAT_ID}" \
        -F caption="GMTSAR8 on ${TELEGRAM_SENDER}: Mean Correlation Stack for subswaths ${swathlist}" \
        -F document="@corr_stack.pdf" \
        "https://api.telegram.org/bot${TELEGRAM_TOKEN}/sendDocument"
    curl \
        -F "chat_id=${TELEGRAM_CHAT_ID}" \
        -F caption="GMTSAR8 on ${TELEGRAM_SENDER}: Std. Dev. Correlation Stack for subswaths ${swathlist}" \
        -F document="@std.pdf" \
        "https://api.telegram.org/bot${TELEGRAM_TOKEN}/sendDocument"
fi
proj_ra2ll.csh trans.dat corr_stack.grd corr_stack_ll.grd
proj_ra2ll.csh trans.dat std.grd std_ll.grd
if [ -n "$TELEGRAM_TOKEN" ]
then
    curl \
        -F "chat_id=${TELEGRAM_CHAT_ID}" \
        -F caption="GMTSAR8 on ${TELEGRAM_SENDER}: Mean Correlation Stack Grid for subswaths ${swathlist}" \
        -F document="@corr_stack_ll.grd" \
        "https://api.telegram.org/bot${TELEGRAM_TOKEN}/sendDocument"
fi
