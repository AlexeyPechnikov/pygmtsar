#!/bin/sh
# Alexey Pechnikov, Aug, 2021, https://github.com/mobigroup/gmtsar

##########################################################################################
# Set ENV variables
##########################################################################################
#export TELEGRAM_CHAT_ID=
#export TELEGRAM_TOKEN=
export TELEGRAM_SENDER=$(hostname)

# recommended range decimation to be 8, azimuth decimation to be 2 (use 4/1, 8/2, 16/4))
export GMTSAR_RANGE_DECIMATION=4
export GMTSAR_AZIMUTH_DECIMATION=1
# set to number of cycles to allow phase jumps
export GMTSAR_DEFOMAX=0
# threshold for well correlated not masked areas (default 0.10)
export GMTSAR_THRESHOLD_SNAPHU=0.12

##########################################################################################
# Set local variables
##########################################################################################
DATADIR=/data/YamchiDamAsc
# select unpacked Sentinel-1 files for 2021 year
DATAPATTERN="S1A*_2021*.SAFE"
WORKDIR=/mnt/GMTSAR
#DEM_BOUNDS="46.4 37.3 49.9 39.6"
DEM=/data/srtm3_egm96.grd
ORBIT=asc
SUBSWATHS=F2
# https://github.com/gmtsar/gmtsar/issues/180
# does not work asc 48.1 37.9 48.0 38.15
# large area - ok asc 48.2 37.8 47.9 38.2
# small area - ok asc 48.1 37.9 48.0 38.2
#export PINS="48.1 37.9 48.0 38.4"
PINS="48.1 37.8 48 38.8"
MASTERDATE=20210213
BASEDAYS=75
BASEMETERS=75
# define to use correlation mask (comment it to disable)
#USEMASK=1
# define to enable landmask (comment it to disable)
USELANDMASK=1
# threshold for poor correlated masked areas (default 0.001)
#export GMTSAR_THRESHOLD_SNAPHU_MASKED=0.001
#-smooth sf           --  smoothing factors, default=0
#-atm ni              --  number of iterations for atmospheric correction, default=0(skip atm correction) 
SBASOPTS="-atm 0 -smooth 5.0"

# notify user via Telegram
if [ -n "$TELEGRAM_TOKEN" ]
then
    curl \
        -F "chat_id=${TELEGRAM_CHAT_ID}" \
        -F caption="GMTSAR on ${TELEGRAM_SENDER}: Start" \
        -F document="@$0" \
        "https://api.telegram.org/bot${TELEGRAM_TOKEN}/sendDocument"
fi
##########################################################################################
# Start GMTSAR processing
##########################################################################################
rm -f GMTSAR.log
GMTSAR0.setup.sh                 "$WORKDIR" "$DEM"                                          2>&1 | tee -a GMTSAR.log
GMTSAR1.data_orbit.sh            "$WORKDIR" "$ORBIT" "$DATADIR" "$DATAPATTERN"              2>&1 | tee -a GMTSAR.log
GMTSAR2.pins_orbit.sh            "$WORKDIR" "$ORBIT"  $PINS                                 2>&1 | tee -a GMTSAR.log
GMTSAR3.baseline_orbit_swaths.sh "$WORKDIR" "$ORBIT" "$SUBSWATHS"                           2>&1 | tee -a GMTSAR.log
GMTSAR4.master_orbit_swaths.sh   "$WORKDIR" "$ORBIT" "$SUBSWATHS" "$MASTERDATE"             2>&1 | tee -a GMTSAR.log
GMTSAR5.pairs_orbit_swaths.sh    "$WORKDIR" "$ORBIT" "$SUBSWATHS" "$BASEDAYS" "$BASEMETERS" 2>&1 | tee -a GMTSAR.log
GMTSAR6.run_orbit_swaths.sh      "$WORKDIR" "$ORBIT" "$SUBSWATHS"                           2>&1 | tee -a GMTSAR.log
GMTSAR7.merge_orbit_swaths.sh    "$WORKDIR" "$ORBIT" "$SUBSWATHS"                           2>&1 | tee -a GMTSAR.log
GMTSAR8.unwrap_orbit.sh          "$WORKDIR" "$ORBIT" "$USEMASK" "$USELANDMASK"              2>&1 | tee -a GMTSAR.log
GMTSAR9.sbas_orbit.sh            "$WORKDIR" "$ORBIT" "$SBASOPTS"                            2>&1 | tee -a GMTSAR.log
##########################################################################################
# End GMTSAR processing
# see /mnt/GMTSAR/asc/SBAS for disp_*_ll.grd and vel_ll.grd & vel_ll.kml & vel_ll.png
##########################################################################################
# notify user via Telegram
if [ -n "$TELEGRAM_TOKEN" ]
then
    curl \
        -F "chat_id=${TELEGRAM_CHAT_ID}" \
        -F caption="GMTSAR on ${TELEGRAM_SENDER}: End" \
        -F document="@GMTSAR.log" \
        "https://api.telegram.org/bot${TELEGRAM_TOKEN}/sendDocument"
fi
