#!/bin/csh -f
#
#  Last updated July 13, 2021 by KG
#  Based on script originally written by X. Xu
#
# ----------------------------------------------
#  This script reads in GNSS displacement data in
#  lon | lat | East | North | Up 
#  and converts it to:
#  range | azimuth | LOS(mm) 
#
  if ($#argv != 4 ) then
errormessage:
    echo ""
    echo "Usage: gnss_enu2los.csh master.PRM master.LED gnss.sllenu dem.grd"
    echo ""
    echo " master.PRM        -  PRM file for the master SAR image"
    echo " master.LED*       -  LED file for the master SAR image"
    echo " gnss.llenu**      -  GNSS displacements (Stat | Lon | Lat | E | N | U) in millimeters"
    echo " dem.grd           -  DEM grid file from your insar processing"
    echo " "
    echo " *Ensure you have the correct LED files listed in the master.PRM file "
    echo "  -- they must match "
    echo " "
    echo " **Assumes the gnss.llenu file is a list of stations with a specific" 
    echo " displacement value (e.g. the displacement between two insar scenes)"
    echo " per station."
    echo " "
    echo "Example: gnss_enu2los.csh master.PRM master.LED gnss_2018-2019.sllenu dem.grd "
    echo ""
    echo "Note: Check out correct_intf_with_gnss.csh to correct your interferogram"
    echo "with GNSS data"
    echo ""
    exit 1
  endif
  echo "gnss_enu2los.csh"
#
# -----------------------------------
# SET VARIABLE NAMES FOR CLARITY
# -----------------------------------
#
set output = "gnss_los.rad"
set PRM = $1
set LED = $2
set gnssenu = $3
set dem = $4
#
# Check the existence of the PRM and LED files
if ( -f $PRM ) then
     echo "PRM file exists..."
else
     echo "PRM does not exist -- this is required "
     exit 1
endif
if ( -f $LED ) then
     echo "LED file exists..."
else
     echo "LED does not exist -- this is required "
     exit 1
endif
#
#
# -----------------------------------
# + EXTRACT ELLIPSOIDAL HEIGHT FROM DEM
# + COMPUTE SATELLITE LOOK DIRECTION AT EACH STATION
# -----------------------------------
#
# pull the lon/lat, pull the height, input into SAT_look to get look direction
 awk '{print $2,$3}' $gnssenu | gmt grdtrack -G$dem -N > tmp.llh 
 SAT_look $PRM < tmp.llh > tmp.lltn
#
# -----------------------------------
# + PERFORM DOT PRODUCT TO CALCULATE LOS DISPLACEMENTS
# -----------------------------------
# perform the dot product with enu displacements and enu look direction
# This creates tmp.lltd which is lon | lat | dem height | LOS          
 paste $gnssenu tmp.lltn | awk '{if ($9 != nan) printf("%.9f %.9f %.9f %.12f\n",$2,$3,$9, ($4*$10)+($5*$11)+($6*$12))}' > tmp.lltd
 echo "LOS displacements calculated!"
#
# -----------------------------------
# CONVERT LON/LAT to RANGE/AZIMUTH
# -----------------------------------
# 
 SAT_llt2rat $PRM 1 < tmp.llh > tmp.ratll
 paste tmp.ratll tmp.lltd | awk '{print $1,$2,$9}' > $output 
 echo "Result created and stored as $output"
#
#
# ---------
# clean up
# ---------
#
 rm *tmp*
