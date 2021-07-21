#!/bin/csh -f
#       $Id$
#
#  DTS - June 22, 2021
#  Edited Last - July 13, 2021 by KG
#
#  Correct the LOS unwrapped phase grid using GPS data. 
#  Everything is in range/azimuth coordinates.
#  Phase grid (interferogram) is LOS in radians 
#
#
  alias rm 'rm -f'
  gmt set IO_NC4_CHUNK_SIZE classic
#

  if ($#argv != 4 ) then
errormessage:
    echo ""
    echo "Usage: correct_insar_with_gnss.csh master.PRM phase.grd gnsslos.rad filter_wavelength"
    echo ""
    echo " master.PRM        -  PRM file for the master SAR image"
    echo " phase.grd         -  phase file to be corrected"
    echo " gnsslos.rad       -  GNSS displacements in LOS in millimeters"
    echo " filter_wavelength -  wavelength of the filter in meters (0.5 gain)"
    echo " "
    echo "Example: correct_insar_with_gnss.csh supermaster.PRM unwrap_dsamp.grd gnss_los.rad 40000 "
    echo ""
    echo "See gnss_enu2los.csh for converting GNSS ENU displacement to LOS"
    echo ""
    exit 1
  endif
  echo "correct_insar_with_gnss.csh"
#
# ----------------------
# SET VARIABLES      
# ----------------------
#
# compute the range (40 degree look angle) and azimuth pixel size in meters
#
  set output = "gnss_corrected_intf.grd"
  set prm = $1
  set insar = $2
  set gnss = $3
  set filterw = $4
  set rng_samp_rate = `grep rng_samp_rate $prm | awk 'NR == 1 {printf("%d", $3)}'`
  set rng_px_size = `grep rng_samp_rate $prm | awk 'NR == 1 {printf("%f", 1.556*299792458/$3/2)}'`
  set vel = `grep SC_vel $prm | awk 'NR == 1 {printf("%f", $3)}'`
  set PRF = `grep PRF $prm | awk 'NR == 1 {printf("%f", $3)}'`
  set azi_px_size = `echo 'scale=3;' $vel '/' $PRF | bc`
# echo $rng_px_size, $azi_px_size
#
#  now get the decimation factors
#
  set x_inc = `gmt grdinfo $insar | grep x_inc | awk 'NR == 1 {printf("%d", $7+.5)}'`
  set y_inc = `gmt grdinfo $insar | grep y_inc | awk 'NR == 1 {printf("%d", $7+.5)}'`
# echo $x_inc, $y_inc
#
# the actual pixel size is
#
  set dx = `echo 'scale=3;' $x_inc '*' $rng_px_size | bc `
  set dy = `echo 'scale=3;' $y_inc '*' $azi_px_size | bc `0
#  echo $dx, $dy
#
# ----------------------
# PREPARE THE GRID CHARACTERISTICS
# ----------------------
#
  echo "Preparing grid characteristics..."
#  downsample the grd-file to make a template for blockmedian and surface
#  make a sample size of the wavelength/16 
#
  set ndx = `echo  $filterw '/' $dx ' / 16' | bc`
  set ndy = `echo  $filterw '/' $dy ' / 16' | bc`
# echo $ndx, $ndy
#
#  subsample the phase grid to make a template for blockmedian and surface
#
  set ndx2 = `echo  $x_inc '*' $ndx | bc `
  set ndy2 = `echo  $y_inc '*' $ndy | bc `
# echo $ndx2,$ndy2
  gmt grdsample $insar  -I$ndx2/$ndy2 -Gtmp.grd
#
# ----------------------
# COMPUTE THE CORRECTION
# ----------------------
  echo "Computing the correction grid..."
# 1 km filter on phase grid to prevent sampling of bad pixels 
  set fsx = ` echo "1000" | awk -v xs="${dx}" '{printf "%d",$1/xs}' | awk '($1 % 2 == 0) {print $1+1} ($1 % 2 != 0) {print $1}' ` 
  set fsy = ` echo "1000" | awk -v ys="${dy}" '{printf "%d",$1/ys}' | awk '($1 % 2 == 0) {print $1+1} ($1 % 2 != 0) {print $1}' `
  gmt grdfilter $insar -Dp -Fg$fsx/$fsy -Gtmp_filt.grd
#
#  Calculate the difference between insar - [your input GNSS data in LOS and mm]
#  this converts gnss to radians, performs the difference, and removes nan values for stations not in your insar scene
#
  set wave = `grep wavelength $prm | awk '{print $3}'`
  awk '{print $1, $2, $3*(4.0*3.141592653)/(-'$wave'*1000.0) }' < $gnss | gmt grdtrack -Gtmp_filt.grd | awk '{print $1,$2,($4-$3)}' | grep -v 'nan' > ins-gps_diff.rad
#
#  now blockmedian and surface the data at this resolution
#
  gmt blockmedian ins-gps_diff.rad -Rtmp.grd > tmp.rad
  gmt surface tmp.rad -Rtmp.grd -Gtmp1.grd 
#
# ----------------------
# FILTER THE CORRECTION GRID
# ----------------------
#
  echo "Applying Gaussian filter to correction grid..."
  echo "...this may take ~1 minute..."
#
#  filter the grid using a Gaussian with wavelength in pixels
#  this length is 16 + 1 pixels
#
  gmt grdfilter tmp1.grd -Dp -Fg17 -Gtmp2.grd
#
# upsample to the original size
#
  gmt grdsample tmp2.grd -Gtmp3.grd
  gmt grd2xyz tmp3.grd | gmt surface -R$insar -T0.5 -Gcorrection.grd
#
# ----------------------
# PERFORM THE CORRECTION
# ----------------------
  echo "Correcting the interferogram..."
#
  gmt grdmath $insar correction.grd SUB = $output 
#
  if ( -f $output ) then
       echo "Interferogram corrected and stored as $output!"
  else
       echo "Something went wrong -- intf not corrected "
       exit 1
  endif
# ----------------------
#  clean up
# ----------------------
#
 rm tmp*
