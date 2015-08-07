#!/bin/csh -f
#
# Original version: Eric Lindsey, May 2015

if ($#argv != 3) then
  echo "  usage: detrend_before_unwrap.csh master #  bounds"
  echo "  master is the stem name for the master image, eg. ALOS2040533050-150222-WBDR1.1__D"
  echo "  # is the frame number (1-5) to use for detrending"
  echo "  bounds (eg. 0/1000/0/1000) is the GMT bounds in range/az coordinates to use"
  exit 1
endif

# correlation threshold for the detrending
set corr_threshold = 0.2

set master = $1

# go to the selected frame directory
cd F$2/intf/*/

if ( -e phasefilt_nodetrend.grd ) then
  mv phasefilt_nodetrend.grd phasefilt.grd
endif

# unwrap a part of the selected frame for detrending
mkdir detrend_part
cd detrend_part
ln -s ../corr.grd .
ln -s ../mask.grd .
ln -s ../phasefilt.grd .

#
## convert to range(meters) / azimuth (seconds) coordinates before fitting the trend
#
#speed of light
set cee = "299792458.0"

# convert coordinates for each frame to absolute range and azimuth (time)
#compute far range in meters
#range_meters = near_range + i * speed_of_light/(2*rng_samp_rate)
set near_range = `grep near_range ../*$master*PRM | head -n1 |awk '{print $3}'`
set rng_samp_rate = `grep rng_samp_rate ../*$master*PRM | head -n1 |awk '{print $3}'`
set numrng = `grep num_rng_bins ../*$master*PRM | head -n1 |awk '{print $3}'`
set far_range = `echo $near_range $numrng $cee $rng_samp_rate |awk '{printf("%.14f", $1 + $2*$3/(2.0*$4))}'`
set drng = `echo $far_range $near_range $numrng |awk '{printf("%.14f", ($1-$2)/$3)}'`

#compute start time relative to the start time of the selected frame, in seconds
#azi_time(n,j) = (SC_clock_start(1) - SC_clock_start(n)) + j / PRF(n)
set numaz = `grep num_lines ../*$master*PRM | head -n1 |awk '{print $3}'`
set PRF = `grep PRF ../*$master*PRM | head -n1 |awk '{print $3}'`
set azi_start = "0.0"
set azi_end   = `echo $numaz $PRF |awk '{printf("%.14f", $1/$2 )}'`
set daz = `echo $azi_end $azi_start $numaz |awk '{printf("%.14f", ($1-$2)/$3)}'`

# unwrap the patch
echo "unwrap the selected region"
#snaphu_interp_lindsey.csh $corr_threshold 0 $3
snaphu.csh $corr_threshold 1 $3

set rng_actual_min = `gmt grdinfo -I- unwrap.grd | sed -e 's:R:/:' | awk -F/ '{print $2}'`
set rng_actual_max = `gmt grdinfo -I- unwrap.grd | sed -e 's:R:/:' | awk -F/ '{print $3}'`
set az_actual_min  = `gmt grdinfo -I- unwrap.grd | sed -e 's:R:/:' | awk -F/ '{print $4}'`
set az_actual_max  = `gmt grdinfo -I- unwrap.grd | sed -e 's:R:/:' | awk -F/ '{print $5}'`

set near_range_crop = `echo $near_range $rng_actual_min $drng |awk '{print $1 + $2*$3}'`
set far_range_crop  = `echo $near_range $rng_actual_max $drng |awk '{print $1 + $2*$3}'`
set azi_start_crop  = `echo $azi_start $az_actual_min $daz |awk '{print $1 + $2*$3}'`
set azi_end_crop    = `echo $azi_start $az_actual_max $daz |awk '{print $1 + $2*$3}'`

# edit the coordinates
echo "edit unwrap.grd coordinates to -R$near_range_crop/$far_range_crop/$azi_start_crop/$azi_end_crop"
gmt grdedit unwrap.grd -R$near_range_crop/$far_range_crop/$azi_start_crop/$azi_end_crop
gmt grd2xyz unwrap.grd -S > unwrap.dat

# compute the trend
echo "find the trend from a simple 3-parameter inversion"
set params = `fit_planar_trend.py unwrap.dat`
set mean   = `echo $params |cut -d\  -f 1`
set rslope = `echo $params |cut -d\  -f 2`
set aslope = `echo $params |cut -d\  -f 3`
echo "found trend parameters (mean, range, az) $mean, $rslope, $aslope"

set SC_clock_start0 = `grep SC_clock_start ../*$master*PRM | head -n1 |awk '{print $3}'`

#back to top directory
cd ../../../../

foreach n ( 1 2 3 4 5 )
  cd F$n/intf/*/
  echo "detrending F$n"
  if ( -e phasefilt_nodetrend.grd ) then
    #don't detrend something that was already detrended
    mv phasefilt_nodetrend.grd phasefilt.grd
  endif
  mv phasefilt.grd phasefilt_nodetrend.grd

  # convert coordinates for each frame to absolute range and azimuth (time)
  #compute far range in meters:  range_meters = near_range + i * (rng_samp_rate)/(2*C)
  set near_range = `grep near_range *$master*PRM | head -n1 |awk '{print $3}'`
  set rng_samp_rate = `grep rng_samp_rate *$master*PRM | head -n1 |awk '{print $3}'`
  set numrng = `grep num_rng_bins *$master*PRM | head -n1 |awk '{print $3}'`
  set far_range = `echo $near_range $numrng $cee $rng_samp_rate |awk '{printf("%.14f", $1 + $2*$3/(2.0*$4))}'`

  #compute start time relative to the start time of the selected frame, in seconds
  set SC_clock_start = `grep SC_clock_start *$master*PRM | head -n1 |awk '{print $3}'`
  set numaz = `grep num_lines *$master*PRM | head -n1 |awk '{print $3}'`
  set PRF = `grep PRF *$master*PRM | head -n1 |awk '{print $3}'`
  set azi_start = `echo $SC_clock_start $SC_clock_start0 |awk '{printf("%.14f", 86400.0*($1-$2))}'`
  set azi_end   = `echo $SC_clock_start $SC_clock_start0 $numaz $PRF |awk '{printf("%.14f", 86400.0*($1-$2) + $3/$4)}'`

  # edit the coordinates
  echo "edit coordinates to -R$near_range/$far_range/$azi_start/$azi_end"
  gmt grdedit phasefilt_nodetrend.grd -R$near_range/$far_range/$azi_start/$azi_end

  # subtract the fitted trend
  gmt grdmath -V phasefilt_nodetrend.grd X $rslope MUL SUB Y $aslope MUL SUB 2 PI MUL MOD PI SUB = phasefilt.grd

  # edit the coordinates back
  gmt grdedit phasefilt.grd -Rcorr.grd
  gmt grdedit phasefilt_nodetrend.grd -Rcorr.grd

  #geocode again if it was done before
  if ( -e phasefilt_mask_ll.grd ) then
    echo "redoing geocoding of phasefilt"
    gmt grdmath phasefilt.grd mask2.grd MUL = phasefilt_mask.grd
    proj_ra2ll.csh trans.dat phasefilt_mask.grd phasefilt_mask_ll.grd
    grd2kml.csh phasefilt_mask_ll phase.cpt
  endif

  cd ../../../
end

