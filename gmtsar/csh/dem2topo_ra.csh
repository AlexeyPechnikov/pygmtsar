#!/bin/csh -f
#       $Id$
# Matt WEI Feb 1 2010
# modified by Xiaopeng Feb 9 2010
# modified by E. Fielding, DST, XT to add TSX data Jan 10 2014
#=======================================================================
#  script to make topography for interferograms 
#  The USGS elevations are height above WGS84 so this is OK.
#
alias rm 'rm -f'
unset noclobber
#
if ($#argv != 3 && $#argv != 2 ) then
 echo " "
 echo "Usage: dem2topo_ra.csh master.PRM dem.grd [interpolation_approach]" 
 echo " "
 echo "    Note: interpolation_approach default is 0:surface with tension"
 echo "          set to 1 for gmt triangulate"
 echo " "
 exit 1
endif 
#
# local variables 
#
  set scale = -JX7i
if ( -f ~/.quiet ) then
    set V = ""
else
	set V = "-V"
endif
#
# set tension
#
set tension = 0.1
#
# set interpolation approach
#
set mode = 0
if ($#argv == 3) then
  set mode = $3
endif
#
#========================Mosaic topo data===============================

#-----------------------------------------------------------------------
# 
#------------------------Get bounds in radar coordinates----------------
set XMAX = `grep num_rng_bins $1 | awk '{print $3}'`
set yvalid = `grep num_valid_az $1 | awk '{print $3}'`
set num_patch = `grep num_patches $1 | awk '{print $3}'`
set YMAX = `echo "$yvalid $num_patch" | awk '{print $1*$2}'`
set SC = `grep SC_identity $1 | awk '{print $3}'`
set PRF = `grep PRF *.PRM | awk 'NR == 1 {printf("%d", $3)}'`
set region = 0/$XMAX/0/$YMAX
echo "Working over $region ... "
#
# look for range sampling rate
#
  set rng_samp_rate = `grep rng_samp_rate $1 | awk 'NR == 1 {printf("%d", $3)}'`
#
# set the range spacing of simulation in units of image range pixel size
#
if($rng_samp_rate > 0 && $rng_samp_rate < 25000000) then
  set rng = 1
else if($rng_samp_rate >= 25000000 && $rng_samp_rate < 72000000 || $SC == 7 ) then
  set rng = 2
else if($rng_samp_rate >= 72000000) then
  set rng = 4
else
   echo "range sampling rate out of bounds"
   exit 0
endif
echo " range decimation is: " $rng

if($SC == 10) then
     gmt grd2xyz --FORMAT_FLOAT_OUT=%lf $2 -s | SAT_llt2rat $1 1 -bod  > trans.dat
  else
     gmt grd2xyz --FORMAT_FLOAT_OUT=%lf $2 -s | SAT_llt2rat $1 0 -bod  > trans.dat
endif
#
# use an azimuth spacing of 2 for low PRF data such as S1 TOPS
#
if ($PRF < 1000) then
  gmt gmtconvert trans.dat -o0,1,2 -bi5d -bo3d | gmt blockmedian -R$region -I$rng/2 -bi3d -bo3d -r $V > temp.rat 
  if ($mode == 0) then
    gmt surface temp.rat -R$region -I$rng/2 -bi3d -T$tension -N1000 -Gpixel.grd -r -Q >& tmp
    set RR = `grep Hint tmp | head -1 | awk '{for(i=1;i<=NF;i++) print $i}' | grep /`
    if ("x$RR" == "x") then
      gmt surface temp.rat -R$region -I$rng/2 -bi3d -T$tension -N1000 -Gpixel.grd -r $V
    else
      gmt surface temp.rat $RR -I$rng/2 -bi3d -T$tension -N1000 -Gpixel.grd -r $V
      gmt grdcut pixel.grd -R$region -Gtmp.grd
      mv tmp.grd pixel.grd
    endif
  else if ($mode == 1) then
    set rng2 = `echo $rng | awk '{print $1*8}'`
    gmt triangulate temp.rat -R$region -I$rng/2 -bi3d -Gtopo_ra_tmp.grd -r $V
    gmt blockmean temp.rat -R$region -I$rng2/16 -bi3d -bo3d -r > mean.rat
    gmt surface mean.rat -R$region -I$rng2/16 -bi3d -T0.1 -N1000 -Gcoarse.grd -r
    gmt grdfill topo_ra_tmp.grd -Agcoarse.grd -Gpixel.grd
    rm topo_ra_tmp.grd coarse.grd mean.rat
  endif
else
  gmt gmtconvert trans.dat -o0,1,2 -bi5d -bo3d | gmt blockmedian -R$region -I$rng/4 -bi3d -bo3d -r $V > temp.rat 
  if ($mode == 0) then
    gmt surface temp.rat -R$region -I$rng/4 -bi3d -T$tension -N1000 -Gpixel.grd -r -Q >& tmp
    set RR = `grep Hint tmp | head -1 | awk '{for(i=1;i<=NF;i++) print $i}' | grep /`
    if ("x$RR" == "x") then
      gmt surface temp.rat -R$region -I$rng/4 -bi3d -T$tension -N1000 -Gpixel.grd -r $V
    else
      gmt surface temp.rat $RR -I$rng/4 -bi3d -T$tension -N1000 -Gpixel.grd -r $V
      gmt grdcut pixel.grd -R$region -Gtmp.grd
      mv tmp.grd pixel.grd
    endif  
  else if ($mode == 1) then
    set rng2 = `echo $rng | awk '{print $1*8}'`
    gmt triangulate temp.rat -R$region -I$rng/4 -bi3d -Gtopo_ra_tmp.grd -r $V
    gmt blockmean temp.rat -R$region -I$rng2/32 -bi3d -bo3d -r > mean.rat
    gmt surface mean.rat -R$region -I$rng2/32 -bi3d -T0.1 -N1000 -Gcoarse.grd -r
    gmt grdfill topo_ra_tmp.grd -Agcoarse.grd -Gpixel.grd
    rm topo_ra_tmp.grd coarse.grd mean.rat
  endif
endif
# 
# flip top to bottom for both ascending and descending passes
#  
  gmt grdmath pixel.grd FLIPUD = topo_ra.grd

#if ($#argv == 3) then
#  set x0 = echo $region | awk -F'/' '{print $1}'
#  set y0 = echo $region | awk -F'/' '{print $3}'
#  set x1 = echo $region | awk -F'/' '{printf("%d",$1 - '$x0')}'
#  set y1 = echo $region | awk -F'/' '{printf("%d",$3 - '$y0')}'
#  gmt grdedit topo_ra.grd -R0/$x1/0/$y1 
#endif

# 
# plotting
# 
  gmt grd2cpt topo_ra.grd -Cgray $V -Z > topo_ra.cpt 
  gmt grdimage topo_ra.grd $scale -P -Ctopo_ra.cpt -Bxaf+lRange -Byaf+lAzimuth -BWSen $V -K > topo_ra.ps
  gmt psscale -Rtopo_ra.grd -J -DJTC+w5i/0.2i+h -Ctopo_ra.cpt -Bxaf -By+lm -O >> topo_ra.ps
  gmt psconvert -Tf -P -A -Z topo_ra.ps
  echo "Topo range/azimuth map: topo_ra.pdf"
#
#  clean up
#
rm pixel.grd temp.rat dem.xyz tmp
rm topo_ra.cpt
