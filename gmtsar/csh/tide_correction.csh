#!/bin/csh -f
# by Xiaohua Xu, Jun 2018
# make solid earth tide correction 
#

if ($#argv < 2) then 
  echo ""
  echo " Usage: tide_correction.csh PRM1 PRM2 dem.grd"
  echo ""
  echo " Output: tide.grd"
  echo ""
  echo ""
  exit 1
endif

set prm = $1
set prm2 = $2
set dem = $3
set inc1 = "0.02/0.02"
set inc2 = "200/50"
#
# get the information from the two PRM files
#
set rng = `grep num_rng_bins $prm | awk '{print $3}'`
set azi = `grep num_lines $prm | awk '{print $3}'`
set wave = `grep wavelength $prm | awk '{print $3}'`
set tt1 = `grep SC_clock_start $prm | awk '{print $3}'`
set tt2 = `grep SC_clock_start $prm2 | awk '{print $3}'`
 
echo ""
echo "Computing tidal correction for dates $tt1 $tt2 ..."
echo ""
#
# downsanple the dem and compute the three components of the solid earth tide as E N and U
#  
gmt grdsample $dem -I$inc1 -Gtmp_dem.grd
gmt grd2xyz tmp_dem.grd | awk '{print $1,$2}' > topo.ll
solid_tide $tt1 < topo.ll > topo1.llenu
solid_tide $tt2 < topo.ll > topo2.llenu
paste topo1.llenu topo2.llenu | awk '{printf("%.6f %.6f %.12f %.12f %.12f \n",$1,$2,$8-$3,$9-$4,$10-$5)}' > topo.llenu
#
# generate the look vector for the downsampled dem and apply to the tide difference vector to make 
# LOS tide difference in the LL coordinates
#
gmt grd2xyz tmp_dem.grd > topo.llt
SAT_look $prm < topo.llt > topo.lltn
paste topo.llenu topo.lltn | awk '{printf("%.12f\n",($3*$9) + ($4*$10) + ($5*$11))}' > topo.d
#
#  project the LOS tide into radar coordinates and convert to phase in the range direction
#
SAT_llt2rat $prm 1 < topo.llt > topo.ratll
paste topo.ratll topo.d | awk '{printf("%.6f %.6f %.12f\n", $1,$2,-$6*2.0*2.0*3.141592653/'$wave')}' > topo.rad
gmt blockmedian topo.rad -R0/$rng/0/$azi -I$inc2 -r > tmp.rad
gmt surface tmp.rad -R0/$rng/0/$azi -I$inc2 -T0.5 -r -Gtide.grd 
#
#  clean up
#
rm tmp_dem.grd tmp.rad topo*
