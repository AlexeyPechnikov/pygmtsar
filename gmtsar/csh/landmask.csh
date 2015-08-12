#!/bin/csh -f
#       $Id$

# make a landmask 

#
if ($#argv != 1) then
  echo ""
  echo "Usage: landmask.csh region_cut[0/10600/0/27648]"
  echo ""
  echo "    make a landmask in radar coordinates "
  echo "NOTE: The region_cut can be specified in batch.config file"
  echo ""
  exit 1
endif

echo ""
echo "MAKE LANDMASK -- START"
echo "REQUIRE FULL RESOLUTION COASTLINE FROM GMT"
echo ""
   
gmt grdlandmask -Glandmask.grd `gmt grdinfo -I- dem.grd` `gmt grdinfo -I dem.grd`  -V -NNaN/1 -Df
proj_ll2ra.csh trans.dat landmask.grd landmask_ra.grd
# if the landmask region is smaller than the region_cut pad with NaN
gmt grd2xyz landmask_ra.grd -bo > landmask_ra.xyz
gmt xyz2grd landmask_ra.xyz -bi -r -R$1 `gmt grdinfo -I landmask_ra.grd` -Gtmp.grd 
mv tmp.grd landmask_ra.grd
gmt grdsample landmask_ra.grd -Gtmp.grd -R$1 -I4/8 -nl+t0.1
mv tmp.grd landmask_ra.grd

# cleanup
rm landmask.grd landmask_ra.xyz

echo "MAKE LANDMASK -- END"
#
