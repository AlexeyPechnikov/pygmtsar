#!/bin/csh -f
#
#  Eric Lindsey, July 2022
#
# Script to create DEM for GMTSAR, relative to WGS84 ellipsoid
#
  if ($#argv != 4) then
    echo ""
    echo "Usage: make_dem.csh W E S N"
    echo "      Uses GMT server to download SRTM 1-arcsec data (@earth_relief_01s)"
    echo "      and removes the EGM96 geoid to make heights relative to WGS84."
    echo ""
    echo "Example: make_dem.csh -115 -112 32 35"
    echo ""
    exit 1
  endif
#
  echo ""
  echo "START: make_dem.csh"
  echo ""
#
# get region in GMT format
#
  set R = "-R$1/$2/$3/$4"
#
# need to set this for the distribution
#
  set sharedir = `gmtsar_sharedir.csh`
#
# get srtm data
#
  gmt grdcut @earth_relief_01s $R -Gdem_ortho.grd -V
#
# resample and remove geoid
#
  gmt grdsample $sharedir/geoid_egm96_icgem.grd -Rdem_ortho.grd -Ggeoid_resamp.grd -V
  gmt grdmath -V dem_ortho.grd geoid_resamp.grd ADD = dem.grd
#
# clean up 
#
  rm geoid_resamp.grd
#
  echo ""
  echo "created dem.grd, heights relative to WGS84 ellipsoid"
  echo ""
  echo "END: make_dem.csh"
  echo ""

