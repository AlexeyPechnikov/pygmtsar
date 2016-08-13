#!/bin/csh -f
#       $Id$
#
#  D. Sandwell FEB 10 2010
#
unset noclobber
#
# script to convert a grd file to a kml file for Google Earth
#
if ($#argv < 2 || $#argv > 3) then
 echo " "
 echo "Usage: grd2kml.csh grd_file_stem cptfile [-R<west>/<east>/<south>/<north>] "
 echo " "
 echo "Example: grd2kml.csh phase phase.cpt "
 echo " "
 exit 1
endif 
#
#
set DX = `gmt grdinfo $1.grd -C | cut -f8`
set DPI = `gmt gmtmath -Q $DX INV RINT = `
#echo $DPI
gmt set COLOR_MODEL = hsv
gmt set PS_MEDIA = letter
#
gmt grdgradient $1.grd -Ggrad.grd -V -Nt0.7 -A60 
if ($#argv == 3) then
  gmt grdimage $1.grd -Igrad.grd -C$2 $3 -Jx1id -P -Y2 -X2 -Q  -V > $1.ps 
else if ($#argv == 2) then
  gmt grdimage $1.grd -Igrad.grd -C$2 -Jx1id -P -Y2 -X2 -Q  -V > $1.ps
endif
#
#   now make the kml and png
#
gmt ps2raster $1.ps -W+k+t"$1" -E$DPI -TG -P -S -V -F$1.png
rm $1.ps
rm grad.grd
rm ps2raster*
rm psconvert*
#
