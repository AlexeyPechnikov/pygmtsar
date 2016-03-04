#!/bin/csh -f
#
#  D. Sandwell FEB 18 2015
#
unset noclobber
#
# script to convert a grd file to a geotiff file
#
if ($#argv < 2 || $#argv > 3) then
 echo " "
 echo "Usage: grd2geotiff.csh grd_file_stem cptfile [-R<west>/<east>/<south>/<north>] "
 echo " "
 echo "Example: grd2geotiff.csh phase phase.cpt "
 echo " "
 exit 1
endif 
#
#
set DX = `gmt grdinfo $1.grd -C | cut -f8`
set DPI = `gmtmath -Q $DX INV RINT = `
echo $DPI
gmtset COLOR_MODEL = hsv
gmtset PAPER_MEDIA = letter
#
grdgradient $1.grd -Ggrad.grd -V -Nt0.7 -A60 
if ($#argv == 3) then
  grdimage $1.grd -Igrad.grd -C$2 $3 -Jx1id -P -Y2 -X2 -Q  -V > $1.ps 
else if ($#argv == 2) then
  grdimage $1.grd -Igrad.grd -C$2 -Jx1id -P -Y2 -X2 -Q  -V > $1.ps
endif
#
#   now make the geotiff 
#
ps2raster $1.ps -W+g+t"$1" -E$DPI -P -S -V
rm $1.ps
rm grad.grd
#
