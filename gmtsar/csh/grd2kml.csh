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
if ( -f ~/.quiet ) then
	set V = ""
	set VS = ""
else
	set V = "-V"
	set VS = "-S -V"
endif

#
set DX = `gmt grdinfo $1.grd -C | cut -f8`
set DPI = `gmt gmtmath -Q $DX INV RINT = `
#echo $DPI
gmt set COLOR_MODEL = hsv
gmt set PS_MEDIA = tabloid
#
if ($#argv == 3) then
  gmt grdimage $1.grd -C$2 $3 -Jx1id -P -Y2i -X2i -Q $V > $1.ps 
else if ($#argv == 2) then
  gmt grdimage $1.grd -C$2 -Jx1id -P -Y2i -X2i -Q $V > $1.ps
endif
#
#   now make the kml and png
#
echo "Make $1.kml and $1.png"
gmt psconvert $1.ps -W+k+t"$1" -E$DPI -TG -P $VS -F$1.png
rm -f $1.ps grad.grd ps2raster* psconvert*
#
