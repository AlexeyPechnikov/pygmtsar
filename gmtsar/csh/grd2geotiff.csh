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
if ( -f ~/.quiet ) then
	set V = ""
	set VS = ""
else
	set V = "-V"
	set VS = "-S -V"
endif

set DX = `gmt grdinfo $1.grd -C | cut -f8`
set DPI = `gmt math -Q $DX INV RINT = `
echo $DPI
gmt set COLOR_MODEL = hsv
gmt set PAPER_MEDIA = tabloid
#
gmt grdgradient $1.grd -Ggrad.grd $V -Nt0.7 -A60 
if ($#argv == 3) then
  gmt grdimage $1.grd -Igrad.grd -C$2 $3 -Jx1id -P -Y2i -X2i -Q $V > $1.ps 
else if ($#argv == 2) then
  gmt grdimage $1.grd -Igrad.grd -C$2 -Jx1id -P -Y2i -X2i -Q $V > $1.ps
endif
#
#   now make the geotiff 
#
echo "Make $1.tiff"
gmt psconvert $1.ps -W+g+t"$1" -E$DPI -P $VS
rm -f $1.ps grad.grd
#
