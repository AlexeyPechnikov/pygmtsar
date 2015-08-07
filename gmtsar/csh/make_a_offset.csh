#!/bin/tcsh -f
#       $Id$
# Script to drive the xcorr to do the azimuthal pixel-tracking
# Originated by Matt Wei, April 19, 2010
# Rewrote by Kang Wang 
# Last update: Aug. 09, 2013
#
# Revised by Xiaohua Xu, Nov 15, 2013: adding some pre_proc, redefined some parameters, add blockmedian, 
#									   add shaded part, delete some redundant lines, easier to use
#


if ($#argv != 6) then
   echo ""
   echo "Usage: make_azi_offset.csh Master.PRM Slave.PRM nx ny xsearch ysearch"
   echo ""
   exit 1
endif

echo "make_azi_offset.csh"

set master = $1
set slave = $2
set nx = $3
set ny = $4
set xsearch = $5
set ysearch = $6

#
# pre_proc needed files

cd intf/*/

#

mkdir azi_offset
cd azi_offset
cp ../../../raw/$1 .
cp ../../../raw/$2 .

ln -s ../../../SLC/*.SLC .
ln -s ../../../topo/trans.dat
ln -s ../../../topo/dem.grd

set PRF = `grep PRF $master |awk -F"=" '{print $2}'`
set SC_vel = `grep SC_vel $master|awk -F"=" '{print $2}'`

echo "Ground velocity: "$SC_vel

set azi_size = `echo $SC_vel $PRF|awk '{printf "%10.3f",$1/$2}' `

update_PRM.csh $master rshift 0
update_PRM.csh $master sub_int_r 0
update_PRM.csh $master stretch_r 0.0
update_PRM.csh $master a_stretch_r 0.0
update_PRM.csh $master ashift 0
update_PRM.csh $master sub_int_a 0.0
update_PRM.csh $master stretch_a 0.0
update_PRM.csh $master a_stretch_a 0.0


update_PRM.csh $slave rshift 0
update_PRM.csh $slave sub_int_r 0
update_PRM.csh $slave stretch_r 0.0
update_PRM.csh $slave a_stretch_r 0.0
update_PRM.csh $slave ashift 0
update_PRM.csh $slave sub_int_a 0.0
update_PRM.csh $slave stretch_a 0.0
update_PRM.csh $slave a_stretch_a 0.0 

#
# make azimuth offset
#
#xcorr $master $slave -nx $nx -ny $ny -xsearch $xsearch -ysearch $ysearch -noshift
awk '{if ($4>-1 && $4<1 && $5>1) print $1,$3,$4}'  freq_xcorr.dat >azi.dat

set xmin = `gmt gmtinfo azi.dat -C |awk '{print $1}'`
set xmax = `gmt gmtinfo azi.dat -C |awk '{print $2}'`
set ymin = `gmt gmtinfo azi.dat -C |awk '{print $3}'`
set ymax = `gmt gmtinfo azi.dat -C |awk '{print $4}'`

set xinc = `echo $xmax $xmin $nx |awk '{printf "%d", ($1-$2)/($3-1)}'`
set yinc = `echo $ymax $ymin $ny |awk '{printf "%d", ($1-$2)/($3-1)}'`

set xinc = `echo $xinc | awk '{print $1*4}'`
set yinc = `echo $yinc | awk '{print $1*4}'`

echo "xinc:"$xinc"  yinc:"$yinc

gmt blockmedian azi.dat -R$xmin/$xmax/$ymin/$ymax  -I$xinc/$yinc > azi_b.dat

mv azi_b.dat azi.dat

gmt xyz2grd azi.dat -R$xmin/$xmax/$ymin/$ymax  -I$xinc/$yinc -Gaoff.grd 

gmt grdmath aoff.grd $azi_size MUL = azi_offset.grd

gmt grd2cpt azi_offset.grd -Cjet -E30 > azioff.cpt

gmt grdimage azi_offset.grd -JX5i -Cazioff.cpt -Q -P > azioff.ps

#
# project to lon/lat coordinates
#

proj_ra2ll_ascii.csh trans.dat azi.dat aoff.llo

set xmin2 = `gmt gmtinfo aoff.llo -C |awk '{print $1}'`
set xmax2 = `gmt gmtinfo aoff.llo -C |awk '{print $2}'`
set ymin2 = `gmt gmtinfo aoff.llo -C |awk '{print $3}'`
set ymax2 = `gmt gmtinfo aoff.llo -C |awk '{print $4}'`

set xinc2 = `echo $xmax2 $xmin2 $nx |awk '{printf "%12.5f", ($1-$2)/($3-1)}'`
set yinc2 = `echo $ymax2 $ymin2 $ny |awk '{printf "%12.5f", ($1-$2)/($3-1)}'`
set xinc2 = `echo $xinc2 | awk '{print $1*4}'`
set yinc2 = `echo $yinc2 | awk '{print $1*4}'`

echo "xinc2:"$xinc2"  yinc2:"$yinc2

gmt xyz2grd aoff.llo -R$xmin2/$xmax2/$ymin2/$ymax2  -I$xinc2/$yinc2 -r -fg -Gaoff_ll.grd

gmt grdmath aoff_ll.grd $azi_size MUL = azi_offset_ll.grd

gmt grdsample dem.grd -Gs_dem.grd -R$xmin2/$xmax2/$ymin2/$ymax2 -I$xinc2/$yinc2 -F

gmt grdgradient s_dem.grd -Gdem_grd.grd -A45 -Nt1
cat > grey_tmp.cpt <<EOF
  0 170 170 170  50 204 204 204
 50 204 204 204 100 238 238 238
B   0   0   0
F 255 255 255
N 255   0   0
EOF

#
# plot the azimuth offset
#

set r_topo = `gmt grdinfo dem.grd -T100`
gmt makecpt -Cgrey_tmp.cpt $r_topo -Z > topo.cpt

set x1 = `gmt grdinfo azi_offset_ll.grd -C |awk '{print $2 }'`
set x2 = `gmt grdinfo azi_offset_ll.grd -C |awk '{print $3 }'`

set y1 = `gmt grdinfo azi_offset_ll.grd -C |awk '{print $4 }'`
set y2 = `gmt grdinfo azi_offset_ll.grd -C |awk '{print $5 }'`

set xlim = `echo $x2 $x1|awk '{print $1-$2}'`
set ylim = `echo $y2 $y1|awk '{print $1-$2}'`

set length = 20 
set width = 16

set scl1 = `echo $ylim $length | awk '{print $2/$1 }' `
set scl2 = `echo $xlim $width | awk '{print $2/$1}'`
set scl = `echo $scl1 $scl2 | awk '{if ($1<$2) {print $1} else {print $2} }'`

set bounds = `gmt grdinfo -I- azi_offset_ll.grd`

gmt gmtdefaults -Ds >.gmtdefaults4
gmt set MAP_FRAME_TYPE plain

gmt psbasemap -B0.25 -Jm$scl"c"  $bounds -K -P  >azioff_ll.ps
gmt pscoast  -J -R -K -O -P -Dh -I1 -W0.5p -S >> azioff_ll.ps

gmt grdimage s_dem.grd -Idem_grd.grd -J -R -Ctopo.cpt -K -P -O   -Q -B >>azioff_ll.ps

gmt grdimage azi_offset_ll.grd -Idem_grd.grd -J -R -Cazioff.cpt -Q -P -K -O >>azioff_ll.ps

gmt psscale -D1.5c/3c/5c/0.35c -Cazioff.cpt -I -E -B1::/:m: -O   >>azioff_ll.ps 
gmt ps2raster azioff_ll.ps -P -Tg

rm aoff_ll.grd aoff.grd aoff.llo azi.dat dem_grd.grd grey_tmp.cpt ps2rast* raln* ralt* s_dem.grd temp.dat topo.cpt
