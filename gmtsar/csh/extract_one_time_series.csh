#!/bin/csh -f
#       $Id$
# Xiaohua Xu, Jul 16, 2021
# 
# Script used to extract time-series from a list of displacement grids
# based on the input longitude and latitude
#
# Require PRM, LED, dem.grd to figure out the range azimuth
#
if ($#argv != 5 && $#argv != 7) then
    echo " "
    echo "Usage: extract_one_time_series.csh longitude latitude PRM dem.grd scene.tab [m_rng m_azi]"
    echo " "
    echo "   m_rng and m_azi is the number of pixels to be averaged for output, default is 5 5"
    echo "   scene.tab is the same input needed by sbas, dem.grd is the elevation model data"
    echo "   output (time_series.dat) will have two colums of data, displacements and standard deviation"
    echo " "
    echo "Example: extract_one_time_series.csh -119.782 36.292 supermaster.PRM dem.grd scene.tab 5 5"
    echo " "
    exit 1
endif

set N = `wc -l $5 | awk '{print $1}'`

set rng = `echo $1 $2 | gmt grdtrack -G$4 | SAT_llt2rat $3 1 | awk '{print $1}'`
set azi = `echo $1 $2 | gmt grdtrack -G$4 | SAT_llt2rat $3 1 | awk '{print $2}'`
echo "extracting at position range $rng and azimuth $azi ..."

set output = time_series.dat
if (-f $output) rm $output

set grid = `awk 'NR==1{print $1}' $5 | awk '{print "disp_"$1".grd"}'` 

set xinc = `gmt grdinfo $grid -C | awk '{print $8}'`
set yinc = `gmt grdinfo $grid -C | awk '{print $9}'`
if ($#argv == 7) then
    set samp_x = `echo $6 | awk '{print int($1/2)}'`
    set samp_y = `echo $7 | awk '{print int($1/2)}'`
else
    set samp_x = 2
    set samp_y = 2
endif
echo "grid xinc: $xinc, yinc: $yinc; average half samp: $samp_x $samp_y ..."
#set xinc = `gmt grdinfo -C ../combine/vel_combined_gs.grd -I | awk '{print substr($0,3,100)}' | awk -F'/' '{print $1}'`
#set yinc = `gmt grdinfo -C ../combine/vel_combined_gs.grd -I | awk '{print substr($0,3,100)}' | awk -F'/' '{print $2}'`
#set RR = `echo $lon $xinc $lat $yinc | awk '{printf("%.6f/%.6f/%.6f/%.6f",$1-2*$2,$1+3*$2,$3-2*$4,$3+3*$4)}'`
set RR = `echo $rng $xinc $azi $yinc | awk '{printf("%d/%d/%d/%d",int($1/$2)*$2-3*$2,int($1/$2)*$2+2*$2,int($3/$4)*$4-3*$4,int($3/$4)*$4+2*$4)}'`
echo "sampling over $RR ..."
#gmt grdcut ../combine/vel_combined_gs.grd -Gtmp_cut.grd -R$RR -Vq
#gmt grdinfo tmp_cut.grd -L2 -C | awk '{print $12,$13}' >> $output

set ii = 1
#set xinc = `gmt grdinfo -C vel.grd -I | awk '{print substr($0,3,100)}' | awk -F'/' '{print $1}'`
#set yinc = `gmt grdinfo -C vel.grd -I | awk '{print substr($0,3,100)}' | awk -F'/' '{print $2}'`
#set RR = `echo $rng $xinc $azi $yinc | awk '{printf("%.6f/%.6f/%.6f/%.6f",$1-2*$2,$1+3*$2,$3-2*$4,$3+3*$4)}'`

#gmt grdcut vel.grd -Gtmp_cut.grd -R$RR -Vq
#gmt grdinfo tmp_cut.grd -L2 -C | awk '{print $12,$13}' >> $output

#echo "" >> $output

while ($ii <= $N)
  #set name = `echo $ii | awk '{printf("disp_%.3d.grd",$1)}'`
  set name = `awk 'NR=='$ii'{print $1}' $5 | awk '{printf("disp_%d.grd",$1)}'`
  echo "Working on date $name"
  gmt grdcut $name -Gtmp_cut.grd -R$RR
  gmt grdinfo tmp_cut.grd -L2 -C | awk '{print $12,$13}' >> $output
  set ii = `echo $ii | awk '{print $1+1}'`
end
rm tmp_cut.grd

