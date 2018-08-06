#!/bin/csh -f
#       $Id$
# Xiaohua Xu and David Sandwell Dec 23 2015
#
#  script to align S1 TOPS mode data 
#
#  1) Make PRM and LED files for both master and slave.
#
#  2) Do geometric back geocoding to make the range and azimuth alignment grids 
#
#  3) Make PRM, LED and SLC files for both master and slave that are aligned
#     at the fractional pixel level. 
#
alias rm 'rm -f'
unset noclobber
#
if ($#argv < 5) then
 echo " "
 echo "Usage: align_tops.csh master_prefix master_orb_file slave_s1a_prefix slave_orb_file dem.grd" 
 echo " "
 echo "Be sure the tiff, xml, orbit and dem files are available in the local directory."
 echo " "
 echo "Example: align_tops.csh s1a-iw3-slc-vv-20150526t014937-20150526t015002-006086-007e23-003 S1A_OPER_AUX_POEORB_OPOD_20150615T155109_V20150525T225944_20150527T005944.EOF.txt s1a-iw3-slc-vv-20150607t014937-20150607t015003-006261-00832e-006 S1A_OPER_AUX_POEORB_OPOD_20150627T155155_V20150606T225944_20150608T005944.EOF.txt dem.grd "
 echo " "
 echo "Output: S1_20150526_F3.PRM S1_20150526_F3.LED S1_20150526_F3.SLC S1_20150607_F3.PRM S1_20150607_F3.LED S1_20150607_F3.SLC "
 echo " "
 exit 1
endif 
#  
#  make sure the files are available
#
if(! -f $1.xml) then
   echo "****** missing file: "$1
   exit
endif
if(! -f $2) then
   echo "****** missing file: "$2
   exit
endif
if(! -f $3.xml) then
   echo "****** missing file: "$3
   exit
endif
if(! -f $4) then
   echo "****** missing file: "$4
   exit
endif
if(! -f $5) then
   echo "****** missing file: "$5
   exit
endif
# 
#  set the full names and create an output prefix
#
set mtiff = ` echo $1.tiff `
set mxml = ` echo $1.xml `
set stiff = ` echo $3.tiff `
set sxml = ` echo $3.xml `
set mpre = ` echo $1 | awk '{ print "S1_"substr($1,16,8)"_"substr($1,25,6)"_F"substr($1,7,1)}'`
set spre = ` echo $3 | awk '{ print "S1_"substr($1,16,8)"_"substr($1,25,6)"_F"substr($1,7,1)}'`
echo $mpre
echo $spre
#
#  1) make PRM and LED files for both master and slave but not the SLC file
#
make_s1a_tops $mxml $mtiff $mpre 0 
make_s1a_tops $sxml $stiff $spre 0 
#
#  replace the LED with the precise orbit
#
ext_orb_s1a $mpre".PRM" $2 $mpre
ext_orb_s1a $spre".PRM" $4 $spre
#
#  calculate the earth radius and make the slave match the master
#
calc_dop_orb $mpre".PRM" tmp 0 0 
cat tmp >> $mpre".PRM"
set earth_radius = `grep earth_radius tmp | awk '{print $3}'`
calc_dop_orb $spre".PRM" tmp2 $earth_radius 0
cat tmp2 >> $spre".PRM"
rm tmp tmp2#
#  2) do a geometric back projection to determine the alignment parameters
#
#  Filter and downsample the topography to 12 seconds or about 360 m
#
gmt grdfilter $5 -D3 -Fg2 -I12s -Ni -Gflt.grd 
gmt grd2xyz --FORMAT_FLOAT_OUT=%lf flt.grd -s > topo.llt
#
# map the topography into the range and azimuth of the master and slave using polynomial refinement
# can do this in parallel
#
# first check whether there are any burst shift
#
set lontie = `SAT_baseline $mpre".PRM" $spre".PRM" | grep lon_tie_point | awk '{print $3}'`
set lattie = `SAT_baseline $mpre".PRM" $spre".PRM" | grep lat_tie_point | awk '{print $3}'`
set tmp_am = `echo $lontie $lattie 0 | SAT_llt2rat $mpre".PRM" 1 | awk '{print $2}'`
set tmp_as = `echo $lontie $lattie 0 | SAT_llt2rat $spre".PRM" 1 | awk '{print $2}'`
set tmp_da = `echo $tmp_am $tmp_as | awk '{printf("%d",$2-$1)}'`
#
# if ther is, modify the master PRM start_time to get a better r/a estimate
#
if ($tmp_da > -1000 && $tmp_da < 1000) then
  SAT_llt2rat $mpre".PRM" 1 < topo.llt > master.ratll &
  SAT_llt2rat $spre".PRM" 1 < topo.llt > slave.ratll &
  wait
else
  echo "Modifying master PRM by $tmp_da lines..."
  cp $mpre".PRM" tmp.PRM
  set prf = `grep PRF tmp.PRM | awk '{print $3}'`
  set ttmp = `grep clock_start tmp.PRM | grep -v SC_clock_start | awk '{print $3}' | awk '{printf ("%.12f",$1 - '$tmp_da'/'$prf'/86400.0)}'`
  update_PRM tmp.PRM clock_start $ttmp
  set ttmp = `grep clock_stop tmp.PRM | grep -v SC_clock_stop | awk '{print $3}' | awk '{printf ("%.12f",$1 - '$tmp_da'/'$prf'/86400.0)}'`
  update_PRM tmp.PRM clock_stop $ttmp
  set ttmp = `grep SC_clock_start tmp.PRM | awk '{print $3}' | awk '{printf ("%.12f",$1 - '$tmp_da'/'$prf'/86400.0)}'`
  update_PRM tmp.PRM SC_clock_start $ttmp
  set ttmp = `grep SC_clock_stop tmp.PRM | awk '{print $3}' | awk '{printf ("%.12f",$1 - '$tmp_da'/'$prf'/86400.0)}'`
  update_PRM tmp.PRM SC_clock_stop $ttmp
#
#  restore the modified lines 
#
  #SAT_llt2rat tmp.PRM 1 < topo.llt > tmp.ratll &
  SAT_llt2rat tmp.PRM 1 < topo.llt > master.ratll &
  SAT_llt2rat $spre".PRM" 1 < topo.llt > slave.ratll &
  wait
  #echo "Restoring $tmp_da lines to master ashifts..."
  #awk '{printf("%.6f %.6f %.6f %.6f %.6f\n",$1,$2-'$tmp_da',$3,$4,$5)}' tmp.ratll > master.ratll
endif
#
#  paste the files and compute the dr and da
#
#paste master.ratll slave.ratll | awk '{printf("%.6f %.6f %.6f %.6f %d\n", $6, $6-$1, $7, $7-$2, "100")}' > tmp.dat
paste master.ratll slave.ratll | awk '{printf("%.6f %.6f %.6f %.6f %d\n", $1, $6 - $1, $2, $7 - $2, "100")}' > tmp.dat
#
#  make sure the range and azimuth are within the bounds of the slave 
#
set rmax = `grep num_rng_bins $spre".PRM" | awk '{print $3}'`
set amax = `grep num_lines $spre".PRM" | awk '{print $3}'`
if ($tmp_da > -1000 && $tmp_da < 1000) then
  awk '{if($1 > 0 && $1 < '$rmax' && $3 > 0 && $3 < '$amax') print $0 }' < tmp.dat > offset.dat
else
  awk '{if($1 > 0 && $1 < '$rmax' && $3 > 0 && $3 < '$amax') print $0 }' < tmp.dat > offset.dat
  awk '{if($1 > 0 && $1 < '$rmax' && $3 > 0 && $3 < '$amax') printf("%.6f %.6f %.6f %.6f %d\n", $1, $2, $3 - '$tmp_da', $4 + '$tmp_da', "100") }' < tmp.dat > offset2.dat
endif
#
#  extract the range and azimuth data
#
awk '{ printf("%.6f %.6f %.6f \n",$1,$3,$2) }' < offset.dat > r.xyz
awk '{ printf("%.6f %.6f %.6f \n",$1,$3,$4) }' < offset.dat > a.xyz

#
#  fit a surface to the range and azimuth offsets
#
gmt blockmedian r.xyz -R0/$rmax/0/$amax -I16/8 -r -bo3d > rtmp.xyz
gmt blockmedian a.xyz -R0/$rmax/0/$amax -I16/8 -r -bo3d > atmp.xyz
gmt surface rtmp.xyz -bi3d -R0/$rmax/0/$amax -I16/8 -T0.3 -Grtmp.grd -N1000  -r &
gmt surface atmp.xyz -bi3d -R0/$rmax/0/$amax -I16/8 -T0.3 -Gatmp.grd -N1000  -r &
wait
gmt grdmath rtmp.grd FLIPUD = r.grd
gmt grdmath atmp.grd FLIPUD = a.grd
#
#  3) make PRM, LED and SLC files for both master and slave that are aligned
#     at the fractional pixel level but still need a integer alignment from 
#     resamp
#  
#  make the new PRM files and SLC
#
make_s1a_tops $mxml $mtiff $mpre 1 
make_s1a_tops $sxml $stiff $spre 1 r.grd a.grd
#
#  resamp the slave and set the aoffset to zero
#
cp $spre".PRM" $spre".PRM0"
if ($tmp_da > -1000 && $tmp_da < 1000) then
  update_PRM $spre".PRM" ashift 0
else
  update_PRM $spre".PRM" ashift $tmp_da
  echo "Restoring $tmp_da lines with resamp..."
endif
resamp $mpre".PRM" $spre".PRM" $spre".PRMresamp" $spre".SLCresamp" 1
mv $spre".SLCresamp" $spre".SLC"
mv $spre".PRMresamp" $spre".PRM"
#
if ($tmp_da > -1000 && $tmp_da < 1000) then
  fitoffset.csh 3 3 offset.dat >> $spre.PRM
else
  fitoffset.csh 3 3 offset2.dat >> $spre.PRM
endif
#
#   re-extract the lED files
#
ext_orb_s1a $mpre".PRM" $2 $mpre
ext_orb_s1a $spre".PRM" $4 $spre
#
#  calculate the earth radius and make the slave match the master
#
calc_dop_orb $mpre".PRM" tmp 0 0
cat tmp >> $mpre".PRM"
set earth_radius = `grep earth_radius tmp | awk '{print $3}'`
calc_dop_orb $spre".PRM" tmp2 $earth_radius 0
cat tmp2 >> $spre".PRM"
rm tmp tmp2
#
rm topo.llt master.ratll slave.ratll *tmp* flt.grd r.xyz a.xyz *.PRM0
