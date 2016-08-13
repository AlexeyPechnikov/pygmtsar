#!/bin/csh -f
#       $Id$
# Xiaohua Xu and David Sandwell Dec 23 2015
#
#  script to align S1A TOPS mode data 
#
#  1) Make PRM and LED files for both master and slave.
#
#  2) Do geometric back projection to determine the alignment parameters.
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
 echo "Output: S1A20150526_F3.PRM S1A20150526_F3.LED S1A20150526_F3.SLC S1A20150607_F3.PRM S1A20150607_F3.LED S1A20150607_F3.SLC "
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
set mpre = ` echo $1 | awk '{ print "S1A"substr($1,16,8)"_"substr($1,25,6)"_F"substr($1,7,1)}'`
set spre = ` echo $3 | awk '{ print "S1A"substr($1,16,8)"_"substr($1,25,6)"_F"substr($1,7,1)}'`
set swath = ` echo $1 | awk '{ print substr($1,7,1)}'`
echo $mpre
echo $spre
#
#  1) make PRM and LED files for both master and slave but not the SLC file
#
make_s1a_tops_6par $mxml $mtiff $mpre 0 0. 0. 0. 0. 0. 0.
make_s1a_tops_6par $sxml $stiff $spre 0 0. 0. 0. 0. 0. 0.
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
rm tmp tmp2
#
#  2) do a geometric back projection to determine the alignment parameters
#
#  Downsample the topography. The topo value is irrelevant so it can be filtered and downsampled.
#
gmt grdfilter $5 -D2 -Fg4 -I30s -Gflt.grd 
gmt grd2xyz --FORMAT_FLOAT_OUT=%lf flt.grd -s > topo.llt
#
# zero the alignment parameters of the slave image
# don'd change the alignment parameters of the master image in case this is a surrogate master
#
update_PRM.csh $spre".PRM" rshift 0
update_PRM.csh $spre".PRM" sub_int_r 0.0
update_PRM.csh $spre".PRM" stretch_r 0.0
update_PRM.csh $spre".PRM" a_stretch_r 0.0
update_PRM.csh $spre".PRM" ashift 0
update_PRM.csh $spre".PRM" sub_int_a 0.0
update_PRM.csh $spre".PRM" stretch_a 0.0
update_PRM.csh $spre".PRM" a_stretch_a 0.0
#
# map the topography into the range and azimuth of the master and slave using polynomial refinement
#
SAT_llt2rat $mpre".PRM" 1 < topo.llt > master.ratll
SAT_llt2rat $spre".PRM" 1 < topo.llt > slave.ratll
#
#  paste the files and compute the dr and da
#
#paste master.ratll slave.ratll | awk '{print( $6, $6-$1, $7, $7-$2, "100")}' > tmp.dat
paste master.ratll slave.ratll | awk '{printf("%.6f %.6f %.6f %.6f %d\n", $6, $6-$1, $7, $7-$2, "100")}' > tmp.dat
#
#  make sure the range and azimuth are within the bounds of the master
#
set rmax = `grep num_rng_bins $spre".PRM" | awk '{print $3}'`
set amax = `grep num_lines $spre".PRM" | awk '{print $3}'`
awk '{if($1 > 0 && $1 < '$rmax' && $3 > 0 && $3 < '$amax') print $0 }' < tmp.dat > offset.dat
#
#  run fitoffset
#
fitoffset.csh 3 3 offset.dat >> $spre".PRM"
#
#  save the rshift and ashift for the end
#
set rshift = `grep rshift $spre".PRM" | tail -1 | awk '{print $3}'`
set ashift = `grep ashift $spre".PRM" | tail -1 | awk '{print $3}'`
#
# clean up the mess
#
rm topo.llt master.ratll slave.ratll tmp.dat offset.dat flt.grd
#
#  3) make PRM, LED and SLC files for both master and slave that are aligned
#     at the fractional pixel level but still need a integer alignment from 
#     resamp
#
set sub_int_r = `grep sub_int_r $spre".PRM" | tail -1 | awk '{print $3}'`
set sub_int_a = `grep sub_int_a $spre".PRM" | tail -1 | awk '{print $3}'`
set stretch_r = `grep stretch_r $spre".PRM" | grep -v a_stretch_r | tail -1 | awk '{print $3}'`
set stretch_a = `grep stretch_a $spre".PRM" | grep -v a_stretch_a | tail -1 | awk '{print $3}'`
set a_stretch_a = `grep a_stretch_a $spre".PRM" | tail -1 | awk '{print $3}'`
set a_stretch_r = `grep a_stretch_r $spre".PRM" | tail -1 | awk '{print $3}'`
#echo $rshift $sub_int_r $stretch_r $a_stretch_r $ashift $sub_int_a $stretch_a $a_stretch_a
#  
#  make the new PRM files and SLC
#
make_s1a_tops_6par $mxml $mtiff $mpre 1 0. 0. 0. 0. 0. 0.
make_s1a_tops_6par $sxml $stiff $spre 1 $sub_int_r $sub_int_a $stretch_a $a_stretch_a $stretch_r $a_stretch_r
#
#   restore the integer rshift and ashift
#
update_PRM.csh $spre".PRM" rshift $rshift
update_PRM.csh $spre".PRM" ashift $ashift
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
# do integer resampling of the slave
#
resamp $mpre".PRM" $spre".PRM" $spre".PRMresamp" $spre".SLCresamp" 1
mv $spre".SLCresamp" $spre".SLC"
mv $spre".PRMresamp" $spre".PRM"
#
rm topo.llt master.ratll slave.ratll *tmp* flt.grd r.xyz a.xyz
