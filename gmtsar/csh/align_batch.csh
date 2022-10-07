#!/bin/tcsh -f
#       $Id$

# Align RAW or SLC images using orbits and DEM
# by Kang Wang in June. 2018
#
if ($#argv != 3) then
   echo ""
   echo "Usage: align_batch_geometry.csh data_type(RAW/SLC) secondary_align(0/1) data.in"
   echo ""
   echo "       align a set of images to a common master using geometric alignment "
   echo "       used the precise orbit and a dem"
   echo "       run this at the top level directory"
   echo "       it will create a new SLC directory "
   echo ""
   echo "  format of data.in is the same as data.in in used in pre_proc_batch.csh:"
   echo "       format of data.in is:"
   echo "         line 1: master_name "
   echo "         line 2 and below: aligned_name"
   echo ""
   echo "       example of data.in for ALOS is:"
   echo "         IMG-HH-ALPSRP096010650-H1.0__A"
   echo "         IMG-HH-ALPSRP089300650-H1.0__A"
   echo "         IMG-HH-ALPSRP236920650-H1.0__A"
   echo ""
   echo "       example of data.in for ERS is:"
   echo "         e1_05783"
   echo "         e1_07787"
   echo "         e1_10292"
   echo ""
   echo "       example of data.in for ENVISAT is:"
   echo "         ENV1_2_127_2925_07195"
   echo "         ENV1_2_127_2925_12706"
   echo "         ENV1_2_127_2925_13207"
   echo ""
   echo "       example of data.in for ENVI_SLC is:"
   echo "         ASA_IMS_1PNESA20030908_175832_000000182019_00399_07968_0000"
   echo "         ASA_IMS_1PNESA20040719_175832_000000182028_00399_12477_0000"
   echo "         ASA_IMS_1PNESA20051121_175837_000000172042_00399_19491_0000"
   echo ""
   echo "Example: align_batch_geometry.csh RAW 1 data.in"
   echo ""
   exit 1
endif

if ($1 != RAW && $1 != SLC) then
  echo "data type must be RAW or SLC"
  exit 1
endif
#
# make working directories
#
  mkdir -p SLC/
#
# clean up
#
  cleanup.csh SLC
  echo ""
  echo "START ALIGN A STACK OF IMAGES"
  echo ""
#
# link the needed files from the raw directory
#
  cd SLC
  ln -s ../raw/*.PRM .
  ln -s ../raw/*.LED .
  if ($1 == RAW) then
    ln -s ../raw/*.raw .
  else
    ln -s ../raw/*.SLC .
  endif
#
# make sure dem.grd is there
#
 if (! -e ../topo/dem.grd) then
   echo "DEM file not found!"
   exit 1
 endif

 echo "Downsample the DEM data"
 gmt grdfilter ../topo/dem.grd -D3 -Fg2 -I12s -Ni -Gflt.grd
 gmt grd2xyz --FORMAT_FLOAT_OUT=%lf flt.grd -s > topo.llt

 set master = `awk '{if (NR==1) print $1}' ../$3`
 awk '{if (NR>1) print $1}' ../$3 > aligned.list

 mv $master.PRM $master.PRM0

 echo "calculating the SAT height for master $master"
 calc_dop_orb $master.PRM0 $master.log 0
 cat $master.PRM0 $master.log > $master.PRM
 echo "fdd1                    = 0" >> $master.PRM
 echo "fddd1                   = 0" >> $master.PRM
 rm -f $master.log
 
 set earth_radius = `grep earth_radius $master.PRM|tail -1 | awk '{print $3}'`
 set fd1 = `grep fd1 $master.PRM |tail -1 |awk '{print $3}'`
 update_PRM $master.PRM earth_radius $earth_radius
 update_PRM $master.PRM fd1 $fd1

 SAT_llt2rat $master".PRM" 1 < topo.llt > master.ratll

  set rmax = `grep num_rng_bins $master".PRM" | awk '{print $3}'`
  set amax = `grep num_lines $master".PRM" | awk '{print $3}'`

if ($1 == RAW) then
 echo "Focusing the master - START"
  sarp.csh $master.PRM
 echo "Focusing the master - END"
endif

set naligned = `cat aligned.list |wc -l`
set ialigned = 0
foreach aligned (`cat aligned.list`)
    @ ialigned = $ialigned + 1
    echo "working on " $aligned  " [ $ialigned / $naligned ]"


    mv $aligned.PRM $aligned.PRM0
    calc_dop_orb $aligned.PRM0 $aligned.log $earth_radius

    cat $aligned.PRM0 $aligned.log > $aligned.PRM
    echo "fdd1                    = 0" >> $aligned.PRM
    echo "fddd1                   = 0" >> $aligned.PRM
    
    update_PRM $aligned.PRM earth_radius $earth_radius
    update_PRM $aligned.PRM fd1 $fd1

    rm -f $aligned.log

    if ($1 == RAW) then
     echo "Focusing the aligned - START"
     sarp.csh $aligned.PRM
     echo "Focusing the aligned - END"
    endif

    SAT_llt2rat $aligned".PRM" 1 < topo.llt > aligned.ratll
    paste master.ratll aligned.ratll | awk '{printf("%.6f %.6f %.6f %.6f %d\n", $1, $6 - $1, $2, $7 - $2, "100")}' > tmp.dat
    awk '{if($1 > 0 && $1 < '$rmax' && $3 > 0 && $3 < '$amax') print $0 }' < tmp.dat > offset.dat
    awk '{ printf("%.6f %.6f %.6f \n",$1,$3,$2) }' < offset.dat > r.xyz
    awk '{ printf("%.6f %.6f %.6f \n",$1,$3,$4) }' < offset.dat > a.xyz

   fitoffset.csh 3 3 offset.dat >> $aligned".PRM"
   resamp $master.PRM $aligned.PRM $aligned.PRMresamp $aligned.SLCresamp 4
   mv $aligned.SLC $aligned.SLC_old
   mv $aligned.SLCresamp $aligned.SLC
   cp $aligned.PRMresamp $aligned.PRM
#
#  do a secondary alignment with xcorr if requested
#
   if ($2 > 0) then
   cp $aligned.PRM tmp.PRM
   update_PRM tmp.PRM rshift 0
   update_PRM tmp.PRM ashift 0
   update_PRM tmp.PRM sub_int_r 0.0
   update_PRM tmp.PRM sub_int_a 0.0
   update_PRM tmp.PRM stretch_r 0.0
   update_PRM tmp.PRM stretch_a 0.0
   update_PRM tmp.PRM a_stretch_r 0.0
   update_PRM tmp.PRM a_stretch_a 0.0
   xcorr $master.PRM tmp.PRM -xsearch 128 -ysearch 256 -nx 15 -ny 30
   fitoffset.csh 1 1 freq_xcorr.dat 20 >> tmp.PRM
   resamp $master.PRM tmp.PRM $aligned.PRMresamp $aligned.SLCresamp 4
   mv $aligned.SLC $aligned.SLC_old
   mv $aligned.SLCresamp $aligned.SLC
# 
#  keep the old $aligned.PRM  because it has nearly the correct shifts needed for phasediff
#
#   cp $aligned.PRMresamp $aligned.PRM
  endif
end
#
#  cleanup the mess
#  
  rm tmp*
  rm *old
  rm *resamp
  rm *.raw *.LED *.PRM0
cd ..
