#!/bin/csh -f
#       $Id$
#
#
#    Xiaohua(Eric) XU, July 7, 2016
#
# Script for merging 3 subswaths TOPS interferograms and then unwrap and geocode for a stack of interferograms. 
#
#

  if ($#argv != 2) then
    echo ""
    echo "Usage: merge_batch.csh inputfile config_file"
    echo ""
    echo "Note: Inputfiles should be as following:"
    echo ""
    echo "      IF1_Swath1_Path:master.PRM:repeat.PRM,IF1_Swath2_Path:master.PRM:repeat.PRM,IF1_Swath3_Path:master.PRM:repeat.PRM"
    echo "      IF2_Swath1_Path:master.PRM:repeat.PRM,IF2_Swath2_Path:master.PRM:repeat.PRM,IF1_Swath3_Path:master.PRM:repeat.PRM"
    echo "      (Use the repeat PRM which contains the shift information.)"
    echo "      e.g. ../F1/intf_all/2015092_2015128/:S1A20150403_ALL_F1.PRM:S1A20150509_ALL_F1.PRM,../F2/intf_all/2015092_2015128/:S1A20150403_ALL_F2.PRM:S1A20150509_ALL_F2.PRM,../F3/intf_all/2015092_2015128/:S1A20150403_ALL_F3.PRM:S1A20150509_ALL_F3.PRM"
    echo ""
    echo "      Script for merging 3 subswaths TOPS interferograms and then unwrap and geocode for a stack of interferograms."
    echo ""
    echo "      Make sure under each path, the processed phasefilt.grd, corr.grd and mask.grd exist."
    echo "      Also make sure the dem.grd is linked. "
    echo "      If trans.dat exits, recomputation of projection matrix will not proceed."
    echo "      The master image of firet line should be the super_master."
    echo ""
    echo "      config_file is the same one used for processing."
    echo ""
    echo "Example: merge_batch.csh filelist batch.config"
    echo ""
    exit 1
  endif

  if (! -f dem.grd) then
    echo "dem.grd is required ..."
    exit 1
  endif

  set master = `grep master_image $2 | awk '{print $3}'`
  set input_file = $1

  
  set now_dir = `pwd`

  
  foreach line (`awk '{print $0}' $input_file`)
    set dir_name = `echo $line | awk -F, '{print $1}' | awk -F: '{print $1}' | awk -F"/" '{print $(NF-1)}'`
    mkdir $dir_name
    cd $dir_name
    echo $line | awk -F, '{for (i=1;i<=NF;i++) print "../"$i}' > tmp.filelist
    
    if (-f ../trans.dat) ln -s ../trans.dat .
    if (-f ../raln.grd) ln -s ../raln.grd .
    if (-f ../ralt.grd) ln -s ../ralt.grd .
    ln -s ../dem.grd .
    ln -s ../$2 .
    
    merge_unwrap_geocode_tops.csh tmp.filelist $2

    if (! -f ../trans.dat) then
      mv trans.dat ../
      ln -s ../trans.dat .
    endif
    if (! -f ../raln.grd) then
      mv raln.grd ../
      ln -s ../raln.grd .
    endif
    if (! -f ../ralt.grd) then
      mv ralt.grd ../
      ln -s ../ralt.grd .
    endif

    cd ..

  end
  

