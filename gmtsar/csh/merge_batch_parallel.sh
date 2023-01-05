#!/bin/bash
#       $Id$
#
#
#    Xiaohua(Eric) XU, July 7, 2016
#    Steffan Davies, December 21, 2022
#
# Script for merging 3 subswaths TOPS interferograms and then unwrap and geocode for a stack of interferograms. 
#
#

if [[ "$#" -ne 2 && "$#" -ne 3 ]]
then
  echo ""
  echo "Usage: merge_batch_parallel.sh inputfile config_file [det_stitch]"
  echo ""
  echo "Note: Inputfiles should be as following:"
  echo ""
  echo "      IF1_Swath1_Path:master.PRM:repeat.PRM,IF1_Swath2_Path:master.PRM:repeat.PRM,IF1_Swath3_Path:master.PRM:repeat.PRM"
  echo "      IF2_Swath1_Path:master.PRM:repeat.PRM,IF2_Swath2_Path:master.PRM:repeat.PRM,IF1_Swath3_Path:master.PRM:repeat.PRM"
  echo "      (Use the repeat PRM which contains the shift information.)"
  echo "      e.g. ../F1/intf_all/2015092_2015128/:S1A20150403_ALL_F1.PRM:S1A20150509_ALL_F1.PRM,../F2/intf_all/2015092_2015128/:S1A20150403_ALL_F2.PRM:S1A20150509_ALL_F2.PRM,../F3/intf_all/2015092_2015128/:S1A20150403_ALL_F3.PRM:S1A20150509_ALL_F3.PRM"
  echo ""
  echo "      Make sure under each path, the processed phasefilt.grd, corr.grd and mask.grd exist."
  echo "      Also make sure the dem.grd is linked. "
  echo "      If trans.dat exits, recomputation of projection matrix will not proceed."
  echo "      The master image of firet line should be the super_master."
  echo ""
  echo "      config_file is the same one used for processing."
  echo ""
  echo "      set det_stitch to 1 if you want to calculate stitching position based on the NaN area in the grids"
  echo ""
  echo "Example: merge_batch.sh filelist batch.config"
  echo ""
  exit 1
fi

if [[ ! -f dem.grd ]]
  then
  echo "dem.grd is required ..."
  exit 1
fi

if [[ "$#" -eq 3 ]]
then
  export det_stitch=1
else
  export det_stitch=""
fi

export input_file=$1
export config_file=$2

awk 'NR==1{print $0}' $input_file | awk -F, '{for (i=1;i<=NF;i++) print "../"$i}' | awk -F: '{print $1$2}'> tmpm.filelist 

export now_dir=`pwd`

parallel_func() {
  line=$1
  dir_name=`echo $line | awk -F, '{print $1}' | awk -F: '{print $1}' | awk -F"/" '{print $(NF-1)}'`
  mkdir $dir_name
  cd $dir_name
  echo $line | awk -F, '{for (i=1;i<=NF;i++) print "../"$i}' > tmp.filelist
  paste ../tmpm.filelist tmp.filelist | awk '{print $1","$2}' > tmp
  rm tmp.filelist

  for f_name in `awk '{print $0}' < tmp`
  do
      mm=`echo $f_name | awk -F, '{print $1}'`
      pth=`echo $f_name | awk -F, '{print $2}' | awk -F: '{print $1}'`
      f1=`echo $f_name | awk -F, '{print $2}' | awk -F: '{print $2}'`
      f2=`echo $f_name | awk -F, '{print $2}' | awk -F: '{print $3}'`
      cp $mm ./supermaster.PRM
      rshift=`grep rshift $pth$f1 | tail -1 | awk '{print $3}'`
      update_PRM supermaster.PRM rshift $rshift
      fs1=`grep first_sample supermaster.PRM | awk '{print $3}'`
      fs2=`grep first_sample $pth$f1 | awk '{print $3}'`
      [[ "$fs2" > "$fs1" ]] && update_PRM supermaster.PRM first_sample $fs2
      cp supermaster.PRM $pth
      echo $pth":supermaster.PRM:"$f2 >> tmp.filelist
  done

  [[ -f ../trans.dat ]] && ln -s ../trans.dat .
  [[ -f ../raln.grd ]] && ln -s ../raln.grd .
  [[ -f ../ralt.grd ]] && ln -s ../ralt.grd .
  [[ -f ../landmask_ra.grd ]] && ln -s ../landmask_ra.grd .
  ln -s ../dem.grd .
  ln -s ../"$config_file" .
  rm tmp

  echo `pwd`
  merge_unwrap_geocode_tops.csh tmp.filelist "$config_file" "$det_stitch"

  if [[ ! -f ../trans.dat && -f trans.dat ]]
  then
    mv trans.dat ../
    ln -s ../trans.dat .
  fi

  if [[ ! -f ../landmask_ra.grd && -f landmask_ra.grd ]]
  then
    mv landmask_ra.grd  ../
    ln -s ../landmask_ra.grd .
  fi

  if [[ ! -f ../raln.grd && -f raln.grd ]]
  then
    mv raln.grd ../
    ln -s ../raln.grd .
  fi

  if [[ ! -f ../ralt.grd && -f raln.grd ]]
  then
    mv ralt.grd ../
    ln -s ../ralt.grd .
  fi

  cd $now_dir
}

export -f parallel_func

# In case there is no trans.dat,
# run only on a single product to create it
# before parallel merging

if [[ ! -f ./trans.dat ]]
then
  cat $input_file | head -n 1 | parallel parallel_func
fi

parallel parallel_func :::: $input_file
