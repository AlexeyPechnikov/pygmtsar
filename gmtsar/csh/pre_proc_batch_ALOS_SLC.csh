#!/bin/csh -f
#       $Id$
#
#  Xiaopeng Tong, Mar 2 2010 
#  modified by D. Sandwell MAR 11 2010
#
#  preprocess all the data based on data.in table file and generate: 
#  1. SLC files
#  2. PRM files 
#  3. time-baseline plot for user to create stacking pairs 

#  format in data.in table file: 
#  	line 1: master_name  
# 	line 2 and below: aligned_name

alias rm 'rm -f'
unset noclobber

#
# check the number of arguments 
# 
  if ($#argv != 2) then 
    echo ""
    echo "Usage: pre_proc_batch_ALOS_SLC.csh data.in batch.config"
    echo "       preprocess a set of images using a common rear_range and radius"
    echo ""
    echo "       format of data.in is:"
    echo "         line 1: master_name "
    echo "         line 2 and below: aligned_name"
    echo ""
    echo "       example of data.in for ALOS is:"
    echo "         IMG-HH-ALPSRP096010650-H1.0__A" 
    echo "         IMG-HH-ALPSRP089300650-H1.0__A"
    echo "         IMG-HH-ALPSRP236920650-H1.0__A"
    echo ""
    echo "Example: pre_proc_batch_ALOS_SLC.csh data.in batch.config"
    echo ""
    exit 1
  endif

#
# read parameters from configuration file
#

  set num_patches = `grep num_patches $2 | awk '{print $3}'`
  set earth_radius = `grep earth_radius $2 | awk '{print $3}'`
  set SLC_factor = `grep SLC_factor $2 | awk '{print $3}'`
  
  set commandline = ""
  if (!($earth_radius == "")) then
    set commandline = "$commandline -radius $earth_radius"
  endif
  if (!($num_patches == "")) then
    set commandline = "$commandline -npatch $num_patches"
  endif
  if (!($SLC_factor == "")) then
    set commandline = "$commandline -SLC_factor $SLC_factor"
  endif

  echo $commandline

#
# open and read data.in table 
#
  echo ""
  echo "START PREPROCESS A STACK OF IMAGES"
  echo ""
  echo "preprocess master image"
  set line1 = `awk 'NR==1 {print $0}' $1`
  set master = `echo $line1[1] | awk '{ print substr($1,8,length($1)-7)}'`
#
# unpack the master if necessary
#
  if(! -f IMG-HH-$master.SLC || ! -f IMG-HH-$master.PRM ) then
    echo "ALOS_pre_process_SLC IMG-HH-$master LED-$master -ALOS1 $commandline"
    ALOS_pre_process_SLC IMG-HH-$master LED-$master -ALOS1 $commandline
  endif

  set RAD = `grep earth_radius IMG-HH-$master.PRM | awk '{print $3}'`
  set SLC_factor = `grep SLC_factor IMG-HH-$master.PRM | awk '{print $3}'`
  set rng_samp_rate_m = `grep rng_samp_rate IMG-HH-$master.PRM | awk 'NR == 1 {printf("%d", $3)}'`
  set npatch = `grep num_patch IMG-HH-$master.PRM | awk '{print $3}'`

  baseline_table.csh IMG-HH-$master.PRM IMG-HH-$master.PRM >! baseline_table.dat
  baseline_table.csh IMG-HH-$master.PRM IMG-HH-$master.PRM GMT >! table.gmt
#
# loop and unpack the aligned image using the same earth radius and near range as the master image
#
  foreach line2 (`awk 'NR>1 {print $0}' $1`)
    echo "pre_proc_batch.csh"
    echo "preprocess aligned images"
     
    set aligned = ` echo $line2 | awk '{ print substr($1,8,length($1)-7)}'`
    if(! -f IMG-HH-$aligned.SLC || ! -f IMG-HH-$aligned.PRM ) then
      echo "ALOS_pre_process_SLC IMG-HH-$aligned LED-$aligned -ALOS1 -radius $RAD "
      ALOS_pre_process_SLC IMG-HH-$aligned LED-$aligned -ALOS1 -radius $RAD 
    endif
#
# check the range sampling rate of the aligned images and do conversion if necessary
#
    set rng_samp_rate_s = `grep rng_samp_rate IMG-HH-$aligned.PRM | awk 'NR == 1 {printf("%d", $3)}'`
    set t = `echo $rng_samp_rate_m $rng_samp_rate_s | awk '{printf("%1.1f\n", $1/$2)}'`
    if ($t == 1.0) then
      echo "The range sampling rate for master and aligned images are: "$rng_samp_rate_m
#
      baseline_table.csh IMG-HH-$master.PRM IMG-HH-$aligned.PRM >> baseline_table.dat
      baseline_table.csh IMG-HH-$master.PRM IMG-HH-$aligned.PRM GMT >> table.gmt
    else if ($t == 2.0) then
      echo "Convert the aligned image from FBD to FBS mode"
      ALOS_fbd2fbs_SLC IMG-HH-$aligned.PRM IMG-HH-$aligned"_"FBS.PRM
#
      baseline_table.csh IMG-HH-$master.PRM IMG-HH-$aligned"_"FBS.PRM >> baseline_table.dat
      baseline_table.csh IMG-HH-$master.PRM IMG-HH-$aligned"_"FBS.PRM GMT >> table.gmt
      echo "Overwriting the old aligned image"
      mv IMG-HH-$aligned"_"FBS.PRM IMG-HH-$aligned.PRM
      update_PRM IMG-HH-$aligned.PRM input_file IMG-HH-$aligned.SLC
      mv IMG-HH-$aligned"_"FBS.SLC IMG-HH-$aligned.SLC

    else if  ($t == 0.5) then
      echo "Use FBS mode image as master"
      exit 1
    else
      echo "The range sampling rate for master and aligned images are not convertable"
      exit 1
    endif

# end of the loop over aligned images
  end

#
# make baseline plots
#

  awk '{print 2006.5+($1-181)/365.25,$2,$7}' < table.gmt > text
#    awk '{print 2006.5+($1-181)/365.25,$2,9,$4,$5,$6,$7}' < table.gmt > text
  set region = `gmt gmtinfo text -C | awk '{print $1-0.5, $2+0.5, $3-500, $4+500}'`
# set region = `minmax text -C | awk '{print $1-0.5, $2+0.5, -1200, 1200}'`
  gmt pstext text -JX8.8i/6.8i -R$region[1]/$region[2]/$region[3]/$region[4] -D0.2/0.2 -X1.5i -Y1i -K -N -F+f8,Helvetica+j5 > stacktable_all.ps
  awk '{print $1,$2}' < text > text2
  gmt psxy text2 -Sp0.2c -G0 -R -JX -Ba1:"year":/a200g00f100:"baseline (m)":WSen -O >> stacktable_all.ps


  echo ""
  echo "END PREPROCESS A STACK OF IMAGES"
  echo ""

#
# clean up the mess
#
  rm text text2 
