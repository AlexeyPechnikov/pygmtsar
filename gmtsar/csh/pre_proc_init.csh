#!/bin/csh -f
#       $Id$
#
#  Xiaopeng Tong, Mar 2 2010 
#  modified by D. Sandwell MAR 11 2010
#
#  preprocess all the data based on data.in table file and generate: 
#  1. raw files
#  2. PRM files 
#  3. time-baseline plot for user to create stacking pairs 

#  format in data.in table file: 
#  	line 1: master_name 
# 	line 2 and below: slave_name

alias rm 'rm -f'
unset noclobber

#
# check the number of arguments 
# 
  if ($#argv != 2) then 
    echo ""
    echo "Usage: pre_proc_init.csh SAT data.in"
    echo ""
    echo "       preprocess a set of images using default Dopper centroid, rear_range and radius, number of patches"
    echo "       a baseline-time plot will be generated"
    echo "       after running this command " 
    echo "       the user should choose an appropriate master image " 
    echo "       the user should also decide a common Dopper centroid, rear_range and radius, number of patches to run batch processing"
    echo "       the data with completely different Doppler centroid or baselines can be omitted from further processing"
    echo ""
    echo " SAT can be ALOS ERS or ENVI(ENVISAT)"
    echo ""
    echo "       format of data.in is:"
    echo "         line 1: master_name"
    echo "         line 2 and below: slave_name"
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
    echo "Example: pre_proc_init.csh ENVI data.in"
    echo ""
    exit 1
  endif

  set SAT = $1
  if ($SAT != ALOS && $SAT != ENVI && $SAT != ERS) then
    echo ""
    echo " SAT can be ALOS ERS or ENVI(ENVISAT)"
    echo ""
    exit 1
  endif
#
# open and read data.in table 
#
  echo ""
  echo "running pre_proc_init.csh"
  echo ""
  echo "preprocess master image"
  set line1 = `awk 'NR==1 {print $0}' $2`
  if ($SAT == ERS || $SAT == ENVI) then
    set master = $line1[1]
  else if ($SAT == ALOS) then
    set master = `echo $line1[1] | awk '{ print substr($1,8,length($1)-7)}'`
  endif
#
# unpack the master if necessary
#
  if ($SAT == ERS || $SAT == ENVI) then

    if(! -f $master.raw || ! -f $master.LED) then
      $1_pre_process $master 0 0 0 
    endif

    echo "before baseline"
    baseline_table.csh $master.PRM $master.PRM >! baseline_table.dat
    baseline_table.csh $master.PRM $master.PRM GMT >! table.gmt
    echo "after baseline"
 
  else if ($SAT == ALOS) then

    if(! -f IMG-HH-$master.raw || ! -f IMG-HH-$master.PRM ) then
      $1_pre_process IMG-HH-$master LED-$master
    endif

    set NEAR = `grep near_range IMG-HH-$master.PRM | awk '{print $3}'`
    set RAD = `grep earth_radius IMG-HH-$master.PRM | awk '{print $3}'`
    set rng_samp_rate_m = `grep rng_samp_rate IMG-HH-$master.PRM | awk 'NR == 1 {printf("%d", $3)}'`
    set FD1 = `grep fd1 IMG-HH-$master.PRM | awk '{print $3}'`
  
    baseline_table.csh IMG-HH-$master.PRM IMG-HH-$master.PRM >! baseline_table.dat
    baseline_table.csh IMG-HH-$master.PRM IMG-HH-$master.PRM GMT >! table.gmt

  endif

#
# loop and unpack the slave image using the same earth radius and near range as the master image
#
  foreach line2 (`awk 'NR>1 {print $0}' $2`)
    echo "pre_proc_batch.csh"
    echo "preprocess slave images"
    if ($SAT == ERS || $SAT == ENVI) then

      set slave = $line2  
      if(! -f $slave.raw || ! -f $slave.LED) then
        $1_pre_process $slave 0 0 0 
      endif
      baseline_table.csh $master.PRM $slave.PRM >> baseline_table.dat
      baseline_table.csh $master.PRM $slave.PRM GMT >> table.gmt
    else if ($SAT == ALOS) then 

      set slave = ` echo $line2 | awk '{ print substr($1,8,length($1)-7)}'`
      if(! -f IMG-HH-$slave.raw || ! -f IMG-HH-$slave.PRM ) then
        $1_pre_process IMG-HH-$slave LED-$slave 
      endif
#
# check the range sampling rate of the slave images and do conversion if necessary
#
      set rng_samp_rate_s = `grep rng_samp_rate IMG-HH-$slave.PRM | awk 'NR == 1 {printf("%d", $3)}'`
      set t = `echo $rng_samp_rate_m $rng_samp_rate_s | awk '{printf("%1.1f\n", $1/$2)}'`
      if ($t == 1.0) then
        echo "The range sampling rate for master and slave images are: "$rng_samp_rate_m
#
        baseline_table.csh IMG-HH-$master.PRM IMG-HH-$slave.PRM >> baseline_table.dat
        baseline_table.csh IMG-HH-$master.PRM IMG-HH-$slave.PRM GMT >> table.gmt
      else if ($t == 2.0) then
        echo "Convert the slave image from FBD to FBS mode"
        ALOS_fbd2fbs IMG-HH-$slave.PRM IMG-HH-$slave"_"FBS.PRM
#
        baseline_table.csh IMG-HH-$master.PRM IMG-HH-$slave"_"FBS.PRM >> baseline_table.dat
        baseline_table.csh IMG-HH-$master.PRM IMG-HH-$slave"_"FBS.PRM GMT >> table.gmt
        echo "Overwriting the old slave image"
        mv IMG-HH-$slave"_"FBS.PRM IMG-HH-$slave.PRM
        update_PRM.csh IMG-HH-$slave.PRM input_file IMG-HH-$slave.raw
        mv IMG-HH-$slave"_"FBS.raw IMG-HH-$slave.raw

      else if  ($t == 0.5) then
	echo "Use the FBS mode image as master"
        exit 1
      else
        echo "The range sampling rate for master and slave images are not convertable"
        exit 1
      endif    
# end of the if ($SAT == ALOS) 
    endif 
# end of the loop over slave images
  end
#
# make baseline plots 
#
  
  if ($SAT == ERS || $SAT == ENVI) then 
    awk '{print 1992+$1/365.25,$2,$7}' < table.gmt > text
#    awk '{print 1992+$1/365.25,$2,7,$4,$5,$6,$7}' < table.gmt > text
#    awk '{print 1992+$1/365.25,$2,7,-45,$5,$6,$7}' < table.gmt > text
  else if ($SAT == ALOS) then
    awk '{print 2006.5+($1-181)/365.25,$2,$7}' < table.gmt > text
#    awk '{print 2006.5+($1-181)/365.25,$2,9,$4,$5,$6,$7}' < table.gmt > text
  endif
  set region = `gmt gmtinfo text -C | awk '{print $1-0.5, $2+0.5, $3-500, $4+500}'`
#  set region = `minmax text -C | awk '{print $1-0.5, $2+0.5, -1200, 1200}'`
  gmt pstext text -JX8.8i/6.8i -R$region[1]/$region[2]/$region[3]/$region[4] -D0.2/0.2 -X1.5i -Y1i -K -N -F+f8,Helvetica+j5 > stacktable_all.ps 
  awk '{print $1,$2}' < text > text2
  gmt psxy text2 -Sp0.2c -G0 -R -JX -Ba1:"year":/a200g00f100:"baseline (m)":WSen -O >> stacktable_all.ps 

#
# clean up the mess
#
  rm text text2 
