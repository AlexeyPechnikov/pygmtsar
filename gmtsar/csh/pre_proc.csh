#!/bin/csh -f
#       $Id$
#
#  Xiaopeng Tong, Mar 2 2010 
#  modified by D. Sandwell MAR 11 2010
#

alias rm 'rm -f'
unset noclobber

#
# check the number of arguments 
# 
  if (!($#argv == 3 || $#argv == 5 || $#argv == 7 || $#argv == 9 || $#argv == 11)) then 
    echo ""
    echo "Usage: pre_proc.csh SAT master_name slave_name [-near near_range] [-radius RE] [-npatch num_patches] [-fd1 DOPP]"
    echo ""
    echo "Example: pre_proc.csh ALOS IMG-HH-ALPSRP099496420-H1.0__A IMG-HH-ALPSRP220276420-H1.0__A"
    echo ""
    exit 1
  endif
#
# parse the command line arguments 
#
  echo "pre_proc.csh"
  set SAT = $1
  if ($SAT != ALOS) then
    echo ""
    echo " SAT must be ALOS - ERS and ENVISAT not yet implemented"
    echo ""
    exit 1
  endif
  set master = ` echo $2 | awk '{ print substr($1,8,length($1)-7)}'`
  set slave = ` echo $3 | awk '{ print substr($1,8,length($1)-7)}'`
  if (! -f IMG-HH-$master || ! -f IMG-HH-$slave || ! -f LED-$master || ! -f LED-$slave) then 
    echo ""
    echo "Error : Can not find input file at current directory!"   
    echo ""
    exit 1
  endif
#
# unpack the master if necessary 
#
  if (! -f IMG-HH-$master.raw || ! -f IMG-HH-$master.PRM ) then 
    echo "pre_process master image"
    $1_pre_process IMG-HH-$master LED-$master $argv[4-$#argv]
  endif

    set NEAR = `grep near_range IMG-HH-$master.PRM | awk '{print $3}'`
    set RAD = `grep earth_radius IMG-HH-$master.PRM | awk '{print $3}'`
    set rng_samp_rate_m = `grep rng_samp_rate IMG-HH-$master.PRM | awk 'NR == 1 {printf("%d", $3)}'`
    set FD1 = `grep fd1 IMG-HH-$master.PRM | awk '{print $3}'`
    set num_patches = `grep num_patches IMG-HH-$master.PRM | awk '{print $3}'`
#
# unpack the slave image using the same earth radius and near range as the master image
#
  echo "pre_process slave image"
  $1_pre_process IMG-HH-$slave LED-$slave -fd1 $FD1 -near $NEAR -radius $RAD -npatch $num_patches -fd1 $FD1
#
# check the range sampling rate of the slave images and do conversion if necessary
#
  set rng_samp_rate_s = `grep rng_samp_rate IMG-HH-$slave.PRM | awk 'NR == 1 {printf("%d", $3)}'`
  set t = `echo $rng_samp_rate_m $rng_samp_rate_s | awk '{printf("%1.1f\n", $1/$2)}'`

  if ($t == 1.0) then
    echo "The range sampling rate for master and slave images are: "$rng_samp_rate_m
  else if ($t == 2.0) then
    echo "Convert the slave image from FBD to FBS mode"
    ALOS_fbd2fbs IMG-HH-$slave.PRM IMG-HH-$slave"_"FBS.PRM
    echo "Overwriting the old slave image"
    mv IMG-HH-$slave"_"FBS.PRM IMG-HH-$slave.PRM
    update_PRM.csh IMG-HH-$slave.PRM input_file IMG-HH-$slave.raw
    mv IMG-HH-$slave"_"FBS.raw IMG-HH-$slave.raw
  else if  ($t == 0.5) then
    echo "Use FBS mode image as master"
    exit 1
  else
    echo "The range sampling rate for master and slave images are not convertable"
    exit 1
  endif    
#
