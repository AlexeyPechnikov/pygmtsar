#!/bin/csh -f
#       $Id$
#
#  D. Sandwell FEB 4 2010
#  M. Wei MAY 4 2010 - ENVISAT
# 
# Align a slave image to a master image and check results
#
alias rm 'rm -f'
unset noclobber
#
# check the number of arguments 
# 
  if ($#argv < 3) then 
    echo ""
    echo "Usage: align.csh SAT master_name slave_name [supermaster_name]"
    echo ""
    echo " The supermaster_namestem is required if this is secondary alignment."
    echo " SAT = ERS or ENVI or ALOS  or generic SAT"
    echo ""
    echo "Example: align.csh ALOS IMG-HH-ALPSRP055750660-H1.0__A IMG-HH-ALPSRP049040660-H1.0__A "
    echo ""
    exit 1
  endif
  set SAT = $1
  if( ($SAT != ALOS) && ($SAT != ENVI) && ($SAT != ERS) && ($SAT != SAT)) then
    echo ""
    echo " SAT must be ERS, ENVI, or ALOS or generic SAT"
    echo ""
    exit 1
  endif
#
# focus the master if necessary
# Do it no matter what for now. Put SLC_file to PRM. Might not be necessary
#
  if(! -f $2.SLC) then
    echo "focussing master"
    sarp.csh $2.PRM 
  else
    update_PRM.csh $2.PRM SLC_file $2.SLC
  endif
#
# focus the slave image
#
# check the range sampling rate 
# 
  set rng_samp_rate_m = `grep rng_samp_rate $2.PRM | awk 'NR == 1 {printf("%d", $3)}'`
  set rng_samp_rate_s = `grep rng_samp_rate $3.PRM | awk 'NR == 1 {printf("%d", $3)}'`
  if ($rng_samp_rate_m != $rng_samp_rate_s) then 
    echo "The range sampling rate for master and slave differ"
    echo "Need to run the interferogram in steps until process2pass.csh is fixed"
    exit 1
  endif 
  echo "align.csh"
  echo "focusing slave"
  sarp.csh $3.PRM 
#
# get the starting alignment parameters and run xcorr
#
  cp $2.PRM $2.PRM0
  cp $3.PRM $3.PRM0
  if($#argv == 4) then
    set RSHIFT = `$1_baseline $4.PRM $3.PRM | grep rshift | awk '{print $3}'`
    set ASHIFT = `$1_baseline $4.PRM $3.PRM | grep ashift | awk '{print $3}'`
#
#   use the PRF of the supermaster in the surrogate master
#
    set PRF = `grep PRF $4.PRM | awk '{print $3}'`
    update_PRM.csh $2.PRM PRM $PRF
  else
    set RSHIFT = `$1_baseline $2.PRM $3.PRM | grep rshift | awk '{print $3}'`
    set ASHIFT = `$1_baseline $2.PRM $3.PRM | grep ashift | awk '{print $3}'`
  endif
  update_PRM.csh $3.PRM rshift $RSHIFT
  update_PRM.csh $3.PRM ashift $ASHIFT
  echo "align.csh"
  echo "correlate master and slave to find offset parameters"
  if( $SAT == "ERS") then
    xcorr $2.PRM $3.PRM -xsearch 128 -ysearch 128
  else
    xcorr $2.PRM $3.PRM -xsearch 64 -ysearch 64
  endif
#
  mv $3.SLC $3.SLC0
  mv $3.PRM junk.PRM
  cp $2.PRM0 $2.PRM
  grep -v shift < junk.PRM > $3.PRM
#
# put in the alignment parameters 
#
  fitoffset.csh 3 3 freq_xcorr.dat 18 >> $3.PRM
  mv freq_xcorr.dat xcorr_$2_$3.dat0
#
# refocus the second image
#
  echo "align.csh"
  echo "refocus slave"
  sarp.csh $3.PRM 
#
rm *SLC0
rm junk*
#done
