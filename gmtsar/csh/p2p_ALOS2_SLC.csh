#!/bin/csh -f
#       $Id$
#
#  Xiaopeng Tong, Feb 10, 2010
#
# Automatically process a single frame of interferogram.
# see instruction.txt for details.
#
# Modification History
#
# Xiaohua Xu, Oct 9, 2013
#
# Changed some code in Part 2 to make it work for one FBD and one FBS files.


alias rm 'rm -f'
unset noclobber
#
  if ($#argv != 3) then
    echo ""
    echo "Usage: p2p_ALOS_SLC.csh master_image slave_image configuration_file "
    echo ""
    echo "Example: p2p_ALOS_SLC.csh IMG-ALPSRP022200660-H1.1__A IMG-ALPSRP028910660-H1.1__A config.alos.slc.txt"
    echo ""
    echo "         Place the raw L1.1 data in a directory called raw and a dem.grd file in "
    echo "         a parallel directory called topo.  Execute this command at the directory"
    echo "         location above raw and topo.  The file dem.grd"
    echo "         is a dem that completely covers the SAR frame - larger is OK."
    echo "         If the dem is omitted then an interferogram will still be created"
    echo "         but there will not be geocoded output."
    echo "         A custom dem.grd can be made at the web site http://topex.ucsd.edu/gmtsar"
    echo ""
    echo ""
    exit 1
  endif

# start 

#
#   make sure the files exist
#
  if(! -f $3 ) then
    echo " no configure file: "$3
    exit
  endif
# 
# read parameters from configuration file
# 
  set SAT = `grep SAT $3 | awk '{print $3}'`
  set stage = `grep proc_stage $3 | awk '{print $3}'`
  set num_patches = `grep num_patches $3 | awk '{print $3}'`
  set SLC_factor = `grep SLC_factor $3 | awk '{print $3}'`
  set earth_radius = `grep earth_radius $3 | awk '{print $3}'`
  set topo_phase = `grep topo_phase $3 | awk '{print $3}'`
  set shift_topo = `grep shift_topo $3 | awk '{print $3}'`
  set switch_master = `grep switch_master $3 | awk '{print $3}'`
#
# if filter wavelength is not set then use a default of 200m
#
  set filter = `grep filter_wavelength $3 | awk '{print $3}'`
  if ( "x$filter" == "x" ) then
  set filter = 200
  echo " "
  echo "WARNING filter wavelength was not set in config.txt file"
  echo "        please specify wavelength (e.g., filter_wavelength = 200)"
  echo "        remove filter1 = gauss_alos_200m"
  endif
  echo $filter
  set dec = `grep dec_factor $3 | awk '{print $3}'` 
  set threshold_snaphu = `grep threshold_snaphu $3 | awk '{print $3}'`
  set threshold_geocode = `grep threshold_geocode $3 | awk '{print $3}'`
  set region_cut = `grep region_cut $3 | awk '{print $3}'`
  set switch_land = `grep switch_land $3 | awk '{print $3}'`
  set defomax = `grep defomax $3 | awk '{print $3}'`
  set near_interp = `grep near_interp $3 | awk '{print $3}'`
  
#
# read file names of raw data
#
  set master = ` echo $1 | awk '{ print substr($1,5,length($1)-4)}'`
  set slave =  ` echo $2 | awk '{ print substr($1,5,length($1)-4)}'`
  set master_led = ` echo $1 | awk '{ print substr($1,8,length($1)-7)}'`
  set slave_led =  ` echo $2 | awk '{ print substr($1,8,length($1)-7)}'`

  if ($switch_master == 0) then
    set ref = $master
    set rep = $slave
  else if ($switch_master == 1) then
    set ref = $slave
    set rep = $master
  else
    echo "Wrong paramter: switch_master "$switch_master
  endif
#
# make working directories
#  
  mkdir -p intf/ SLC/

#############################
# 1 - start from preprocess #
#############################

  if ($stage == 1) then
#
# first clean up 
# 
    rm raw/*.PRM*
    rm raw/*.SLC
    rm raw/*.LED
# 
# preprocess the raw data
#
    echo " "
    echo "PREPROCESS - START"
    cd raw
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
#
# make sure the files exits
#
    echo $commandline
    if(! -f $1 ) then
      echo " no file  "$1
      exit
    endif
    if(! -f $2 ) then
      echo " no file  "$2
      exit
    endif
#   ALOS_pre_process_SLC IMG-$master LED-$master_led $commandline -rbias -68. 
#   ALOS_pre_process_SLC IMG-$slave LED-$slave_led $commandline -rbias -68. 
    ALOS_pre_process_SLC IMG-$master LED-$master_led $commandline 
    ALOS_pre_process_SLC IMG-$slave LED-$slave_led $commandline 

#
# special code to filter the SLC data in range does not seem to work
#
    cd ..
    echo "PREPROCESS.CSH - END"
  endif

#############################################
# 2 - start from focus and align SLC images #
#############################################
  
  if ($stage <= 2) then
# 
# clean up 
#
    cleanup.csh SLC
#
# align SLC images 
# 
    echo " "
    echo "ALIGN - START"
    cd SLC
    cp ../raw/*.PRM .
    ln -s ../raw/IMG-$master.SLC . 
    ln -s ../raw/IMG-$slave.SLC . 
    ln -s ../raw/IMG-$master.LED . 
    ln -s ../raw/IMG-$slave.LED . 
    

# if images do not match, convert the FBD image to FBS
    set rng_samp_rate_m = `grep rng_samp_rate IMG-$master.PRM | awk 'NR == 1 {printf("%d", $3)}'`
    set rng_samp_rate_s = `grep rng_samp_rate IMG-$slave.PRM | awk 'NR == 1 {printf("%d", $3)}'`
    set t = `echo $rng_samp_rate_m $rng_samp_rate_s | awk '{printf("%1.1f\n", $1/$2)}'`
    if ($t == 1.0) then
        echo "The range sampling rate for master and slave images are: "$rng_samp_rate_m
    else if ($t == 2.0) then
        echo "Convert the slave image from FBD to FBS mode"
	ALOS_fbd2fbs_SLC IMG-$slave.PRM IMG-$slave"_"FBS.PRM
        echo "Overwriting the old slave image"
        mv IMG-$slave"_"FBS.PRM IMG-$slave.PRM
	update_PRM IMG-$slave.PRM input_file IMG-$slave.SLC
        mv IMG-$slave"_"FBS.SLC IMG-$slave.SLC
    else if  ($t == 0.5) then
	echo "Convert the master image from FBD to FBS mode"
	ALOS_fbd2fbs_SLC IMG-$master.PRM IMG-$master"_"FBS.PRM
        echo "Overwriting the old master image"
        mv IMG-$master"_"FBS.PRM IMG-$master.PRM
	update_PRM IMG-$master.PRM input_file IMG-$master.SLC
        mv IMG-$master"_"FBS.SLC IMG-$master.SLC
    else
	echo "The range sampling rate for master and slave images are not convertable"
    	exit 1
    endif

    cp IMG-$slave.PRM IMG-$slave.PRM0
    SAT_baseline IMG-$master.PRM IMG-$slave.PRM0 >> IMG-$slave.PRM
    xcorr IMG-$master.PRM IMG-$slave.PRM -xsearch 64 -ysearch 64 -nx 32 -ny 64
    fitoffset.csh 2 2 freq_xcorr.dat 18 >> IMG-$slave.PRM
    resamp IMG-$master.PRM IMG-$slave.PRM IMG-$slave.PRMresamp IMG-$slave.SLCresamp 4
    rm IMG-$slave.SLC
    mv IMG-$slave.SLCresamp IMG-$slave.SLC
    cp IMG-$slave.PRMresamp IMG-$slave.PRM
        
    cd ..
    echo "ALIGN - END"
  endif

##################################
# 3 - start from make topo_ra    #
##################################

  if ($stage <= 3) then
#
# clean up
#
    cleanup.csh topo
#
# make topo_ra if there is dem.grd
#
    if ($topo_phase == 1) then 
      echo " "
      echo "DEM2TOPO_RA.CSH - START"
      echo "USER SHOULD PROVIDE DEM FILE"
      cd topo
      cp ../SLC/IMG-$master.PRM master.PRM 
      ln -s ../raw/IMG-$master.LED . 
      if (-f dem.grd) then 
        dem2topo_ra.csh master.PRM dem.grd 
      else 
        echo "no DEM file found: " dem.grd 
        exit 1
      endif
      cd .. 
      echo "DEM2TOPO_RA.CSH - END"
# 
# shift topo_ra
# 
      if ($shift_topo == 1) then 
        echo " "
        echo "OFFSET_TOPO - START"
        cd SLC 
        slc2amp.csh IMG-$master.PRM 2 amp-$master.grd 
        cd ..
        cd topo
        ln -s ../SLC/amp-$master.grd . 
        offset_topo amp-$master.grd topo_ra.grd 0 0 7 topo_shift.grd 
        cd ..
        echo "OFFSET_TOPO - END"
      else if ($shift_topo == 0) then 
        echo "NO TOPO_RA SHIFT "
      else 
        echo "Wrong paramter: shift_topo "$shift_topo
        exit 1
      endif

      else if ($topo_phase == 0) then 
      echo "NO TOPO_RA IS SUBSTRACTED"
    else 
      echo "Wrong paramter: topo_phase "$topo_phase
      exit 1
    endif
  endif

##################################################
# 4 - start from make and filter interferograms  #
##################################################

  if ($stage <= 4) then
#
# clean up
#
    cleanup.csh intf
# 
# make and filter interferograms
# 
    echo " "
    echo "INTF.CSH, FILTER.CSH - START"
    cd intf/
    set ref_id  = `grep SC_clock_start ../SLC/IMG-$ref.PRM | awk '{printf("%d",int($3))}' `
    set rep_id  = `grep SC_clock_start ../SLC/IMG-$rep.PRM | awk '{printf("%d",int($3))}' `
    mkdir $ref_id"_"$rep_id
    cd $ref_id"_"$rep_id
    ln -s ../../raw/IMG-$ref.LED . 
    ln -s ../../raw/IMG-$rep.LED .
    ln -s ../../SLC/IMG-$ref.SLC . 
    ln -s ../../SLC/IMG-$rep.SLC .
    cp ../../SLC/IMG-$ref.PRM . 
    cp ../../SLC/IMG-$rep.PRM .
    if($topo_phase == 1) then
      if ($shift_topo == 1) then
        ln -s ../../topo/topo_shift.grd .
        intf.csh IMG-$ref.PRM IMG-$rep.PRM -topo topo_shift.grd  
        filter.csh IMG-$ref.PRM IMG-$rep.PRM $filter $dec 
      else 
        ln -s ../../topo/topo_ra.grd . 
        intf.csh IMG-$ref.PRM IMG-$rep.PRM -topo topo_ra.grd 
        filter.csh IMG-$ref.PRM IMG-$rep.PRM $filter $dec 
      endif
    else
      intf.csh IMG-$ref.PRM IMG-$rep.PRM
      filter.csh IMG-$ref.PRM IMG-$rep.PRM $filter $dec 
    endif
    cd ../..
    echo "INTF.CSH, FILTER.CSH - END"
  endif

################################
# 5 - start from unwrap phase  #
################################

  if ($stage <= 5 ) then
    if ($threshold_snaphu != 0 ) then
      cd intf
      set ref_id  = `grep SC_clock_start ../SLC/IMG-$ref.PRM | awk '{printf("%d",int($3))}' `
      set rep_id  = `grep SC_clock_start ../SLC/IMG-$rep.PRM | awk '{printf("%d",int($3))}' `
      cd $ref_id"_"$rep_id
      if ((! $?region_cut) || ($region_cut == "")) then
        set region_cut = `gmt grdinfo phase.grd -I- | cut -c3-20`
      endif

#
# landmask
#
      if ($switch_land == 1) then
        cd ../../topo
        if (! -f landmask_ra.grd) then
          landmask.csh $region_cut
        endif
        cd ../intf
        cd $ref_id"_"$rep_id
        ln -s ../../topo/landmask_ra.grd .
      endif

      echo " "
      echo "SNAPHU.CSH - START"
      echo "threshold_snaphu: $threshold_snaphu"
      
      if ($near_interp == 1) then
        snaphu_interp.csh $threshold_snaphu $defomax $region_cut
      else
        snaphu.csh $threshold_snaphu $defomax $region_cut
      endif

      echo "SNAPHU.CSH - END"
      cd ../..
    else 
      echo ""
      echo "SKIP UNWRAP PHASE"
    endif
  endif

###########################
# 6 - start from geocode  #
###########################

  if ($stage <= 6) then
    cd intf
    set ref_id  = `grep SC_clock_start ../SLC/IMG-$ref.PRM | awk '{printf("%d",int($3))}' `
    set rep_id  = `grep SC_clock_start ../SLC/IMG-$rep.PRM | awk '{printf("%d",int($3))}' `
    cd $ref_id"_"$rep_id
    echo " "
    echo "GEOCODE.CSH - START"
    rm raln.grd ralt.grd
    if ($topo_phase == 1) then
      rm trans.dat
      ln -s  ../../topo/trans.dat . 
      echo "threshold_geocode: $threshold_geocode"
      geocode.csh $threshold_geocode
    else 
      echo "topo_ra is needed to geocode"
      exit 1
    endif
    echo "GEOCODE.CSH - END"
    cd ../..
  endif

# end

