#!/bin/csh -f
#
#  David Dandwell, December 29, 2015
#
# process Sentinel-1A TOPS data
# Automatically process a single frame of interferogram.
# see instruction.txt for details.
#

alias rm 'rm -f'
unset noclobber
#
  if ($#argv < 3) then
    echo ""
    echo "Usage: p2p_S1A_TOPS.csh master_image slave_image configuration_file"
    echo ""
    echo "Example: p2p_S1A_TOPS.csh S1A20150526_F1 S1A20150607_F1 config.tsx.slc.txt "
    echo ""
    echo "         Place the pre-processed data in a directory called raw and a dem.grd file in "
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
 if((! -f raw/$1.PRM) || (! -f raw/$1.LED) || (! -f raw/$1.SLC)) then
   echo " missing input files  raw/"$1
   exit
 endif
 if((! -f raw/$2.PRM) || (! -f raw/$2.LED) || (! -f raw/$2.SLC)) then
   echo " missing input files  raw/"$2
   exit
 endif
  if(! -f $3 ) then
    echo " no configure file: "$3
    exit
  endif
# 
# read parameters from configuration file
# 
  set stage = `grep proc_stage $3 | awk '{print $3}'`
  set earth_radius = `grep earth_radius $3 | awk '{print $3}'`
  if ((! $?earth_radius) || ($earth_radius == "")) then
    set earth_radius = 0
  endif
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
  set range_dec = `grep range_dec $3 | awk '{print $3}'`
  set azimuth_dec = `grep azimuth_dec $3 | awk '{print $3}'`
#
# read file names of raw data
#
  set master = $1 
  set slave = $2 

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
# preprocess the raw data
#
    echo " "
    echo "PREPROCESS - START"
    cd raw
#
# preprocess the raw data make the raw data and copy the PRM to PRM00
# in case the script is run a second time
#
#   make_raw.com
#
    if(-e $master.PRM00) then
       cp $master.PRM00 $master.PRM
       cp $slave.PRM00 $slave.PRM
    else
       cp $master.PRM $master.PRM00
       cp $slave.PRM $slave.PRM00
    endif
#
# set the num_lines to be the min of the master and slave
#
    @ m_lines  = `grep num_lines ../raw/$master.PRM | awk '{printf("%d",int($3))}' `
    @ s_lines  = `grep num_lines ../raw/$slave.PRM | awk '{printf("%d",int($3))}' `
    if($s_lines <  $m_lines) then
      update_PRM.csh $master.PRM num_lines $s_lines
      update_PRM.csh $master.PRM num_valid_az $s_lines
      update_PRM.csh $master.PRM nrows $s_lines
    else
      update_PRM.csh $slave.PRM num_lines $m_lines
      update_PRM.csh $slave.PRM num_valid_az $m_lines
      update_PRM.csh $slave.PRM nrows $m_lines
    endif
#
#   set the higher Doppler terms to zerp to be zero
#
    update_PRM.csh $master.PRM fdd1 0
    update_PRM.csh $master.PRM fddd1 0
#
    update_PRM.csh $slave.PRM fdd1 0
    update_PRM.csh $slave.PRM fddd1 0
#
    rm *.log
    rm *.PRM0

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
    ln -s ../raw/$master.SLC .
    ln -s ../raw/$slave.SLC .
    ln -s ../raw/$master.LED . 
    ln -s ../raw/$slave.LED .
    
#    cp $slave.PRM $slave.PRM0
#    resamp $master.PRM $slave.PRM $slave.PRMresamp $slave.SLCresamp 1
#    rm $slave.SLC
#    mv $slave.SLCresamp $slave.SLC
#    cp $slave.PRMresamp $slave.PRM
    cd ..
    echo "ALIGN - END"
  endif

##################################
# 3 - start from make topo_ra    #
##################################
#if (6 == 9) then
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
      cp ../SLC/$master.PRM master.PRM 
      ln -s ../raw/$master.LED . 
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
#
#  make sure the range increment of the amplitude image matches the topo_ra.grd
#
        set rng = `grdinfo topo/topo_ra.grd | grep x_inc | awk '{print $7}'`
        cd SLC 
        echo " range decimation is:  " $rng
        slc2amp.csh $master.PRM $rng amp-$master.grd
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
    set ref_id  = `grep SC_clock_start ../raw/$master.PRM | awk '{printf("%d",int($3))}' `
    set rep_id  = `grep SC_clock_start ../raw/$slave.PRM | awk '{printf("%d",int($3))}' `
    mkdir $ref_id"_"$rep_id
    cd $ref_id"_"$rep_id
    ln -s ../../raw/$ref.LED . 
    ln -s ../../raw/$rep.LED .
    ln -s ../../SLC/$ref.SLC . 
    ln -s ../../SLC/$rep.SLC .
    cp ../../SLC/$ref.PRM . 
    cp ../../SLC/$rep.PRM .

    if($topo_phase == 1) then
      if ($shift_topo == 1) then
        ln -s ../../topo/topo_shift.grd .
        intf.csh $ref.PRM $rep.PRM -topo topo_shift.grd  
        filter.csh $ref.PRM $rep.PRM $filter $dec $range_dec $azimuth_dec 
      else 
        ln -s ../../topo/topo_ra.grd . 
        intf.csh $ref.PRM $rep.PRM -topo topo_ra.grd 
        filter.csh $ref.PRM $rep.PRM $filter $dec $range_dec $azimuth_dec
      endif
    else
      intf.csh $ref.PRM $rep.PRM
      filter.csh $ref.PRM $rep.PRM $filter $dec $range_dec $azimuth_dec
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
      set ref_id  = `grep SC_clock_start ../SLC/$master.PRM | awk '{printf("%d",int($3))}' `
      set rep_id  = `grep SC_clock_start ../SLC/$slave.PRM | awk '{printf("%d",int($3))}' `
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
      
      snaphu.csh $threshold_snaphu $defomax $region_cut

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
    if ($threshold_geocode != 0 ) then
      cd intf
      set ref_id  = `grep SC_clock_start ../SLC/$master.PRM | awk '{printf("%d",int($3))}' `
      set rep_id  = `grep SC_clock_start ../SLC/$slave.PRM | awk '{printf("%d",int($3))}' `
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
    else
      echo ""
      echo "SKIP GEOCODE"
      echo ""
    endif
  endif

# end

