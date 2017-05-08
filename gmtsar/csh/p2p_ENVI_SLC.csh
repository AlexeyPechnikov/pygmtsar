#!/bin/csh -f
#       $Id$
#
#  Xiaopeng Tong, Jan 14, 2014
#  Modified by Anders Hogrelius, April 18, 2017
#
# process envisat formatted  L1 data
# Automatically process a single frame of interferogram.
# see instruction.txt for details.
#


alias rm 'rm -f'
unset noclobber
#
  if ($#argv != 3) then
    echo ""
    echo "Usage: p2p_ENVI_SLC.csh master_image slave_image configuration_file "
    echo ""
    echo "Example: p2p_ENVI_SLC.csh ASA_IMS_1PNESA20060130_054255_000000182044_00392_20486_0000 ASA_IMS_1PNESA20070709_054254_000000172059_00392_28001_0000 config.envi.slc.txt"
    echo ""
    echo "         Place the L.1 data in a directory called raw and a dem.grd file in "
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
#
#
 if(! -f $3 ) then
   echo " no configure file: "$3
   exit
 endif


# 
# read parameters from configuration file
# 
  set stage = `grep proc_stage $3 | awk '{print $3}'`
  set near_range = `grep near_range $3 | awk '{print $3}'`
  if ((! $?near_range) || ($near_range == "")) then
    set near_range = 0
  endif
  set earth_radius = `grep earth_radius $3 | awk '{print $3}'`
  if ((! $?earth_radius) || ($earth_radius == "")) then
    set earth_radius = 0
  endif
  set npatch = `grep num_patches $3 | awk '{print $3}'`
  if ((! $?npatch) || ($npatch == "")) then
    set npatch = 0
  endif
  set fd = `grep fd1 $3 | awk '{print $3}'`
  set topo_phase = `grep topo_phase $3 | awk '{print $3}'`
  set shift_topo = `grep shift_topo $3 | awk '{print $3}'`
  set switch_master = `grep switch_master $3 | awk '{print $3}'`
#
# if filter wavelength is not set then use a default of 100m
#
  set filter = `grep filter_wavelength $3 | awk '{print $3}'`
  if ( "x$filter" == "x" ) then
  set filter = 100
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
# first clean up
#
    cleanup.csh raw
#
# preprocess the raw data
#
    echo ""
    echo " PREPROCESS Envisat SLC DATA  -- START"
    cd raw
    echo "ENVI_SLC_pre_process $master $earth_radius "
    ENVI_SLC_pre_process $master $earth_radius 
    set NEAR = `grep near_range $master.PRM | awk '{print $3}'`
    set RAD = `grep earth_radius $master.PRM | awk '{print $3}'`
    echo "ENVI_SLC_pre_process $slave $RAD "
    ENVI_SLC_pre_process $slave $RAD 
#
#   check patch number, if different, use the smaller one
#
    set pch1 = `grep patch $master.PRM | awk '{printf("%d ",$3)}'`
    set pch2 = `grep patch $slave.PRM | awk '{printf("%d ",$3)}'`
    echo "Different number of patches: $pch1 $pch2"
    if ($pch1 != $pch2) then
      if ($pch1 < $pch2) then
        update_PRM.csh $slave.PRM num_patches $pch1
        echo "Number of patches is set to $pch1"
      else
        update_PRM.csh $master.PRM num_patches $pch2
        echo "Number of patches is set to $pch2"
      endif
    endif
  endif
#


#############################################
# 2 - start from align SLC images #
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
    
    cp $slave.PRM $slave.PRM0
    SAT_baseline $master.PRM $slave.PRM0 >> $slave.PRM
    xcorr $master.PRM $slave.PRM -xsearch 128 -ysearch 128
    fitoffset.csh 2 2 freq_xcorr.dat >> $slave.PRM
    resamp $master.PRM $slave.PRM $slave.PRMresamp $slave.SLCresamp 4
    rm $slave.SLC
    mv $slave.SLCresamp $slave.SLC
    cp $slave.PRMresamp $slave.PRM
        
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
        set rng = `gmt grdinfo topo/topo_ra.grd | grep x_inc | awk '{print $7}'`
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
        filter.csh $ref.PRM $rep.PRM $filter $dec 
      else 
        ln -s ../../topo/topo_ra.grd . 
        intf.csh $ref.PRM $rep.PRM -topo topo_ra.grd 
        filter.csh $ref.PRM $rep.PRM $filter $dec 
      endif
    else
      intf.csh $ref.PRM $rep.PRM
      filter.csh $ref.PRM $rep.PRM $filter $dec 
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
  endif

# end

