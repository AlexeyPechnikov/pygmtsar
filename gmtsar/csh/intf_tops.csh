#!/bin/csh -f
#       $Id$

# generate interferograms for tops stacks
# used for time series analysis

# Xiaohua(Eric) Xu, Jan 20 2016
#

  if ($#argv != 2) then
    echo ""
    echo "Usage: intf_tops.csh intf.in batch_tops.config"
    echo "  generate interferograms for a set of tops images in intf.in, dem required in ./topo"
    echo "  supermaster's name required in batch_tops.config"
    echo ""
    echo "  format of data.in:"
    echo "    master_image_stem:slave_image_stem"
    echo ""
    echo "  example of data.in"
    echo "    S1A20150628_ALL_F1:S1A20150720_ALL_F1"
    echo "    S1A20150720_ALL_F1:S1A20150809_ALL_F1"
    echo ""
    echo "  outputs:"
    echo "    to ./intf_all"
    echo ""
    exit 1
  endif

#
# make sure the files exist
#
  if (! -f $1) then
    echo "no input file:" $2
    exit
  endif
  
  if (! -f $2) then
    echo "no config file:" $2
    exit
  endif
#
# read parameters from config file
#

  set stage = `grep proc_stage $2 | awk '{print $3}'`
  set master = `grep master_image $2 | awk '{print $3}'`
#
# if filter wavelength is not set then use a default of 200m
#
  set filter = `grep filter_wavelength $2 | awk '{print $3}'`
  if ( "x$filter" == "x" ) then
  set filter = 200
  echo " "
  echo "WARNING filter wavelength was not set in config.txt file"
  echo "        please specify wavelength (e.g., filter_wavelength = 200)"
  echo "        remove filter1 = gauss_alos_200m"
  endif
  set dec = `grep dec_factor $2 | awk '{print $3}'`
  set topo_phase = `grep topo_phase $2 | awk '{print $3}'`
  set shift_topo = `grep shift_topo $2 | awk '{print $3}'`
  set threshold_snaphu = `grep threshold_snaphu $2 | awk '{print $3}'`
  set threshold_geocode = `grep threshold_geocode $2 | awk '{print $3}'`
  set region_cut = `grep region_cut $2 | awk '{print $3}'`
  set switch_land = `grep switch_land $2 | awk '{print $3}'`
  set defomax = `grep defomax $2 | awk '{print $3}'`
  set range_dec = `grep range_dec $2 | awk '{print $3}'`
  set azimuth_dec = `grep azimuth_dec $2 | awk '{print $3}'`

##################################
# 1 - start from make topo_ra  #
##################################

if ($stage <= 1) then
#
# clean up
#
  cleanup.csh topo
#
# make topo_ra
#
  if ($topo_phase == 1) then
    echo " "
    echo "DEM2TOPOPHASE.CSH - START"
    echo "USER SHOULD PROVIDE DEM FILE"
    cd topo
    cp ../raw/$master.PRM ./master.PRM
    ln -s ../raw/$master.LED .
    if (-f dem.grd) then
      dem2topo_ra.csh master.PRM dem.grd
    else
      echo "no DEM file found: " dem.grd
      exit 1
    endif
    cd ..
    echo "DEM2TOPOPHASE.CSH - END"

#
#  shift topo_ra
#  
    if ($shift_topo == 1) then
      echo " "
      echo "OFFSET_TOPO - START"
      cd topo
      ln -s ../raw/$master.SLC .
      slc2amp.csh master.PRM 4 amp-$master.grd
      offset_topo amp-$master.grd topo_ra.grd 0 0 7 topo_shift.grd
      cd ..
      echo  "OFFSET_TOPO - END"
    else if ($shift_topo == 0) then
      echo "NO TOPOPHASE SHIFT "
    else
      echo "Wrong paramter: shift_topo "$shift_topo
      exit 1
    endif
  else if ($topo_phase == 0) then
    echo "NO TOPOPHASE IS SUBSTRACTED"
  else
    echo "Wrong paramter: topo_phase "$topo_phase
    exit 1
  endif
endif

##################################################
# 2 - start from make and filter interferograms  #
#                unwrap phase and geocode        #
##################################################

if ($stage <= 2) then
#
# make working directories
#
  echo ""
  echo "START FORM A STACK OF INTERFEROGRAMS"
  echo ""
  mkdir -p intf/
  mkdir -p intf_all/
#
# loop over intf.in
#
  foreach line (`awk '{print $0}' $1`)
    set ref = `echo $line | awk -F: '{print $1}'`
    set rep = `echo $line | awk -F: '{print $2}'`
    set ref_id  = `grep SC_clock_start ./raw/$ref.PRM | awk '{printf("%d",int($3))}' `
    set rep_id  = `grep SC_clock_start ./raw/$rep.PRM | awk '{printf("%d",int($3))}' `

    echo ""
    echo "INTF.CSH, FILTER.CSH - START"
    cd intf
    mkdir $ref_id"_"$rep_id
    cd $ref_id"_"$rep_id
    ln -s ../../raw/$ref.LED .
    ln -s ../../raw/$rep.LED .
    ln -s ../../raw/$ref.SLC .
    ln -s ../../raw/$rep.SLC .
    cp ../../raw/$ref.PRM .
    cp ../../raw/$rep.PRM .
    
    if($topo_phase == 1) then
      if($shift_topo == 1) then
        ln -s ../../topo/topo_shift.grd .
        intf.csh $ref.PRM $rep.PRM -topo topo_shift.grd
      else
        ln -s ../../topo/topo_ra.grd .
        intf.csh $ref.PRM $rep.PRM -topo topo_ra.grd
      endif
    else
      intf.csh $ref.PRM $rep.PRM
    endif
    filter.csh $ref.PRM $rep.PRM $filter $dec $range_dec $azimuth_dec
    echo "INTF.CSH, FILTER.CSH - END"

#
# unwrapping
#
    if ($region_cut == "") then
      set region_cut = `gmt grdinfo phase.grd -I- | cut -c3-20`
    endif
    
    if ($threshold_snaphu != 0 ) then
      if ($switch_land == 1) then
        cd ../../topo
        if (! -f landmask_ra.grd) then
          landmask.csh $region_cut
        endif
        cd ../intf
        cd $ref_id"_"$rep_id
        ln -s ../../topo/landmask_ra.grd .
      endif

      echo ""
      echo "SNAPHU.CSH - START"
      echo "threshold_snaphu: $threshold_snaphu"
      snaphu.csh $threshold_snaphu $defomax $region_cut
      echo "SNAPHU.CSH - END"
    else
      echo ""
      echo "SKIP UNWRAP PHASE"
    endif
#
# geocoding
#
    echo ""
    echo "GEOCODE.CSH - START"
    rm raln.grd ralt.grd
    if ($topo_phase == 1 && $threshold_geocode != 0) then
      rm trans.dat
      ln -s  ../../topo/trans.dat .
      echo "threshold_geocode: $threshold_geocode"
      geocode.csh $threshold_geocode
    else if($topo_phase == 1 && $threshold_geocode == 0) then
      echo "SKIP GEOCODING"
    else
      echo "topo_ra is needed to geocode"
      exit 1
    endif
  
  cd ../..
  if(-f intf_all/$ref_id"_"$rep_id) rm -rf intf_all/$ref_id"_"$rep_id 
  mv intf/$ref_id"_"$rep_id intf_all/$ref_id"_"$rep_id

  end
endif

echo ""
echo "END STACK OF TOPS INTERFEROGRAMS"
echo ""


