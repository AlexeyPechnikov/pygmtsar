#!/bin/csh -f
#       $Id$
#
#  Xiaohua Xu, Jan, 2018
#
#  Automatically perform two-path processing on raw(1.0)/SLC(1.1) data
#  

  if ($#argv != 3 && $#argv != 4) then
    echo ""
    echo "Usage: p2p_processing.csh SAT master_image slave_image [configuration_file] "
    echo ""
    echo "Example: p2p_processing.csh ALOS IMG-HH-ALPSRP055750660-H1.0__A IMG-HH-ALPSRP049040660-H1.0__A [config.alos.txt]"
    echo ""
    echo "    Put the data and orbit files in the raw folder, put DEM in the topo folder"
    echo "    The SAT needs to be specified, choices with in ERS, ENVI, ALOS, ALOS_SLC, ALOS2, ALOS2_Scan"
    echo "    S1_Strip, S1_TOPS, ENVI_SLC, CSK_Raw, CSK_SLC, TSX, RS2"
    echo ""
    echo "    Make sure the files from the same date have the same stem, e.g. aaaa.tif aaaa.xml aaaa.cos aaaa.EOF, etc"
    echo ""
    echo "    If the configuration file is left blank, the program will generate one "
    echo "    with default parameters "
    exit 1
  endif
    
# start 
#  Make sure the config exist
  if ($#argv == 4) then 
    if(! -f $4 ) then
      echo " no configure file: "$4
      echo " Leave it blank to generate config file with default values."
      exit 1
    endif
  endif
  
#
#  Read parameters from the configure file
#

  set SAT = `echo $1`
  if ($#argv == 4) then
    set conf = `echo $4`
  else
    pop_config.csh $SAT > config.$SAT.txt
    set conf = `echo "config.$SAT.txt"`
  endif
  # conf may need to be changed later on
  set stage = `grep proc_stage $conf | awk '{print $3}'`
  set num_patches = `grep num_patches $conf | awk '{print $3}'`
  set near_range = `grep near_range $conf | awk '{print $3}'`
  set earth_radius = `grep earth_radius $conf | awk '{print $3}'`
  set fd = `grep fd1 $conf | awk '{print $3}'`
  set topo_phase = `grep topo_phase $conf | awk '{print $3}'`
  set shift_topo = `grep shift_topo $conf | awk '{print $3}'`
  set switch_master = `grep switch_master $conf | awk '{print $3}'`
  set filter = `grep filter_wavelength $conf | awk '{print $3}'` 
  #if ( "x$filter" == "x" ) then 
  #  set filter = 200
  #  echo " "
  #  echo "WARNING filter wavelength was not set in config.txt file"
  #  echo "        please specify wavelength (e.g., filter_wavelength = 200)"
  #  echo "        remove filter1 = gauss_alos_200m"
  #endif
  set dec = `grep dec_factor $conf | awk '{print $3}'` 
  set threshold_snaphu = `grep threshold_snaphu $conf | awk '{print $3}'`
  set threshold_geocode = `grep threshold_geocode $conf | awk '{print $3}'`
  set region_cut = `grep region_cut $conf | awk '{print $3}'`
  set mask_water = `grep mask_water $conf | awk '{print $3}'`
  set defomax = `grep defomax $conf | awk '{print $3}'`
  set range_dec = `grep range_dec $conf | awk '{print $3}'`
  set azimuth_dec = `grep azimuth_dec $conf | awk '{print $3}'`
  set SLC_factor = `grep SLC_factor $conf | awk '{print $3}'`
  set master = ` echo $2 `
  set slave =  ` echo $3 `
  echo ""
  echo "Working on images $master $slave ..."
  

#
#  combine preprocess parameters
#  

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


#############################
# 1 - start from preprocess #
#############################
#
#   make sure the files exist
#
  if ($stage == 1) then
    echo ""
    echo "PREPROCESS - START"
    echo ""
    if ($SAT == "ALOS" || $SAT == "ALOS2" || $SAT == "ALOS_SLC" || $SAT == "ALOS2_Scan") then
      if(! -f raw/$master ) then
        echo " no file  raw/"$master
        exit
      endif
      if(! -f raw/$slave ) then
        echo " no file  raw/"$slave
        exit
      endif
    else if ($SAT == "ENVI_SLC") then
      if(! -f raw/$master.N1 ) then
        echo " no file  raw/"$master.N1
        exit
      endif
      if(! -f raw/$slave.N1 ) then
        echo " no file  raw/"$slave.N1
        exit
      endif
    else if ($SAT == "ERS") then
      if(! -f raw/$master.dat ) then
        echo " no file  raw/"$master.dat
        exit
      endif
      if(! -f raw/$slave.dat ) then
        echo " no file  raw/"$slave.dat
        exit
      endif
      if(! -f raw/$master.ldr ) then
        echo " no file  raw/"$master.ldr
        exit
      endif
      if(! -f raw/$slave.ldr ) then
        echo " no file  raw/"$slave.ldr
        exit
      endif
    else if ($SAT == "ENVI") then
      if(! -f raw/$master.baq ) then
        echo " no file  raw/"$master.baq
        exit
      endif
      if(! -f raw/$slave.baq ) then
        echo " no file  raw/"$slave.baq
        exit
      endif
    else if ($SAT == "S1_Strip" || $SAT == "S1_TOPS") then
      if(! -f raw/$master.xml ) then
        echo " no file  raw/"$master".xml"
        exit
      endif
      if(! -f raw/$master.tiff ) then
        echo " no file  raw/"$master".tiff"
        exit
      endif
      if(! -f raw/$slave.xml ) then
        echo " no file  raw/"$slave".xml"
        exit
      endif
      if(! -f raw/$slave.tiff ) then
        echo " no file  raw/"$slave".tiff"
        exit
      endif
      if ($SAT == "S1_TOPS") then
        if(! -f raw/$master.EOF ) then
          echo " no file  raw/"$master".EOF"
        endif
        if(! -f raw/$slave.EOF ) then
          echo " no file  raw/"$slave".EOF"
        endif
      endif
    else if ($SAT == "CSK_Raw" || $SAT == "CSK_SLC") then
      if(! -f raw/$master.h5 ) then
        echo " no file  raw/"$master".h5"
        exit
      endif
      if(! -f raw/$slave.h5 ) then
        echo " no file  raw/"$slave".h5"
        exit
      endif
    else if ($SAT == "RS2") then
      if(! -f raw/$master.xml ) then
        echo " no file  raw/"$master".xml"
        exit
      endif
      if(! -f raw/$master.tif ) then
        echo " no file  raw/"$master".tif"
        exit
      endif
      if(! -f raw/$slave.xml ) then
        echo " no file  raw/"$slave".xml"
        exit
      endif
      if(! -f raw/$slave.tif ) then
        echo " no file  raw/"$slave".tif"
        exit
      endif
    else if ($SAT == "TSX") then
      if(! -f raw/$master.xml ) then
        echo " no file  raw/"$master".xml"
        exit
      endif
      if(! -f raw/$slave.xml ) then
        echo " no file  raw/"$slave".xml"
        exit
      endif
      if(! -f raw/$master.cos ) then
        echo " no file  raw/"$master".cos"
        exit
      endif
      if(! -f raw/$slave.cos ) then
        echo " no file  raw/"$slave".cos"
        exit
      endif
    endif
  
#
#  Start preprocessing
#
  
  
    rm raw/*.PRM*
    rm raw/*.SLC
    rm raw/*.LED
    cd raw
    
    #echo "pre_proc.csh $SAT $master $slave $commandline"
    
    pre_proc.csh $SAT $master $slave $commandline   
        
    cd ..
    echo " "
    echo "PREPROCESS - END"
    echo ""
  endif
 
#############################################
# 2 - start from focus and align SLC images #
#############################################
# 

  mkdir -p intf SLC

  if ($stage <= 2) then 
    cleanup.csh SLC
#
# focus and align SLC images 
# 
    echo " "
    echo "ALIGN.CSH - START"
    echo ""
    cd SLC
    if ($SAT == "ERS" || $SAT == "ENVI" || $SAT == "ALOS" || $SAT == "CSK_Raw") then
      cp ../raw/*.PRM .
      ln -s ../raw/$master.raw . 
      ln -s ../raw/$slave.raw . 
      ln -s ../raw/$master.LED . 
      ln -s ../raw/$slave.LED . 
      if ($SAT == "ERS" || $SAT == "ENVI") then
        align.csh $SAT $master $slave 
      else
        align.csh SAT $master $slave
      endif
    else if ($SAT == "ALOS2" || $SAT == "ALOS_SLC" || $SAT == "ALOS2_Scan" || $SAT == "S1_Strip" || $SAT == "CSK_SLC" || $SAT == "RS2" || $SAT == "ENVI_SLC" || $SAT == "TSX") then
      cp ../raw/*.PRM .
      ln -s ../raw/$master.SLC . 
      ln -s ../raw/$slave.SLC . 
      ln -s ../raw/$master.LED . 
      ln -s ../raw/$slave.LED .
      cp $slave.PRM $slave.PRM0
      SAT_baseline $master.PRM $slave.PRM0 >> $slave.PRM
      if ($SAT == "ALOS2_Scan") then
        xcorr $master.PRM $slave.PRM -xsearch 32 -ysearch 256 -nx 32 -ny 128
        awk '{print $4}' < freq_xcorr.dat > tmp.dat
        set amedian = `sort -n tmp.dat | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'`
        set amax = `echo $amedian | awk '{print $1+3}'`
        set amin = `echo $amedian | awk '{print $1-3}'`
        awk '{if($4 > '$amin' && $4 < '$amax') print $0}' < freq_xcorr.dat > tmp2.dat
        fitoffset.csh 2 3 tmp2.dat >> $slave.PRM 10
        rm tmp2.dat
      else
        xcorr $master.PRM $slave.PRM -xsearch 128 -ysearch 128
        fitoffset.csh 2 2 freq_xcorr.dat 10 >> $slave.PRM
      endif
      resamp $master.PRM $slave.PRM $slave.PRMresamp $slave.SLCresamp 4
      rm $slave.SLC
      mv $slave.SLCresamp $slave.SLC
      cp $slave.PRMresamp $slave.PRM
    else if ($SAT == "S1_TOPS") then
      set master = `echo $master | awk '{ print "S1A"substr($1,16,8)"_"substr($1,25,6)"_F"substr($1,7,1)}'`
      set slave = `echo $slave | awk '{ print "S1A"substr($1,16,8)"_"substr($1,25,6)"_F"substr($1,7,1)}'`
      cp ../raw/*.PRM .
      ln -s ../raw/$master.SLC . 
      ln -s ../raw/$slave.SLC . 
      ln -s ../raw/$master.LED . 
      ln -s ../raw/$slave.LED .
    endif
    
    cd ..
    echo ""
    echo "ALIGN.CSH - END"
    echo ""  
  
  endif
##################################
# 3 - start from make topo_ra  #
##################################
#
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
        cd SLC 
        set rng_samp_rate = `grep rng_samp_rate $master.PRM | awk 'NR == 1 {printf("%d", $3)}'`
        set rng = `gmt grdinfo ../topo/topo_ra.grd | grep x_inc | awk '{print $7}'`
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
#
  if ($stage <= 4) then
#
# clean up
#
    cleanup.csh intf
#
# select the master
#    
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
# make and filter interferograms
# 
    echo " "
    echo "INTF.CSH, FILTER.CSH - START"
    cd intf/
    set ref_id  = `grep SC_clock_start ../raw/$ref.PRM | awk '{printf("%d",int($3))}' `
    set rep_id  = `grep SC_clock_start ../raw/$rep.PRM | awk '{printf("%d",int($3))}' `
    mkdir $ref_id"_"$rep_id
    cd $ref_id"_"$rep_id
    ln -s ../../SLC/$ref.LED . 
    ln -s ../../SLC/$rep.LED .
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
      echo "NO TOPOGRAPHIC PHASE REMOVAL PORFORMED"
      intf.csh $ref.PRM $rep.PRM
      filter.csh $ref.PRM $rep.PRM $filter $dec $range_dec $azimuth_dec
    endif
    cd ../..
    echo "INTF.CSH, FILTER.CSH - END"
  endif

################################
# 5 - start from unwrap phase  #
################################
#
  if ($stage <= 5 ) then
    if ($threshold_snaphu != 0 ) then
      cd intf
      set ref_id  = `grep SC_clock_start ../raw/$ref.PRM | awk '{printf("%d",int($3))}' `
      set rep_id  = `grep SC_clock_start ../raw/$rep.PRM | awk '{printf("%d",int($3))}' `
      cd $ref_id"_"$rep_id
      if ((! $?region_cut) || ($region_cut == "")) then
        set region_cut = `gmt grdinfo phase.grd -I- | cut -c3-20`
      endif
#
# landmask
#
      if ($mask_water == 1) then
        cd ../../topo
        if (! -f landmask_ra.grd) then
          landmask.csh $region_cut
        endif
        cd ../intf
        cd $ref_id"_"$rep_id
        ln -s ../../topo/landmask_ra.grd .
      endif
#
      echo " "
      echo "SNAPHU.CSH - START"
      echo "threshold_snaphu: $threshold_snaphu"
#
      snaphu.csh $threshold_snaphu $defomax $region_cut
#
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
#
  if ($stage <= 6) then
    if ($threshold_geocode != 0 ) then
      cd intf
      set ref_id  = `grep SC_clock_start ../raw/$ref.PRM | awk '{printf("%d",int($3))}' `
      set rep_id  = `grep SC_clock_start ../raw/$rep.PRM | awk '{printf("%d",int($3))}' `
      cd $ref_id"_"$rep_id
      echo " "
      echo "GEOCODE.CSH - START"
      #if (-f raln.grd) rm raln.grd 
      #if (-f ralt.grd) rm ralt.grd
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
#
# end  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
