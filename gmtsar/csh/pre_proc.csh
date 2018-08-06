#!/bin/csh -f
#       $Id$
#
#  Xiaohua Xu, Jan, 2018 
#
#  Preprocess the data to Raw or SLC or aligned SLC depending on which satellite to use
#

alias rm 'rm -f'
unset noclobber

#
# check the number of arguments 
# 
  if (!($#argv == 3 || $#argv == 5 || $#argv == 7 || $#argv == 9 || $#argv == 11)) then 
    echo ""
    echo "Usage: pre_proc.csh SAT master_stem slave_stem [-near near_range] [-radius RE] [-npatch num_patches] [-fd1 DOPP]"
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
  echo "0" $argv[4-$#argv] | awk '{for (i=1;i<=12;i++) print $i}' > tmp
  set a = `grep near tmp -n | awk -F: '{print $1}'`
  set NEAR = `awk 'NR=='$a'+1{print $1}' tmp`
  set a = `grep near tmp -n | awk -F: '{print $1}'`
  set RAD = `awk 'NR=='$a'+1{print $1}' tmp`
  set a = `grep npatch tmp -n | awk -F: '{print $1}'`
  set num_patches = `awk 'NR=='$a'+1{print $1}' tmp`
  set a = `grep fd1 tmp -n | awk -F: '{print $1}'`
  set FD1 = `awk 'NR=='$a'+1{print $1}' tmp`
  rm tmp
  set master = `echo $2`
  set slave = `echo $3`
  set commandline = `echo $argv[4-$#argv]`

  if ($SAT == "ALOS") then
    echo ""
    echo "Pre-Process ALOS data - START"
    echo ""
    
    set master = ` echo $2 | awk '{ print substr($1,8,length($1)-7)}'`
    set slave = ` echo $3 | awk '{ print substr($1,8,length($1)-7)}'`
    if (! -f IMG-HH-$master || ! -f IMG-HH-$slave || ! -f LED-$master || ! -f LED-$slave) then 
      echo ""
      echo "$master $slave"
      echo "Error : Can not find input file at current directory!"   
      echo ""
      exit 1
    endif
#
# unpack the master if necessary 
#
    if (! -f IMG-HH-$master.raw || ! -f IMG-HH-$master.PRM ) then 
      echo "pre_process master image"
      ALOS_pre_process IMG-HH-$master LED-$master $argv[4-$#argv]
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
    ALOS_pre_process IMG-HH-$slave LED-$slave -fd1 $FD1 -near $NEAR -radius $RAD -npatch $num_patches -fd1 $FD1
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
      update_PRM IMG-HH-$slave.PRM input_file IMG-HH-$slave.raw
      mv IMG-HH-$slave"_"FBS.raw IMG-HH-$slave.raw
    else if  ($t == 0.5) then
      echo "Use FBS mode image as master"
      ALOS_fbd2fbs IMG-HH-$master.PRM IMG-HH-$master"_"FBS.PRM
      echo "Overwriting the old master image"
      mv IMG-HH-$master"_"FBS.PRM IMG-HH-$master.PRM
      update_PRM IMG-HH-$master.PRM input_file IMG-HH-$master.raw
      mv IMG-HH-$master"_"FBS.raw IMG-HH-$master.raw
      exit 1
    else
      echo "The range sampling rate for master and slave images are not convertable"
      exit 1
    endif 
       
    echo ""
    echo " Pre-Process ALOS data - END"
    echo ""   
       
  else if ($SAT == "ERS") then
    echo ""
    echo " Pre-Process ERS data - START"
    echo ""
  
# set 0 for master to use it's own value
    echo 
    ERS_pre_process $master $NEAR $RAD $num_patches $FD1
    set NEAR = `grep near_range $master.PRM | awk '{print $3}'`
    set RAD = `grep earth_radius $master.PRM | awk '{print $3}'`
    set FD1 = `grep fd1 $master.PRM | awk '{print $3}'`
    set num_patches = `grep num_patches $master.PRM | awk '{print $3}'`
    ERS_pre_process $slave $NEAR $RAD $num_patches $FD1
#   
#   check patch number, if different, use the smaller one
# 
    set pch1 = `grep patch $master.PRM | awk '{printf("%d ",$3)}'`
    set pch2 = `grep patch $slave.PRM | awk '{printf("%d ",$3)}'`
    echo "Different number of patches: $pch1 $pch2"
    if ($pch1 != $pch2) then
      if ($pch1 < $pch2) then
        update_PRM $slave.PRM num_patches $pch1
        echo "Number of patches is set to $pch1"
      else
        update_PRM $master.PRM num_patches $pch2
        echo "Number of patches is set to $pch2"
      endif
    endif
#
#   set the Doppler to be the average of the two
#
    grep fd1 $master.PRM | awk '{printf("%f ",$3)}' > temp
    grep fd1 $slave.PRM | awk '{printf("%f",$3)}' >> temp
    set fda = `cat temp | awk '{print( ($1 + $2)/2.)}'`
    echo " use average Doppler $fda "
    update_PRM $master.PRM fd1 $fda
    update_PRM $slave.PRM fd1 $fda
    rm -r temp
  
    echo ""
    echo " Pre-Process ERS data - END"
    echo ""
  
  else if ($SAT == "ENVI") then
    echo ""
    echo " Pre-Process ENVISAT data - START"
    echo ""
    
    ENVI_pre_process $master $NEAR $RAD $num_patches $FD1
    set NEAR = `grep near_range $master.PRM | awk '{print $3}'`
    set RAD = `grep earth_radius $master.PRM | awk '{print $3}'`
    ENVI_pre_process $slave $NEAR $RAD $num_patches $FD1
#   
#   check patch number, if different, use the smaller one
# 
    set pch1 = `grep patch $master.PRM | awk '{printf("%d ",$3)}'`
    set pch2 = `grep patch $slave.PRM | awk '{printf("%d ",$3)}'`
    echo "Different number of patches: $pch1 $pch2"
    if ($pch1 != $pch2) then
      if ($pch1 < $pch2) then
        update_PRM $slave.PRM num_patches $pch1
        echo "Number of patches is set to $pch1"
      else
        update_PRM $master.PRM num_patches $pch2
        echo "Number of patches is set to $pch2"
      endif
    endif
#
#   set the Doppler to be the average of the two
#
    grep fd1 $master.PRM | awk '{printf("%f ",$3)}' > temp
    grep fd1 $slave.PRM | awk '{printf("%f",$3)}' >> temp
    set fda = `cat temp | awk '{print( ($1 + $2)/2.)}'`
    echo " use average Doppler $fda "
    update_PRM $master.PRM fd1 $fda
    update_PRM $slave.PRM fd1 $fda
    rm -r temp
    
    echo ""
    echo " Pre-Process ENVISAT data - END"
    echo ""
  else if ($SAT == "ENVI_SLC") then
    echo ""
    echo " Pre-Process ENVISAT SLC data - START"
    echo ""
    
    ENVI_SLC_pre_process $master $RAD 
    set NEAR = `grep near_range $master.PRM | awk '{print $3}'`
    set RAD = `grep earth_radius $master.PRM | awk '{print $3}'`
    echo "ENVI_SLC_pre_process $slave $RAD "
    ENVI_SLC_pre_process $slave $RAD 
    
    echo ""
    echo " Pre-Process ENVISAT SLC data - END"
    echo ""  
  else if ($SAT == "ALOS_SLC" || $SAT == "ALOS2" || $SAT == "ALOS2_SCAN") then
    echo ""
    echo " Pre-Process ALOS SLC data - START"
    echo ""
    
    set master = ` echo $2 | awk '{ print substr($1,5,length($1)-4)}'`
    set slave =  ` echo $3 | awk '{ print substr($1,5,length($1)-4)}'`
    set master_led = ` echo $2 | awk '{ print substr($1,8,length($1)-7)}'`
    set slave_led =  ` echo $3 | awk '{ print substr($1,8,length($1)-7)}'`

#echo $master $slave $master_led $slave_led
    if (! -f IMG-$master || ! -f IMG-$slave || ! -f LED-$master_led || ! -f LED-$slave_led) then 
      echo ""
      echo "Error : Can not find input file at current directory!"   
      echo ""
      exit 1
    endif

    
    ALOS_pre_process_SLC IMG-$master LED-$master_led $commandline 
    ALOS_pre_process_SLC IMG-$slave LED-$slave_led $commandline 
    
    # make FBD FBS conversion
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
    
    echo ""
    echo " Pre-Process ALOS SLC data - END"
    echo ""
  else if ($SAT == "CSK_RAW") then
    echo ""
    echo " Pre-Process CSK Raw data - START"
    echo ""
    
    make_raw_csk $master.h5 $master
    make_raw_csk $slave.h5 $slave

#
#   calculate SC_vel and SC_height
#
    mv $master.PRM $master.PRM0
    calc_dop_orb $master.PRM0 $master.log $RAD $FD1
    cat $master.PRM0 $master.log > $master.PRM
    echo "fdd1                    = 0" >> $master.PRM
    echo "fddd1                   = 0" >> $master.PRM
#
    mv $slave.PRM $slave.PRM0
    calc_dop_orb $slave.PRM0 $slave.log $RAD $FD1
    cat $slave.PRM0 $slave.log > $slave.PRM
    echo "fdd1                    = 0" >> $slave.PRM
    echo "fddd1                   = 0" >> $slave.PRM
    rm *.log
    rm *.PRM0
#   
#   check patch number, if different, use the smaller one
# 
    set pch1 = `grep patch $master.PRM | awk '{printf("%d ",$3)}'`
    set pch2 = `grep patch $slave.PRM | awk '{printf("%d ",$3)}'`
    echo "Different number of patches: $pch1 $pch2"
    if ($pch1 != $pch2) then
      if ($pch1 < $pch2) then
        update_PRM $slave.PRM num_patches $pch1
        echo "Number of patches is set to $pch1"
      else
        update_PRM $master.PRM num_patches $pch2
        echo "Number of patches is set to $pch2"
      endif
    endif
#
#   set the Doppler to be the average of the two
#
    grep fd1 $master.PRM | awk '{printf("%f ",$3)}' > temp
    grep fd1 $slave.PRM | awk '{printf("%f",$3)}' >> temp
    set fda = `cat temp | awk '{print( ($1 + $2)/2.)}'`
    echo " use average Doppler $fda "
    update_PRM $master.PRM fd1 $fda
    update_PRM $slave.PRM fd1 $fda
    rm -r temp
    
    echo ""
    echo " Pre-Process CSK Raw data - END"
    echo ""
  else if ($SAT == "CSK_SLC" || $SAT == "TSX" || $SAT == "S1_STRIP" || $SAT == "RS2") then
    echo ""
    echo " Pre-Process CSK/TSX/RS2/S1_STRIP SLC data - START"
    echo ""
    if ($SAT == "CSK_SLC") then     
      make_slc_csk $master.h5 $master
      make_slc_csk $slave.h5 $slave
    else if ($SAT == "TSX") then 
      make_slc_tsx $master.xml $master.cos $master
      make_slc_tsx $slave.xml $slave.cos $slave      
    else if ($SAT == "RS2") then
      make_slc_rs2 $master.xml $master.tif $master
      make_slc_rs2 $slave.xml $slave.tif $slave
      mv $master.LED save-$master.LED
      extend_orbit save-$master.LED $master.LED 3
      mv $slave.LED save-$slave.LED
      extend_orbit save-$slave.LED $slave.LED 3
    else 
      make_slc_s1a $master.xml $master.tiff $master
      mv $master.LED save-$master.LED
      extend_orbit save-$master.LED $master.LED 3
      make_slc_s1a $slave.xml $slave.tiff $slave
      mv $slave.LED save-$slave.LED
      extend_orbit save-$slave.LED $slave.LED 3
    endif
#
# set the num_lines to be the min of the master and slave
#
    @ m_lines  = `grep num_lines ../raw/$master.PRM | awk '{printf("%d",int($3))}' `
    @ s_lines  = `grep num_lines ../raw/$slave.PRM | awk '{printf("%d",int($3))}' `
    if($s_lines <  $m_lines) then
      update_PRM $master.PRM num_lines $s_lines
      update_PRM $master.PRM num_valid_az $s_lines
      update_PRM $master.PRM nrows $s_lines
    else
      update_PRM $slave.PRM num_lines $m_lines
      update_PRM $slave.PRM num_valid_az $m_lines
      update_PRM $slave.PRM nrows $m_lines
    endif
#
#   calculate SC_vel and SC_height
#   set the Doppler to be zero
#
    cp $master.PRM $master.PRM0
    calc_dop_orb $master.PRM0 $master.log $RAD 0
    cat $master.PRM0 $master.log > $master.PRM
    echo "fdd1                    = 0" >> $master.PRM
    echo "fddd1                   = 0" >> $master.PRM
#
    cp $slave.PRM $slave.PRM0
    calc_dop_orb $slave.PRM0 $slave.log $RAD 0
    cat $slave.PRM0 $slave.log > $slave.PRM
    echo "fdd1                    = 0" >> $slave.PRM
    echo "fddd1                   = 0" >> $slave.PRM
    rm *.log
    rm *.PRM0
  
    echo ""
    echo " Pre-Process SLC data - END"
    echo ""
  else if ($SAT == "S1_TOPS") then
    echo ""
    echo " Pre-Process S1_TOPS data - START"
    echo ""
    if (! -f ../topo/dem.grd) then
      echo "missing file ../topo/dem.grd"
      exit 1
    endif
    ln -s ../topo/dem.grd .
    align_tops.csh $master $master.EOF $slave $slave.EOF dem.grd
    echo ""
    echo " Pre-Process S1_TOPS data - END"
    echo ""
  
  endif
     
#
