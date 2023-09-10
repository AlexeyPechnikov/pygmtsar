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
  if (!($#argv == 3 || $#argv == 5 || $#argv == 7 || $#argv == 9 || $#argv == 11 || $#argv == 13 )) then 
    echo ""
    echo "Usage: pre_proc.csh SAT master_stem aligned_stem [-near near_range] [-radius RE] [-npatch num_patches] [-fd1 DOPP] [-ESD mode] [-skip_master 1]"
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
  echo "0" $argv[4-$#argv] | awk '{for (i=1;i<=20;i++) print $i}' > tmp
  set a = `grep near tmp -n | awk -F: '{print $1}'`
  set NEAR = `awk 'NR=='$a'+1{print $1}' tmp`
  set a = `grep near tmp -n | awk -F: '{print $1}'`
  set RAD = `awk 'NR=='$a'+1{print $1}' tmp`
  set a = `grep npatch tmp -n | awk -F: '{print $1}'`
  set num_patches = `awk 'NR=='$a'+1{print $1}' tmp`
  set a = `grep fd1 tmp -n | awk -F: '{print $1}'`
  set FD1 = `awk 'NR=='$a'+1{print $1}' tmp`
  set a = `grep SLC_factor tmp -n | awk -F: '{print $1}'`
  set SLC_factor = `awk 'NR=='$a'+1{print $1}' tmp`
  set a = `grep ESD tmp -n | awk -F: '{print $1}'`
  set ESD_mode = `awk 'NR=='$a'+1{print $1}' tmp`
  set a = `grep skip_master tmp -n | awk -F: '{print $1}'`
  set skip_master = `awk 'NR=='$a'+1{print $1}' tmp`

  set commandline = ""
  if (!($NEAR == 0)) then 
    set commandline = "$commandline -near $NEAR"
  endif
  if (!($RAD == 0)) then 
    set commandline = "$commandline -radius $RAD"
  endif
  if (!($num_patches == 0)) then  
    set commandline = "$commandline -npatch $num_patches"
  endif
  if (!($FD1 == 0)) then  
    set commandline = "$commandline -fd1 $FD1"
  endif
  if (!($SLC_factor == 0)) then  
    set commandline = "$commandline -SLC_factor $SLC_factor"
  endif 

echo $commandline

  if ($skip_master != 0) then
    echo ""
    echo "skip_master set to $skip_master "
    echo ""
  endif

  rm tmp
  set master = `echo $2`
  set aligned = `echo $3`

  if ($SAT == "ALOS") then
    echo ""
    echo "Pre-Process ALOS data - START"
    echo ""
    
    set master = ` echo $2 | awk '{ print substr($1,8,length($1)-7)}'`
    set aligned = ` echo $3 | awk '{ print substr($1,8,length($1)-7)}'`
    if (! -f IMG-HH-$master || ! -f IMG-HH-$aligned || ! -f LED-$master || ! -f LED-$aligned) then 
      echo ""
      echo "$master $aligned"
      echo "Error : Can not find input file at current directory!"   
      echo ""
      exit 1
    endif
#
# unpack the master if necessary 
#
    if (! -f IMG-HH-$master.raw || ! -f IMG-HH-$master.PRM ) then 
      echo "pre_process master image"
      if ($skip_master == 0 || $skip_master == 2) then
        ALOS_pre_process IMG-HH-$master LED-$master $commandline
      endif
    endif

    set NEAR = `grep near_range IMG-HH-$master.PRM | awk '{print $3}'`
    set RAD = `grep earth_radius IMG-HH-$master.PRM | awk '{print $3}'`
    set FD1 = `grep fd1 IMG-HH-$master.PRM | awk '{print $3}'`
    set num_patches = `grep num_patches IMG-HH-$master.PRM | awk '{print $3}'`
#
# unpack the aligned image using the same earth radius and near range as the master image
#
    echo "pre_process aligned image"
    if ($skip_master == 0 || $skip_master == 1) then
      ALOS_pre_process IMG-HH-$aligned LED-$aligned -fd1 $FD1 -near $NEAR -radius $RAD -npatch $num_patches -fd1 $FD1
    endif
#
# check the range sampling rate of the aligned images and do conversion if necessary
#
    if ($skip_master == 0 || $skip_master == 1) then
      set rng_samp_rate_m = `grep rng_samp_rate IMG-HH-$master.PRM | awk 'NR == 1 {printf("%d", $3)}'`
      set rng_samp_rate_s = `grep rng_samp_rate IMG-HH-$aligned.PRM | awk 'NR == 1 {printf("%d", $3)}'`
      set t = `echo $rng_samp_rate_m $rng_samp_rate_s | awk '{printf("%1.1f\n", $1/$2)}'`

      if ($t == 1.0) then
        echo "The range sampling rate for master and aligned images are: "$rng_samp_rate_m
      else if ($t == 2.0) then
        echo "Convert the aligned image from FBD to FBS mode"
        ALOS_fbd2fbs IMG-HH-$aligned.PRM IMG-HH-$aligned"_"FBS.PRM
        echo "Overwriting the old aligned image"
        mv IMG-HH-$aligned"_"FBS.PRM IMG-HH-$aligned.PRM
        update_PRM IMG-HH-$aligned.PRM input_file IMG-HH-$aligned.raw
        mv IMG-HH-$aligned"_"FBS.raw IMG-HH-$aligned.raw
        echo "IMG-HH-$aligned is converted to FBS mode" > ALOS_fbd2fbs_log_$aligned
      else if  ($t == 0.5) then
        echo "Use FBS mode image as master"
        ALOS_fbd2fbs IMG-HH-$master.PRM IMG-HH-$master"_"FBS.PRM
        echo "Overwriting the old master image"
        mv IMG-HH-$master"_"FBS.PRM IMG-HH-$master.PRM
        update_PRM IMG-HH-$master.PRM input_file IMG-HH-$master.raw
        mv IMG-HH-$master"_"FBS.raw IMG-HH-$master.raw
        echo "IMG-HH-$master is converted to FBS mode" > ALOS_fbd2fbs_log_$master
        exit 1
      else
        echo "The range sampling rate for master and aligned images are not convertable"
        exit 1
      endif 
    endif
       
    echo ""
    echo " Pre-Process ALOS data - END"
    echo ""   
       
  else if ($SAT == "ERS") then
    echo ""
    echo " Pre-Process ERS data - START"
    echo ""
  
# set 0 for master to use it's own value
    if ($skip_master == 0 || $skip_master == 2) then
      ERS_pre_process $master $NEAR $RAD $num_patches $FD1
    endif
    set NEAR = `grep near_range $master.PRM | awk '{print $3}'`
    set RAD = `grep earth_radius $master.PRM | awk '{print $3}'`
    set FD1 = `grep fd1 $master.PRM | awk '{print $3}'`
    set num_patches = `grep num_patches $master.PRM | awk '{print $3}'`
    if ($skip_master == 0 || $skip_master == 1) then
      ERS_pre_process $aligned $NEAR $RAD $num_patches $FD1
    endif
#   
#   check patch number, if different, use the smaller one
# 
    if ($skip_master == 0) then
      set pch1 = `grep patch $master.PRM | awk '{printf("%d ",$3)}'`
      set pch2 = `grep patch $aligned.PRM | awk '{printf("%d ",$3)}'`
      echo "Different number of patches: $pch1 $pch2"
      if ($pch1 != $pch2) then
        if ($pch1 < $pch2) then
          update_PRM $aligned.PRM num_patches $pch1
          echo "Number of patches is set to $pch1"
        else
          update_PRM $master.PRM num_patches $pch2
          echo "Number of patches is set to $pch2"
        endif
      endif
    else
      if ($skip_master == 1) then
        set pch1 = `grep patch $master.PRM | awk '{printf("%d ",$3)}'`
        update_PRM $aligned.PRM num_patches $pch1
        echo "Number of patches is set to $pch1"
      endif
    endif
#
#   set the Doppler to be the average of the two
#
    if ($skip_master == 0) then
      grep fd1 $master.PRM | awk '{printf("%f ",$3)}' > temp
      grep fd1 $aligned.PRM | awk '{printf("%f",$3)}' >> temp
      set fda = `cat temp | awk '{print( ($1 + $2)/2.)}'`
      echo " use average Doppler $fda "
      update_PRM $master.PRM fd1 $fda
      update_PRM $aligned.PRM fd1 $fda
      rm -r temp
    else
      if ($skip_master == 1) then
        set fda = `grep fd1 $master.PRM | awk '{print $3}'`
        update_PRM $aligned.PRM fd1 $fda
      endif
    else
      
    endif
  
    echo ""
    echo " Pre-Process ERS data - END"
    echo ""
  
  else if ($SAT == "ENVI") then
    echo ""
    echo " Pre-Process ENVISAT data - START"
    echo ""
    if ($skip_master == 0 || $skip_master == 2) then
      ENVI_pre_process $master $NEAR $RAD $num_patches $FD1
    endif
    set NEAR = `grep near_range $master.PRM | awk '{print $3}'`
    set RAD = `grep earth_radius $master.PRM | awk '{print $3}'`
    if ($skip_master == 0 || $skip_master == 1) then
      ENVI_pre_process $aligned $NEAR $RAD $num_patches $FD1
    endif
#   
#   check patch number, if different, use the smaller one
# 
    if ($skip_master == 0) then
      set pch1 = `grep patch $master.PRM | awk '{printf("%d ",$3)}'`
      set pch2 = `grep patch $aligned.PRM | awk '{printf("%d ",$3)}'`
      echo "Different number of patches: $pch1 $pch2"
      if ($pch1 != $pch2) then
        if ($pch1 < $pch2) then
          update_PRM $aligned.PRM num_patches $pch1
          echo "Number of patches is set to $pch1"
        else
          update_PRM $master.PRM num_patches $pch2
          echo "Number of patches is set to $pch2"
        endif
      endif
    else
      if ($skip_master == 1) then
        set pch1 = `grep patch $master.PRM | awk '{printf("%d ",$3)}'`
        update_PRM $aligned.PRM num_patches $pch1
        echo "Number of patches is set to $pch1"
      endif
    endif
#
#   set the Doppler to be the average of the two
#
    if ($skip_master == 0) then 
      grep fd1 $master.PRM | awk '{printf("%f ",$3)}' > temp
      grep fd1 $aligned.PRM | awk '{printf("%f",$3)}' >> temp
      set fda = `cat temp | awk '{print( ($1 + $2)/2.)}'`
      echo " use average Doppler $fda "
      update_PRM $master.PRM fd1 $fda
      update_PRM $aligned.PRM fd1 $fda
      rm -r temp
    else
      if ($skip_master == 1) then
        set fda = `grep fd1 $master.PRM | awk '{print $3}'`
        update_PRM $aligned.PRM fd1 $fda
      endif
    endif
    
    echo ""
    echo " Pre-Process ENVISAT data - END"
    echo ""
  else if ($SAT == "ENVI_SLC") then
    echo ""
    echo " Pre-Process ENVISAT SLC data - START"
    echo ""
    
    if ($skip_master == 0 || $skip_master == 2) then 
      ENVI_SLC_pre_process $master $RAD 
    endif
    set NEAR = `grep near_range $master.PRM | awk '{print $3}'`
    set RAD = `grep earth_radius $master.PRM | awk '{print $3}'`
    echo "ENVI_SLC_pre_process $aligned $RAD "
    if ($skip_master == 0 || $skip_master == 1) then
      ENVI_SLC_pre_process $aligned $RAD 
    endif
    
    echo ""
    echo " Pre-Process ENVISAT SLC data - END"
    echo ""  
  else if ($SAT == "ALOS_SLC" || $SAT == "ALOS2" || $SAT == "ALOS2_SCAN") then
    echo ""
    echo " Pre-Process ALOS SLC data - START"
    echo ""
    
    set master = ` echo $2 | awk '{ print substr($1,5,length($1)-4)}'`
    set aligned =  ` echo $3 | awk '{ print substr($1,5,length($1)-4)}'`
    set master_led = ` echo $2 | awk '{ print substr($1,8,length($1)-7)}'`
    set aligned_led =  ` echo $3 | awk '{ print substr($1,8,length($1)-7)}'`

#echo $master $aligned $master_led $aligned_led
    if (! -f IMG-$master || ! -f IMG-$aligned || ! -f LED-$master_led || ! -f LED-$aligned_led) then 
      echo ""
      echo "Error : Can not find input file at current directory!"   
      echo ""
      exit 1
    endif

    if ($SAT == "ALOS_SLC") then 
      if ($skip_master == 0 || $skip_master == 2) then
        ALOS_pre_process_SLC IMG-$master LED-$master_led $commandline -ALOS1
      endif
      if ($skip_master == 0 || $skip_master == 1) then
        ALOS_pre_process_SLC IMG-$aligned LED-$aligned_led $commandline -ALOS1
      endif
    else
      if ($skip_master == 0 || $skip_master == 2) then
        ALOS_pre_process_SLC IMG-$master LED-$master_led $commandline 
      endif
      if ($skip_master == 0 || $skip_master == 1) then
        ALOS_pre_process_SLC IMG-$aligned LED-$aligned_led $commandline 
      endif
    endif
    
    # make FBD FBS conversion
    if ($skip_master == 0 || $skip_master == 1) then
      set rng_samp_rate_m = `grep rng_samp_rate IMG-$master.PRM | awk 'NR == 1 {printf("%d", $3)}'`
      set rng_samp_rate_s = `grep rng_samp_rate IMG-$aligned.PRM | awk 'NR == 1 {printf("%d", $3)}'`
      set t = `echo $rng_samp_rate_m $rng_samp_rate_s | awk '{printf("%1.1f\n", $1/$2)}'`
      if ($t == 1.0) then
        echo "The range sampling rate for master and aligned images are: "$rng_samp_rate_m
      else if ($t == 2.0) then
        echo "Convert the aligned image from FBD to FBS mode"
	    ALOS_fbd2fbs_SLC IMG-$aligned.PRM IMG-$aligned"_"FBS.PRM
        echo "Overwriting the old aligned image"
        mv IMG-$aligned"_"FBS.PRM IMG-$aligned.PRM
	    update_PRM IMG-$aligned.PRM input_file IMG-$aligned.SLC
        mv IMG-$aligned"_"FBS.SLC IMG-$aligned.SLC
        echo "IMG-HH-$aligned is converted to FBS mode" > ALOS2_fbd2fbs_log_"$aligned"
      else if  ($t == 0.5) then
	    echo "Convert the master image from FBD to FBS mode"
	    ALOS_fbd2fbs_SLC IMG-$master.PRM IMG-$master"_"FBS.PRM
        echo "Overwriting the old master image"
        mv IMG-$master"_"FBS.PRM IMG-$master.PRM
	    update_PRM IMG-$master.PRM input_file IMG-$master.SLC
        mv IMG-$master"_"FBS.SLC IMG-$master.SLC
        echo "IMG-HH-$master is converted to FBS mode" > ALOS2_fbd2fbs_log_"$master"
      else
	    echo "The range sampling rate for master and aligned images are not convertable"
        exit 1
      endif
    endif
    
    echo ""
    echo " Pre-Process ALOS SLC data - END"
    echo ""
  else if ($SAT == "CSK_RAW") then
    echo ""
    echo " Pre-Process CSK Raw data - START"
    echo ""
    if ($skip_master == 0 || $skip_master == 2) then
      make_raw_csk $master.h5 $master
      # calculate SC_vel and SC_height
      mv $master.PRM $master.PRM0
      calc_dop_orb $master.PRM0 $master.log $RAD $FD1
      cat $master.PRM0 $master.log > $master.PRM
      echo "fdd1                    = 0" >> $master.PRM
      echo "fddd1                   = 0" >> $master.PRM
    endif
    if ($skip_master == 0 || $skip_master == 1) then
      make_raw_csk $aligned.h5 $aligned
      # calculate SC_vel and SC_height
      #mv $master.PRM $master.PRM0
      mv $aligned.PRM $aligned.PRM0
      calc_dop_orb $aligned.PRM0 $aligned.log $RAD $FD1
      cat $aligned.PRM0 $aligned.log > $aligned.PRM
      echo "fdd1                    = 0" >> $aligned.PRM
      echo "fddd1                   = 0" >> $aligned.PRM
      rm *.log
      rm *.PRM0
    endif
#
#   
#   check patch number, if different, use the smaller one
# 
    if ($skip_master == 0) then
      set pch1 = `grep patch $master.PRM | awk '{printf("%d ",$3)}'`
      set pch2 = `grep patch $aligned.PRM | awk '{printf("%d ",$3)}'`
      echo "Different number of patches: $pch1 $pch2"
      if ($pch1 != $pch2) then
        if ($pch1 < $pch2) then
          update_PRM $aligned.PRM num_patches $pch1
          echo "Number of patches is set to $pch1"
        else
          update_PRM $master.PRM num_patches $pch2
          echo "Number of patches is set to $pch2"
        endif
      endif
    else
      if ($skip_master == 1) then
        set pch1 = `grep patch $master.PRM | awk '{printf("%d ",$3)}'`
        update_PRM $aligned.PRM num_patches $pch1
        echo "Number of patches is set to $pch1"
      endif
    endif
#
#   set the Doppler to be the average of the two
#
    if ($skip_master == 0) then
      grep fd1 $master.PRM | awk '{printf("%f ",$3)}' > temp
      grep fd1 $aligned.PRM | awk '{printf("%f",$3)}' >> temp
      set fda = `cat temp | awk '{print( ($1 + $2)/2.)}'`
      echo " use average Doppler $fda "
      update_PRM $master.PRM fd1 $fda
      update_PRM $aligned.PRM fd1 $fda
      rm -r temp
    else
      if ($skip_master == 1) then
        set fda = `grep fd1 $master.PRM | awk '{print $3}'`
        update_PRM $aligned.PRM fd1 $fda
      endif
    endif
    
    echo ""
    echo " Pre-Process CSK Raw data - END"
    echo ""
  else if ($SAT == "CSK_SLC" || $SAT == "TSX" || $SAT == "S1_STRIP" || $SAT == "RS2" || $SAT == "GF3" || $SAT == "LT1") then
    echo ""
    echo " Pre-Process CSK/TSX/RS2/S1_STRIP/GF3/LT1 SLC data - START"
    echo ""
    if ($SAT == "CSK_SLC") then     
      if ($skip_master == 0 || $skip_master == 2) then
        make_slc_csk $master.h5 $master
      endif
      if ($skip_master == 0 || $skip_master == 1) then
        make_slc_csk $aligned.h5 $aligned
      endif
    else if ($SAT == "TSX") then 
      if ($skip_master == 0 || $skip_master == 2) then
        make_slc_tsx $master.xml $master.cos $master
      endif
      if ($skip_master == 0 || $skip_master == 1) then
        make_slc_tsx $aligned.xml $aligned.cos $aligned      
      endif
    else if ($SAT == "RS2") then
      if ($skip_master == 0 || $skip_master == 2) then
        make_slc_rs2 $master.xml $master.tif $master
        mv $master.LED save-$master.LED
        extend_orbit save-$master.LED $master.LED 3
      endif
      if ($skip_master == 0 || $skip_master == 1) then
        make_slc_rs2 $aligned.xml $aligned.tif $aligned
        mv $aligned.LED save-$aligned.LED
        extend_orbit save-$aligned.LED $aligned.LED 3
      endif
    else if ($SAT == "GF3") then
      if ($skip_master == 0 || $skip_master == 2) then
        make_slc_gf3 $master.xml $master.tiff $master
      endif
      if ($skip_master == 0 || $skip_master == 1) then
        make_slc_gf3 $aligned.xml $aligned.tiff $aligned
      endif
    else if ($SAT == "LT1") then
      if ($skip_master == 0 || $skip_master == 2) then
        make_slc_lt1 $master.xml $master.tiff $master
      endif
      if ($skip_master == 0 || $skip_master == 1) then
        make_slc_lt1 $aligned.xml $aligned.tiff $aligned
      endif
    else 
      if ($skip_master == 0 || $skip_master == 2) then
        make_slc_s1a $master.xml $master.tiff $master
        mv $master.LED save-$master.LED
        extend_orbit save-$master.LED $master.LED 3
      endif
      if ($skip_master == 0 || $skip_master == 1) then
        make_slc_s1a $aligned.xml $aligned.tiff $aligned
        mv $aligned.LED save-$aligned.LED
        extend_orbit save-$aligned.LED $aligned.LED 3
      endif
    endif
#
# set the num_lines to be the min of the master and aligned
#
    if ($skip_master == 0) then 
      @ m_lines  = `grep num_lines ../raw/$master.PRM | awk '{printf("%d",int($3))}' `
      @ s_lines  = `grep num_lines ../raw/$aligned.PRM | awk '{printf("%d",int($3))}' `
      if($s_lines <  $m_lines) then
        update_PRM $master.PRM num_lines $s_lines
        update_PRM $master.PRM num_valid_az $s_lines
        update_PRM $master.PRM nrows $s_lines
      else
        update_PRM $aligned.PRM num_lines $m_lines
        update_PRM $aligned.PRM num_valid_az $m_lines
        update_PRM $aligned.PRM nrows $m_lines
      endif
    else
      if ($skip_master == 1) then
        @ m_lines  = `grep num_lines ../raw/$master.PRM | awk '{printf("%d",int($3))}' `
        update_PRM $aligned.PRM num_lines $m_lines
        update_PRM $aligned.PRM num_valid_az $m_lines
        update_PRM $aligned.PRM nrows $m_lines 
      endif
    endif

#
#   calculate SC_vel and SC_height
#   set the Doppler to be zero
#
    if ($skip_master == 0 || $skip_master == 2) then
      cp $master.PRM $master.PRM0
      calc_dop_orb $master.PRM0 $master.log $RAD 0
      cat $master.PRM0 $master.log > $master.PRM
      echo "fdd1                    = 0" >> $master.PRM
      echo "fddd1                   = 0" >> $master.PRM
    endif
#
    if ($skip_master == 0 || $skip_master == 1) then
      cp $aligned.PRM $aligned.PRM0
      calc_dop_orb $aligned.PRM0 $aligned.log $RAD 0
      cat $aligned.PRM0 $aligned.log > $aligned.PRM
      echo "fdd1                    = 0" >> $aligned.PRM
      echo "fddd1                   = 0" >> $aligned.PRM
      rm *.log
      rm *.PRM0
    endif
  
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

    if ($ESD_mode == 0) then
      if ($skip_master == 0) then
        align_tops.csh $master $master.EOF $aligned $aligned.EOF dem.grd
      else if ($skip_master == 1) then
        align_tops.csh $master 0  $aligned $aligned.EOF dem.grd
      else if ($skip_master == 2) then
        align_tops.csh $master $master.EOF $aligned 0 dem.grd
      else
        echo "[ERROR]: Wrong skip_master parameter used for TOPS data, please check"
        exit 1
      endif
    else
      if ($skip_master == 0) then
        align_tops_esd.csh $master $master.EOF $aligned $aligned.EOF dem.grd $ESD_mode
      else if ($skip_master == 1) then
        align_tops_esd.csh $master 0 $aligned $aligned.EOF dem.grd $ESD_mode
      else if ($skip_master == 2) then
        align_tops_esd.csh $master $master.EOF $aligned 0 dem.grd $ESD_mode
      else
        echo "[ERROR]: Wrong skip_master parameter used for TOPS data, please check"
        exit 1
      endif
    endif

    echo ""
    echo " Pre-Process S1_TOPS data - END"
    echo ""
  
  endif
     
#
