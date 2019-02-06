#!/bin/csh -f
#       $Id$
#
#  Xiaohua Xu, Jan, 2019
#
#  Automatically generate Multi Aperture Interferogram
#
  if ($#argv != 3 && $#argv != 4) then
    echo ""
    echo "Usage: MAI_processing.csh SAT master_image slave_image [configuration_file] "
    echo ""
    echo "Example: MAI_processing.csh ALOS IMG-HH-ALPSRP055750660-H1.0__A IMG-HH-ALPSRP049040660-H1.0__A [config.alos.txt]"
    echo ""
    exit 1
  endif

  if ($#argv == 4) then 
    if(! -f $4 ) then
      echo " no configure file: "$4
      echo " Leave it blank to generate config file with default values."
      exit 1
    endif
  endif 
 
  set SAT = `echo $1`
  if ($#argv == 4) then
    set conf = `echo $4`
  else
    pop_config.csh $SAT > config.$SAT.txt
    set conf = `echo "config.$SAT.txt"`
  endif
  set master = ` echo $2 `
  set slave =  ` echo $3 `
  set region_cut = `grep region_cut $conf | awk '{print $3}'`
  set filter = `grep filter_wavelength $conf | awk '{print $3}'`
  set dec = `grep dec_factor $conf | awk '{print $3}'`
  set range_dec = `grep range_dec $conf | awk '{print $3}'`
  set azimuth_dec = `grep azimuth_dec $conf | awk '{print $3}'`
  set near_interp = `grep near_interp $conf | awk '{print $3}'`
  set threshold_snaphu = `grep threshold_snaphu $conf | awk '{print $3}'`
  set threshold_geocode = `grep threshold_geocode $conf | awk '{print $3}'`
  set mask_water = `grep mask_water $conf | awk '{print $3}'`
  set switch_land = `grep switch_land $conf | awk '{print $3}'`
  set defomax = `grep defomax $conf | awk '{print $3}'`
  set near_interp = `grep near_interp $conf | awk '{print $3}'`

  sed "s/.*skip_stage.*/skip_stage = 3,4,5,6/g" $conf > tmp_config
  p2p_processing.csh $SAT $master $slave tmp_config

  cd raw
  split_aperture $master.PRM > MAI_m.rec
  split_aperture $slave.PRM > MAI_s.rec
  cd ..

  mkdir SLC_F
  cd SLC_F
  ln -s ../raw/$master.LED .
  cp ../raw/$master.PRM .
  ln -s ../raw/$master.SLCF ./$master.SLC
  ln -s ../raw/$slave.LED .
  cp ../raw/$slave.PRM .
  ln -s ../raw/$slave.SLCF ./$slave.SLC
  cp $slave.PRM $slave.PRM0
  if ($SAT == "ALOS2_SCAN") then
    ln -s ../SLC/freq_alos2.dat
    fitoffset.csh  2 3 freq_alos2.dat 10 >> $slave.PRM
  else if ($SAT == "ERS" || $SAT == "ENVI" || $SAT == "ALOS" || $SAT == "CSK_RAW") then
    ln -s ../SLC/freq_xcorr.dat .
    fitoffset.csh 3 3 freq_xcorr.dat 18 >> $slave.PRM
  else
    ln -s ../SLC/freq_xcorr.dat .
    fitoffset.csh 2 2 freq_xcorr.dat 18 >> $slave.PRM
  endif
  resamp $master.PRM $slave.PRM $slave.PRMresamp $slave.SLCresamp 4
  rm $slave.SLC
  mv $slave.SLCresamp $slave.SLC
  cp $slave.PRMresamp $slave.PRM
  if ($region_cut != "") then
    echo "Cutting SLC image to $region_cut"
    cut_slc $master.PRM junk1 $region_cut
    cut_slc $slave.PRM junk2 $region_cut
    mv junk1.PRM $master.PRM
    mv junk2.PRM $slave.PRM
    mv junk1.SLC $master.SLC
    mv junk2.SLC $slave.SLC
  endif
  intf.csh $master.PRM $slave.PRM
  cd ..

  mkdir SLC_B
  cd SLC_B
  ln -s ../raw/$master.LED .
  cp ../raw/$master.PRM .
  ln -s ../raw/$master.SLCB ./$master.SLC
  ln -s ../raw/$slave.LED .
  cp ../raw/$slave.PRM .
  ln -s ../raw/$slave.SLCB ./$slave.SLC
  cp $slave.PRM $slave.PRM0
  if ($SAT == "ALOS2_SCAN") then
    ln -s ../SLC/freq_alos2.dat
    fitoffset.csh  2 3 freq_alos2.dat 10 >> $slave.PRM
  else if ($SAT == "ERS" || $SAT == "ENVI" || $SAT == "ALOS" || $SAT == "CSK_RAW") then
    ln -s ../SLC/freq_xcorr.dat .
    fitoffset.csh 3 3 freq_xcorr.dat 18 >> $slave.PRM
  else
    ln -s ../SLC/freq_xcorr.dat .
    fitoffset.csh 2 2 freq_xcorr.dat 18 >> $slave.PRM
  endif
  resamp $master.PRM $slave.PRM $slave.PRMresamp $slave.SLCresamp 4
  rm $slave.SLC
  mv $slave.SLCresamp $slave.SLC
  cp $slave.PRMresamp $slave.PRM
  if ($region_cut != "") then
    echo "Cutting SLC image to $region_cut"
    cut_slc $master.PRM junk1 $region_cut
    cut_slc $slave.PRM junk2 $region_cut
    mv junk1.PRM $master.PRM
    mv junk2.PRM $slave.PRM
    mv junk1.SLC $master.SLC
    mv junk2.SLC $slave.SLC
  endif
  intf.csh $master.PRM $slave.PRM
  cd .. 

  cd raw

  set prf = `grep PRF $master.PRM | awk '{print $3}'`
  set height = `grep SC_height $master.PRM | head -1 | awk '{print $3}'`
  set radius = `grep earth_radius $master.PRM | head -1 | awk '{print $3}'`
  set SC_vel = `grep SC_vel $master.PRM | head -1 | awk '{print $3}'`
  set spec_sep = `grep frequency_separation MAI_m.rec | awk '{print $3}'`
  cd ..

  set pix_size = `echo $prf $height $radius $SC_vel | awk '{printf("%.6f",$4/($2+$3)*$3/$1)}'`
  mkdir MAI_intf
  cd MAI_intf
  gmt grdmath ../SLC_F/real.grd=bf ../SLC_B/real.grd=bf MUL ../SLC_F/imag.grd=bf ../SLC_B/imag.grd=bf MUL ADD 1e7 MUL = real.grd=bf
  gmt grdmath ../SLC_F/imag.grd=bf ../SLC_B/real.grd=bf MUL ../SLC_F/real.grd=bf ../SLC_B/imag.grd=bf MUL -1 MUL ADD 1e7 MUL = imag.grd=bf
  cp ../SLC/$master.PRM .
  cp ../SLC/$slave.PRM .
  ln -s ../SLC/$master.LED .
  ln -s ../SLC/$slave.LED .
  ln -s ../SLC/$master.SLC .
  ln -s ../SLC/$slave.SLC .
  filter.csh $master.PRM $slave.PRM $filter $dec $range_dec $azimuth_dec

  gmt grdmath phase.grd 2 PI MUL DIV $spec_sep DIV $prf MUL $pix_size MUL = MAI_intf.grd
  gmt grdmath phasefilt.grd 2 PI MUL DIV $spec_sep DIV $prf MUL $pix_size MUL = MAI_intf_filt.grd
  
  #gmt grdmath ../intf_f/*/phase.grd ../intf_b/*/phase.grd SUB PI ADD 2 PI MUL MOD PI SUB 2 PI MUL DIV $spec_sep DIV $prf MUL $pix_size MUL = MAI_intf.grd
  if ($threshold_geocode != 0) then
    ln -s ../topo/trans.dat .
    if ($mask_water != 0 || $switch_land != 0) then
      set region = `gmt grdinfo MAI_intf.grd -I- | awk -F'R' '{print $2}'`
      cd ../topo
      landmask.csh $region
      cd ../MAI_intf
      ln -s ../topo/landmask_ra.grd .
      gmt grdsample landmask_ra.grd -RMAI_intf.grd -Glandmask_ra_patch.grd
      gmt grdmath MAI_intf.grd landmask_ra_patch.grd MUL = tmp.grd
      mv tmp.grd MAI_intf.grd
    endif
    proj_ra2ll.csh trans.dat MAI_intf.grd MAI_intf_ll.grd
    proj_ra2ll.csh trans.dat MAI_intf_filt.grd MAI_intf_filt_ll.grd
    gmt makecpt -Cjet -T-$pix_size/$pix_size/0.1 -D -Z > MAI.cpt
    grd2kml.csh MAI_intf_ll MAI.cpt
    grd2kml.csh MAI_intf_filt_ll MAI.cpt
  endif
  cd ..






  

