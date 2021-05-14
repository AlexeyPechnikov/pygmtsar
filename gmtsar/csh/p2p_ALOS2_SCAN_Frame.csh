#!/bin/csh -f
#
#  Xiaohua XU, Jan 10, 2019
#
# process ALOS2 ScanSAR data
# Automatically process a single frame of interferogram.
#

alias rm 'rm -f'
unset noclobber
#
  if ($#argv !=  4) then
    echo ""
    echo "Usage: p2p_ALOS2_Scan_Frame.csh Master_stem Aligned_stem config.alos2.txt parallel"
    echo ""
    echo "Example: p2p_ALOS2_Scan_Frame.csh IMG-HH-ALOS2234653650-180927-WBDR1.1__D IMG-HH-ALOS2236723650-181011-WBDR1.1__D config.alos2.txt 1" 
    echo ""
    echo "	Place the IMG-files and LED-files in the raw folder, DEM in the topo folder."
    echo "	During processing, F1 - F5 and merge folder will be generated."
    echo "  All SLCs will be upsampled to PRF 3350 and then processed."
    echo "	Final results will be placed in the merge folder, with phase"	
    echo "	corr [unwrapped phase]."
    echo "	parallel = 0-sequential  1-parallel "
    echo ""
    exit 1
  endif

# start 
#
# set processing mode seq
#
  set seq = $4
  echo "Processing 0-sequential  1-parallel [$seq] ..."

  set master = `echo $1 | awk '{print substr($1,8,length($1)-7)}'`
  set aligned = `echo $2 | awk '{print substr($1,8,length($1)-7)}'`

#if ( 6 == 9 ) then

#
# determine file names
#
  set pth = `pwd`
  foreach swath (1 2 3 4 5)
    echo "Linking files for Subswath $swath ..."
    mkdir F$swath
    cd F$swath
    mkdir raw topo
    cd topo 
    ln -s ../../topo/dem.grd .
    cd ../raw
    ln -s ../../raw/$1"-F"$swath .
    ln -s ../../raw/$2"-F"$swath .
    ln -s ../../raw/LED-$master ./LED-$master"-F"$swath
    ln -s ../../raw/LED-$aligned ./LED-$aligned"-F"$swath
    cd ..
    sed "s/.*threshold_geocode.*/threshold_geocode = 0/g" ../$3 | sed "s/.*threshold_snaphu.*/threshold_snaphu = 0/g" | sed "s/.*iono_skip_est.*/iono_skip_est = 1/g"> $3
    cd ..
  end

# 
# process data
#   
  if ($seq == 0) then
    foreach swath (1 2 3 4 5)
      cd F$swath
      sed "s/.*skip_stage.*/skip_stage = 2,3,4,5,6/g" $3 > tmp_config
      p2p_processing.csh ALOS2_SCAN $1"-F"$swath $2"-F"$swath tmp_config
      cd raw 
      samp_slc.csh $1"-F"$swath 3350 0
      samp_slc.csh $2"-F"$swath 3350 0
      cd ..
      sed "s/.*skip_stage.*/skip_stage = 1/g" $3 > tmp_config
      p2p_processing.csh ALOS2_SCAN $1"-F"$swath $2"-F"$swath tmp_config
      rm tmp_config
      cd ..
    end
  else if ($seq == 1) then
    foreach swath (1 2 3 4 5)
      cd F$swath
      sed "s/.*skip_stage.*/skip_stage = 2,3,4,5,6/g" $3 > tmp_config
      p2p_processing.csh ALOS2_SCAN $1"-F"$swath $2"-F"$swath tmp_config >&log&
      cd ..
    end
    wait 

    foreach swath (1 2 3 4 5)
      cd F$swath/raw
      samp_slc.csh $1"-F"$swath 3350 0 
      samp_slc.csh $2"-F"$swath 3350 0
      cd ../..
    end 
    wait 

    foreach swath (1 2 3 4 5)
      cd F$swath
      sed "s/.*skip_stage.*/skip_stage = 1/g" $3 > tmp_config
      p2p_processing.csh ALOS2_SCAN $1"-F"$swath $2"-F"$swath tmp_config >&log&
      cd ..
    end 
    wait 

  else
    echo "Invalid parallel mode"
    exit 1
  endif


#endif


#
# merge_unwrap_geocode
#
  mkdir merge
  cd merge
  ln -s ../topo/dem.grd .
  ln -s ../F1/intf/*/gauss* .
  set pth1 = `ls ../F1/intf/*/*PRM | awk NR==1'{print $1}' | awk -F"/" '{for (i=1;i<NF;i++) printf("%s/",$i)}'`
  set prm1m = `ls ../F1/intf/*/*PRM | awk NR==1'{print $1}' | awk -F"/" '{print $NF}'`
  set pth2 = `ls ../F2/intf/*/*PRM | awk NR==1'{print $1}' | awk -F"/" '{for (i=1;i<NF;i++) printf("%s/",$i)}'`
  set prm2m = `ls ../F2/intf/*/*PRM | awk NR==1'{print $1}' | awk -F"/" '{print $NF}'`
  set pth3 = `ls ../F3/intf/*/*PRM | awk NR==1'{print $1}' | awk -F"/" '{for (i=1;i<NF;i++) printf("%s/",$i)}'`
  set prm3m = `ls ../F3/intf/*/*PRM | awk NR==1'{print $1}' | awk -F"/" '{print $NF}'`
  set pth4 = `ls ../F4/intf/*/*PRM | awk NR==1'{print $1}' | awk -F"/" '{for (i=1;i<NF;i++) printf("%s/",$i)}'`
  set prm4m = `ls ../F4/intf/*/*PRM | awk NR==1'{print $1}' | awk -F"/" '{print $NF}'`
  set pth5 = `ls ../F5/intf/*/*PRM | awk NR==1'{print $1}' | awk -F"/" '{for (i=1;i<NF;i++) printf("%s/",$i)}'`
  set prm5m = `ls ../F5/intf/*/*PRM | awk NR==1'{print $1}' | awk -F"/" '{print $NF}'`
  
  echo $pth1$prm1m":"$pth1"phasefilt.grd" > phaselist
  echo $pth2$prm2m":"$pth2"phasefilt.grd" >> phaselist
  echo $pth3$prm3m":"$pth3"phasefilt.grd" >> phaselist
  echo $pth4$prm4m":"$pth4"phasefilt.grd" >> phaselist
  echo $pth5$prm5m":"$pth5"phasefilt.grd" >> phaselist

  head -3 phaselist > first.txt
  merge_swath first.txt first.grd first
  echo "first.PRM:first.grd" > second.txt
  tail -2 phaselist >> second.txt
  merge_swath second.txt second.grd second
  mv second.PRM $prm1m
  mv second.grd phasefilt.grd
  rm first* second*


  echo $pth1$prm1m":"$pth1"corr.grd" > corrlist
  echo $pth2$prm2m":"$pth2"corr.grd" >> corrlist
  echo $pth3$prm3m":"$pth3"corr.grd" >> corrlist
  echo $pth4$prm4m":"$pth4"corr.grd" >> corrlist
  echo $pth5$prm5m":"$pth5"corr.grd" >> corrlist

  head -3 corrlist > first.txt
  merge_swath first.txt first.grd first
  echo "first.PRM:first.grd" > second.txt
  tail -2 corrlist >> second.txt
  merge_swath second.txt second.grd second
  mv second.grd corr.grd
  rm first* second*


  echo $pth1$prm1m":"$pth1"mask.grd" > masklist
  echo $pth2$prm2m":"$pth2"mask.grd" >> masklist
  echo $pth3$prm3m":"$pth3"mask.grd" >> masklist
  echo $pth4$prm4m":"$pth4"mask.grd" >> masklist
  echo $pth5$prm5m":"$pth5"mask.grd" >> masklist

  head -3 masklist > first.txt
  merge_swath first.txt first.grd first
  echo "first.PRM:first.grd" > second.txt
  tail -2 masklist >> second.txt
  merge_swath second.txt second.grd second
  mv second.grd mask.grd
  rm first* second*

  if (! -f trans.dat) then
    if (! -f dem.grd) then
      echo "ERROR: missing dem.grd ... (link from the topo folder)"
      exit 1
    endif
    set led = `grep led_file $prm1m | awk '{print $3}'`
    cp $pth1$led .
    gmt grd2xyz --FORMAT_FLOAT_OUT=%lf dem.grd -s | SAT_llt2rat $prm1m 1 -bod > trans.dat
  endif

    # Read in parameters
  set threshold_snaphu = `grep threshold_snaphu $2 | awk '{print $3}'`
  set threshold_geocode = `grep threshold_geocode $2 | awk '{print $3}'`
  set region_cut = `grep region_cut $2 | awk '{print $3}'`
  set mask_water = `grep mask_water $2 | awk '{print $3}'`
  set defomax = `grep defomax $2 | awk '{print $3}'`
  set near_interp = `grep near_interp $2 | awk '{print $3}'`

  # Unwrapping
  if ($region_cut == "") then
    set region_cut = `gmt grdinfo phasefilt.grd -I- | cut -c3-20`
  endif
  if ($threshold_snaphu != 0 ) then
    if ($mask_water == 1) then
      if (! -f landmask_ra.grd) then
        landmask.csh $region_cut
      endif
    endif

    echo ""
    echo "SNAPHU.CSH - START"
    echo "threshold_snaphu: $threshold_snaphu"
    if ($near_interp == 1) then
      snaphu_interp.csh $threshold_snaphu $defomax $region_cut
    else
      snaphu.csh $threshold_snaphu $defomax $region_cut
    endif
    echo "SNAPHU.CSH - END"
  else
    echo ""
    echo "SKIP UNWRAP PHASE"
  endif

  if ($threshold_geocode != 0) then
    echo ""
    echo "GEOCODE-START"
    proj_ra2ll.csh trans.dat phasefilt.grd phasefilt_ll.grd
    proj_ra2ll.csh trans.dat corr.grd corr_ll.grd
    gmt makecpt -Crainbow -T-3.15/3.15/0.05 -Z > phase.cpt
    set BT = `gmt grdinfo -C corr.grd | awk '{print $7}'`
    gmt makecpt -Cgray -T0/$BT/0.05 -Z > corr.cpt
    grd2kml.csh phasefilt_ll phase.cpt
    grd2kml.csh corr_ll corr.cpt

    if (-f unwrap.grd) then
      gmt grdmath unwrap.grd mask.grd MUL = unwrap_mask.grd
      proj_ra2ll.csh trans.dat unwrap.grd unwrap_ll.grd
      proj_ra2ll.csh trans.dat unwrap_mask.grd unwrap_mask_ll.grd
      set BT = `gmt grdinfo -C unwrap.grd | awk '{print $7}'`
      set BL = `gmt grdinfo -C unwrap.grd | awk '{print $6}'`
      gmt makecpt -T$BL/$BT/0.5 -Z > unwrap.cpt
      grd2kml.csh unwrap_mask_ll unwrap.cpt
      grd2kml.csh unwrap_ll unwrap.cpt
    endif

    echo "GEOCODE END"
  endif 













