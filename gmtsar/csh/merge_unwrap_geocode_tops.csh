#!/bin/csh -f
#       $Id$
#
#
#    Xiaohua(Eric) XU, July 7, 2016
#
# Script for merging 3 subswaths TOPS interferograms and then unwrap and geocode. 
#
  if ($#argv != 2 && $#argv != 3) then
    echo ""
    echo "Usage: merge_unwrap_geocode_tops.csh inputfile config_file [det_stitch]"
    echo ""
    echo "Note: Inputfiles should be as following:"
    echo ""
    echo "      Swath1_Path:Swath1_master.PRM:Swath1_repeat.PRM"
    echo "      Swath2_Path:Swath2_master.PRM:Swath2_repeat.PRM"
    echo "      Swath3_Path:Swath3_master.PRM:Swath3_repeat.PRM"
    echo "      (Use the repeat PRM which contains the shift information.)"
    echo "      e.g. ../F1/intf/2015016_2015030/:S1A20151012_134357_F1.PRM"
    echo ""
    echo "      Make sure under each path, the processed phasefilt.grd, corr.grd and mask.grd exist."
    echo "      Also make sure the dem.grd is linked. "
    echo ""
    echo "      config_file is the same one used for processing."
    echo ""
    echo "      set det_stitch to 1 if you want to calculate stitching position based on the NaN area in the grids"
    echo ""
    echo "Example: merge_unwrap_geocode_tops.csh filelist batch.config"
    echo ""
    exit 1
  endif

  if (-f tmp_phaselist) rm tmp_phaselist
  if (-f tmp_corrlist) rm tmp_corrlist
  if (-f tmp_masklist) rm tmp_masklist

  if (! -f dem.grd ) then
    echo "Please link dem.grd to current folder"
    exit 1
  endif

  if ($#argv == 3) then
    set det_stitch = 1
  else
    set det_stitch = 0
    set n1 = 0
    set n2 = 0
  endif

  set region_cut = `grep region_cut $2 | awk '{print $3}'`

  # Creating inputfiles for merging
  foreach line (`awk '{print $0}' $1`)
    set now_dir = `pwd`
    set pth = `echo $line | awk -F: '{print $1}'`
    set prm = `echo $line | awk -F: '{print $2}'`
    set prm2 = `echo $line | awk -F: '{print $3}'`
    cd $pth
    set rshift = `grep rshift $prm2 | tail -1 | awk '{print $3}'`
    set fs1 = `grep first_sample $prm | awk '{print $3}'`
    set fs2 = `grep first_sample $prm2 | awk '{print $3}'`
    cp $prm tmp.PRM
    if ($fs2 > $fs1) then
      update_PRM tmp.PRM first_sample $fs2
    endif
    update_PRM tmp.PRM rshift $rshift
    cd $now_dir

    echo $pth"tmp.PRM:"$pth"phasefilt.grd" >> tmp_phaselist
    echo $pth"tmp.PRM:"$pth"corr.grd" >> tmp_corrlist
    echo $pth"tmp.PRM:"$pth"mask.grd" >> tmp_masklist
  end 

  set pth = `awk -F: 'NR==1 {print $1}' $1`
  set stem = `awk -F: 'NR==1 {print $2}' $1 | awk -F"." '{print $1}'`
  #echo $pth $stem

  echo ""
  echo "Merging START"
  if ($det_stitch == 1) then
    echo "Calculating valid starting columns of data ..."
    set nl = `wc -l $1 | awk '{print $1}'`
    if ($nl == 2) then
      set pth2 = `head -1 $1 | awk -F: '{print $1}'`
      gmt grdcut $pth2"phasefilt.grd" -Z+N -Gtmp.grd
      set xm1 = `gmt grdinfo $pth2"phasefilt.grd" -C | awk '{print $3}'`
      set xc1 = `gmt grdinfo tmp.grd -C | awk '{print $3}'`
      set incx = `gmt grdinfo tmp.grd -C | awk '{print $8}'`
      set n12 = `echo $xm1 $xc1 $incx | awk '{printf("%d",($1-$2)/$3)}'`
 
      set pth2 = `tail -1 $1 | awk -F: '{print $1}'`
      gmt grdcut $pth2"phasefilt.grd" -Z+N -Gtmp.grd
      set x01 = `gmt grdinfo tmp.grd -C | awk '{print $2}'`
      set incx = `gmt grdinfo tmp.grd -C | awk '{print $8}'`
      set n21 = `echo $x01 $incx | awk '{printf("%d",$1/$2)}'`
      set n1 = `echo $n12 $n21 | awk '{printf("%d",($1+$2)/2)}'`
      set n2 = 0
      rm tmp.grd
    else if ($nl == 3) then
      set pth2 = `head -1 $1 | awk -F: '{print $1}'`
      gmt grdcut $pth2"phasefilt.grd" -Z+N -Gtmp.grd
      set xm1 = `gmt grdinfo $pth2/phasefilt.grd -C | awk '{print $3}'`
      set xc1 = `gmt grdinfo tmp.grd -C | awk '{print $3}'`
      set incx = `gmt grdinfo tmp.grd -C | awk '{print $8}'`
      set n12 = `echo $xm1 $xc1 $incx | awk '{printf("%d",($1-$2)/$3)}'`

      set pth2 = `head -2 $1 | tail -1 | awk -F: '{print $1}'`
      gmt grdcut $pth2"phasefilt.grd" -Z+N -Gtmp.grd
      set x02 = `gmt grdinfo tmp.grd -C | awk '{print $2}'`
      set incx = `gmt grdinfo tmp.grd -C | awk '{print $8}'`
      set n21 = `echo $x02 $incx | awk '{printf("%d",$1/$2)}'`
      set n1 = `echo $n12 $n21 | awk '{printf("%d",($1+$2)/2)}'`
      set xm2 = `gmt grdinfo $pth2/phasefilt.grd -C | awk '{print $3}'`
      set xc2 = `gmt grdinfo tmp.grd -C | awk '{print $3}'`
      set n22 = `echo $xm1 $xc1 $incx | awk '{printf("%d",($1-$2)/$3)}'`

      set pth2 = `tail -1 $1 | awk -F: '{print $1}'`
      gmt grdcut $pth2"phasefilt.grd" -Z+N -Gtmp.grd
      set x03 = `gmt grdinfo tmp.grd -C | awk '{print $2}'`
      set incx = `gmt grdinfo tmp.grd -C | awk '{print $8}'`
      set n31 = `echo $x03 $incx | awk '{printf("%d",$1/$2)}'`
      set n2 = `echo $n22 $n31 | awk '{printf("%d",($1+$2)/2)}'`
      rm tmp.grd
    else
      echo "Incorrect number of records in input filelist .."
      exit 1
    endif
    echo "Stitching postitions set to $n1 $n2"
  endif
  
  if ($n1 > 5 && $n2 > 5) then
    merge_swath tmp_phaselist phasefilt.grd $stem $n1 $n2> merge_log
    merge_swath tmp_corrlist corr.grd $n1 $n2 > merge_log_corr
    merge_swath tmp_masklist mask.grd $n1 $n2 > merge_log_mask
  else
    merge_swath tmp_phaselist phasefilt.grd $stem > merge_log
    merge_swath tmp_corrlist corr.grd  > merge_log_corr
    merge_swath tmp_masklist mask.grd  > merge_log_mask
  endif
    
  echo "Merging END"
  echo ""

  set iono = `grep correct_iono $2 | awk '{print $3}'`
  set skip_iono = `grep iono_skip_est $2 | awk '{print $3}'`
  if ($iono != 0 & $skip_iono == 0) then
    if (! -f ph_iono_orig.grd) then
      echo "Need ph_iono_orig.grd to correct ionosphere ..."
    else
      echo "Correcting ionosphere ..."
      gmt grdsample ph_iono_orig.grd -Rphasefilt.grd -Gtmp.grd
      gmt grdmath phasefilt.grd tmp.grd SUB PI ADD 2 PI MUL MOD PI SUB = tmp2.grd
      mv phasefilt.grd phasefilt_orig.grd
      mv tmp2.grd phasefilt.grd
      rm tmp.grd
    endif
  endif

  
  # This step is essential, cut the DEM so it can run faster.
  if (! -f trans.dat) then
    set led = `grep led_file $pth$stem".PRM" | awk '{print $3}'`
    cp $pth$led .
    echo "Recomputing the projection LUT..."
  # Need to compute the geocoding matrix with supermaster.PRM with rshift set to 0
    set rshift = `grep rshift $stem".PRM" | tail -1 | awk '{print $3}'`
    update_PRM $stem".PRM" rshift 0
    gmt grd2xyz --FORMAT_FLOAT_OUT=%lf dem.grd -s | SAT_llt2rat $stem".PRM" 1 -bod > trans.dat
  # Set rshift back for other usage
    update_PRM $stem".PRM" rshift $rshift
  endif

  # Read in parameters
  set threshold_snaphu = `grep threshold_snaphu $2 | awk '{print $3}'`
  set threshold_geocode = `grep threshold_geocode $2 | awk '{print $3}'`
  set region_cut = `grep region_cut $2 | awk '{print $3}'`
  set switch_land = `grep switch_land $2 | awk '{print $3}'`
  set defomax = `grep defomax $2 | awk '{print $3}'`
  set near_interp = `grep near_interp $2 | awk '{print $3}'`
  set mask_water = `grep mask_water $2 | awk '{print $3}'`

  # Unwrapping
  if ($region_cut == "") then
    set region_cut = `gmt grdinfo phasefilt.grd -I- | cut -c3-20`
  endif
  if ($threshold_snaphu != 0 ) then
    if ($mask_water == 1 || $switch_land == 1) then
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

  # Geocoding 
  #if (-f raln.grd) rm raln.grd
  #if (-f ralt.grd) rm ralt.grd
 
  if ($threshold_geocode != 0) then
    echo ""
    echo "GEOCODE-START"
    proj_ra2ll.csh trans.dat phasefilt.grd phasefilt_ll.grd
    gmt grdmath corr.grd $threshold_geocode GE 0 NAN mask.grd MUL = mask2.grd
    gmt grdmath phasefilt.grd mask2.grd MUL = phasefilt_mask.grd
    proj_ra2ll.csh trans.dat phasefilt_mask.grd phasefilt_mask_ll.grd
    proj_ra2ll.csh trans.dat corr.grd corr_ll.grd
    gmt makecpt -Crainbow -T-3.15/3.15/0.05 -Z > phase.cpt
    set BT = `gmt grdinfo -C corr.grd | awk '{print $7}'`
    gmt makecpt -Cgray -T0/$BT/0.05 -Z > corr.cpt
    grd2kml.csh phasefilt_ll phase.cpt
    grd2kml.csh phasefilt_mask_ll phase.cpt
    grd2kml.csh corr_ll corr.cpt

    if (-f unwrap.grd) then
      gmt grdmath unwrap.grd mask2.grd MUL = unwrap_mask.grd
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

  rm tmp_phaselist tmp_corrlist tmp_masklist *.eps *.bb
