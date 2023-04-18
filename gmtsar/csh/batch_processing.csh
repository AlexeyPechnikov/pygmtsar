#!/bin/csh -f

  if ($#argv != 5 && $#argv != 4) then
    echo ""
    echo "Usage: batch_processing.csh SAT master_image inputfile step [configuration_file] "
    echo ""
    echo "Example: batch_processing.csh ALOS IMG-HH-ALPSRP055750660-H1.0__A imagelist 1 [config.alos.txt]"
    echo ""
    echo "    Put the data and orbit files in the raw folder, put DEM in the topo folder"
    echo "    The SAT needs to be specified, choices with in ERS, ENVI, ALOS, ALOS_SLC, ALOS2, ALOS2_SCAN"
    echo "    S1_STRIP, S1_TOPS, ENVI_SLC, CSK_RAW, CSK_SLC, TSX, RS2, GF3"
    echo ""
    echo "    Make sure the files from the same date have the same stem, e.g. stem.tif stem.xml stem.cos stem.EOF, etc"
    echo ""
    echo "    If the configuration file is left blank, the program will generate one "
    echo "    with default parameters "
    echo ""
    echo "Details: The step could be 1-preprocessing, 2-alignment, 3-backgeocoding, 4-interferometry, "
    echo "    5-phase unwrapping, and 6-geocoding. For TOPS data, step 1 and 2 are done together." 
    echo ""
    echo "    step 1: preprocessing, inputfile should be a list of inputdata with or without orbit files, master_image should be in the first line"
    echo "    NEED at least THREE records to run "
    echo "    e.g., for ALOS data, the list should be: (pick FBS to be master image)"
    echo "        IMG-HH-ALPSRP160702940-H1.0__D "
    echo "        IMG-HH-ALPSRP200962940-H1.0__D "
    echo "        IMG-HH-ALPSRP268062940-H1.0__D "
    echo "    and for Sentinel-1 TOPS data, the list should be like:"
    echo "        s1a-iw1-slc-vv-20150109t134413-20150109t134421-004095-004f4a-001:S1A_OPER_AUX_POEORB_OPOD_20210305T105546_V20150108T225944_20150110T005944.EOF "
    echo "        s1a-iw1-slc-vv-20150121t134413-20150121t134421-004270-005317-001:S1A_OPER_AUX_POEORB_OPOD_20210305T143838_V20150120T225944_20150122T005944.EOF "
    echo "        s1a-iw1-slc-vv-20150226t134412-20150226t134420-004795-005f58-001:S1A_OPER_AUX_POEORB_OPOD_20210306T014915_V20150225T225944_20150227T005944.EOF "
    echo ""
    echo "    step 2: alignment, inputfile should be a list of inputdata with master image in the first line "
    echo "    e.g., for ALOS data, the list should be like:"
    echo "    NEED at least THREE records to run "
    echo "        IMG-HH-ALPSRP160702940-H1.0__D "
    echo "        IMG-HH-ALPSRP200962940-H1.0__D "
    echo "        IMG-HH-ALPSRP268062940-H1.0__D "
    echo "    and for Sentinel-1 TOPS data, this step could be skipped."
    echo "    For secondary alignment, run this step multiple times with -skip_master = 1, for secondary and tertiary alignment."
    echo ""
    echo "    step 3: backgeocoding, make sure DEM is in the topo directory, and PRM and LED file of the master image exist in the raw directory. "
    echo "    This only need to be run once for a stack of data."
    echo ""
    echo "    step 4: interferometry, inputfile should be a list of interfer "
    echo "    e.g., for ALOS data the list should be like:"
    echo "        IMG-HH-ALPSRP160702940-H1.0__D:IMG-HH-ALPSRP200962940-H1.0__D "
    echo "        IMG-HH-ALPSRP160702940-H1.0__D:IMG-HH-ALPSRP268062940-H1.0__D "
    echo "        IMG-HH-ALPSRP200962940-H1.0__D:IMG-HH-ALPSRP268062940-H1.0__D"
    echo "    and for Sentinel-1 TOPS data, the list should be like:"
    echo "        S1_20150121_ALL_F1:S1_20150310_ALL_F1"
    echo "        S1_20150121_ALL_F1:S1_20150403_ALL_F1"
    echo "        S1_20150310_ALL_F1:S1_20150403_ALL_F1"
    echo ""
    echo "    setp 5: phase unwrapping, give the same input as step 4, and the script will go through every directory in intf_all and run "
    echo "    phase unwrapping using parameters specified in the config file"
    echo ""
    echo "    step 6: geocoding, give the same input as step 4, and the script will go through every directory in intf_all and run geocoding"
    echo ""
    echo "Note: for data that comes with multiple subswaths, this script only works with one subswath. One'll need to perform merging first "
    echo "    and then unwrap or geocode."
    echo "    Also, this script does data processing only. For time-series analysis, refer to the sbas or sbas_parallel program, or use external"
    echo "    software to construct time-series."
    echo ""
    echo ""
    exit 1
  endif


# start 
#  Make sure the config exist
  if ($#argv == 5) then
    if(! -f $5 ) then
      echo " [ERROR]: no configure file: "$5
      echo " Leave it blank to generate config file with default values."
      exit 1
    endif
  endif

  if ($4 != 1 && $4 != 2 && $4 != 3 && $4 != 4 && $4 != 5 && $4 != 6) then
    echo "[ERROR]: Wrong step input"
    exit 1
  endif

#
#  Read parameters from the configure file
#
  set SAT = `echo $1`
  if ($#argv == 5) then
    set conf = `echo $5`
  else
    pop_config.csh $SAT > config.$SAT.txt
    set conf = `echo "config.$SAT.txt"`
  endif

  set stage = $4
    set num_patches = `grep num_patches $conf | awk '{print $3}'`
  set near_range = `grep near_range $conf | awk '{print $3}'`
  set earth_radius = `grep earth_radius $conf | awk '{print $3}'`
  set fd = `grep fd1 $conf | awk '{print $3}'`
  set topo_phase = `grep topo_phase $conf | awk '{print $3}'`
  set topo_interp_mode = `grep topo_interp_mode $conf | awk '{print $3}'`
  if ( "x$topo_interp_mode" == "x" ) then
    set topo_interp_mode = 0
  endif
  set shift_topo = `grep shift_topo $conf | awk '{print $3}'`
  set switch_master = `grep switch_master $conf | awk '{print $3}'`
  set filter = `grep filter_wavelength $conf | awk '{print $3}'`
  set compute_phase_gradient = `grep compute_phase_gradient $conf | awk '{print $3}'`
  set iono = `grep correct_iono $conf | awk '{print $3}'`
  if ( "x$iono" == "x" ) then
    set iono = 0
  endif
  set iono_filt_rng = `grep iono_filt_rng $conf | awk '{print $3}'`
  set iono_filt_azi = `grep iono_filt_azi $conf | awk '{print $3}'`
  set iono_dsamp = `grep iono_dsamp $conf | awk '{print $3}'`
  set iono_skip_est = `grep iono_skip_est $conf | awk '{print $3}'`
  set spec_div = `grep spec_div $conf | awk '{print $3}'`
  if ( "x$spec_div" == "x" ) then
    set spec_div = 0
  endif
  set spec_mode = `grep spec_mode $conf | awk '{print $3}'`
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
  set switch_land = `grep switch_land $conf | awk '{print $3}'`
  set defomax = `grep defomax $conf | awk '{print $3}'`
  set range_dec = `grep range_dec $conf | awk '{print $3}'`
  set azimuth_dec = `grep azimuth_dec $conf | awk '{print $3}'`
  set SLC_factor = `grep SLC_factor $conf | awk '{print $3}'`
  set near_interp = `grep near_interp $conf | awk '{print $3}'`
  set geometric_coreg = `grep geometric_coreg $conf | awk '{print $3}'`
  set master = ` echo $2`
  set inputlist =  ` echo $3 `
  echo ""


#############################
# 1 - start from preprocess #
#############################

  if ($stage == 1) then

    echo ""
    echo "PREPROCESS BATCH - START"
    echo ""
    echo "  DO MAKE SURE Every needed image and orbit file exist in the raw directory"
    echo ""
    echo "  Cleaning raw directory, removing all .raw, .PRM, .LED and .SLC files"
    echo ""

    cd raw
    rm *.PRM* *.raw *.LED *.SLC
    cd ..

    set nl = `wc -l $inputlist | awk '{print $1}'`
    if ($nl < 3) then
      echo ""
      echo "NEED at least THREE records to run batch preprocessing"
      echo ""
      exit 1
    endif
#
#  Start preprocessing
#
#    if ($SAT == "S1_TOPS") then
#      ln -s ../topo/dem.grd .
#      preproc_batch_tops.csh ../$inputlist dem.grd 1
#      echo "Finished preprocessing $nl $SAT images "
#      exit 1
#    endif

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
    if (!($spec_div == 0)) then
      set commandline = "$commandline -ESD $spec_mode"
    endif


    set ii = 0
    foreach aligned (`cat $inputlist`)
      if ($SAT == "S1_TOPS") then
        cd raw
        set orb = `echo $aligned | awk -F':' '{print $2}'`
        set aligned = `echo $aligned | awk -F':' '{print $1}'`
        ln -s $orb ./$aligned".EOF"
        cd ..
      endif
      if ($ii == 0) then
      else if ($ii == 1) then
        sed "s/.*proc_stage.*/proc_stage = 1/g" $conf | sed "s/.*skip_stage.*/skip_stage = 2,3,4,5,6/g" | sed "s/.*skip_master.*/skip_master = 0/g" > tmp_conf_$aligned
        p2p_processing.csh $SAT $master $aligned tmp_conf_$aligned
        #echo "pre_proc.csh $SAT $master $aligned $commandline -skip_master 0"
        #pre_proc.csh $SAT $master $aligned $commandline -skip_master 0
        rm tmp_conf_$aligned
      else
        sed "s/.*proc_stage.*/proc_stage = 1/g" $conf | sed "s/.*skip_stage.*/skip_stage = 2,3,4,5,6/g" | sed "s/.*skip_master.*/skip_master = 1/g" > tmp_conf_$aligned
        p2p_processing.csh $SAT $master $aligned tmp_conf_$aligned
        #echo "pre_proc.csh $SAT $master $aligned $commandline -skip_master 1"
        #pre_proc.csh $SAT $master $aligned $commandline -skip_master 1
        rm tmp_conf_$aligned
      endif
      set ii = `echo $ii | awk '{print $ii+1}'`
    end
    cd raw
    ls *.PRM > prmlist
    if ($SAT == "S1_TOPS") then
      set mmaster = ` echo $master | awk '{ print "S1_"substr($1,16,8)"_"substr($1,25,6)"_F"substr($1,7,1)}'`
    endif
    get_baseline_table.csh prmlist $mmaster".PRM"
    cd ..

    echo ""
    echo "PREPROCESS BATCH - END"
    echo ""
    echo "Finished preprocessing $ii $SAT images "
    echo ""
    exit 1
  endif

#############################################
# 2 - start from focus and align SLC images #
#############################################
  
  if ($stage == 2) then

#
# focus and align SLC images 
# 
    echo ""
    echo "ALIGN BATCH - START"
    echo ""

    set ii = 0
    if ($geometric_coreg == 0 || $SAT == "S1_TOPS" ) then
      foreach aligned (`cat $inputlist`)
        if ($SAT == "S1_TOPS") then
          set orb = `echo $aligned | awk -F':' '{print $2}'`
          set aligned = `echo $aligned | awk -F':' '{print $1}'`
          echo "Skipping TOPS data $aligned, as alignment is done in the first step.."
        endif
        if ($ii == 0) then
        else if ($ii == 1) then
          sed "s/.*proc_stage.*/proc_stage = 2/g" $conf | sed "s/.*skip_stage.*/skip_stage = 3,4,5,6/g" | sed "s/.*skip_master.*/skip_master = 0/g" > tmp_conf_$aligned
          p2p_processing.csh $SAT $master $aligned tmp_conf_$aligned
          rm tmp_conf_$aligned
        else 
          sed "s/.*proc_stage.*/proc_stage = 2/g" $conf | sed "s/.*skip_stage.*/skip_stage = 3,4,5,6/g" | sed "s/.*skip_master.*/skip_master = 1/g" > tmp_conf_$aligned
          p2p_processing.csh $SAT $master $aligned tmp_conf_$aligned
          rm tmp_conf_$aligned
        endif
        set ii = `echo $ii | awk '{print $ii+1}'`
      end
    else
      if ($SAT == "ERS" || $SAT == "ENVI" || $SAT == "ALOS" || $SAT == "CSK_RAW") then
        align_batch.csh RAW 1 $inputlist
      else 
        align_batch.csh SLC 1 $inputlist
      endif
    endif

    echo ""
    echo "ALIGN BATCH - END"
    echo ""
    echo "Finished aligning $ii $SAT images"
    echo ""
    exit 1
  endif

##################################
# 3 - start from make topo_ra    #
##################################

  if ($stage == 3) then

#
# back-geocode dem.grd and make topo_ra
#
    echo ""
    echo "BACKGEOCODING - START"
    echo ""
    sed "s/.*proc_stage.*/proc_stage = 3/g" $conf | sed "s/.*skip_stage.*/skip_stage = 4,5,6/g" > tmp_conf_$master
    p2p_processing.csh $SAT $master $master tmp_conf_$master
    rm tmp_conf_$master
    echo ""
    echo "BACKGEOCODING - END"
    echo ""
  endif

##################################################
# 4 - start from make and filter interferograms  #
##################################################

  if ($stage == 4) then

    echo ""
    echo "INTERFEROMETRY BATCH - START"
    echo "" 
    mkdir -p intf_all
    sed "s/.*proc_stage.*/proc_stage = 4/g" $conf > tmp_conf_$master
    foreach pair (`cat $inputlist`)
      set ref = `echo $pair | awk -F: '{print $1}'`
      set rep = `echo $pair | awk -F: '{print $2}'`
      set ref_id  = `grep SC_clock_start ./raw/$ref.PRM | awk '{printf("%d",int($3))}' `
      set rep_id  = `grep SC_clock_start ./raw/$rep.PRM | awk '{printf("%d",int($3))}' `
      if ($SAT == "S1_TOPS") then
        set tref = `ls raw/$ref".PRM" | awk -F'/' '{print substr($2,4,8)}'`
        set ref = `ls raw/*$tref*.xml | awk -F'/' '{print substr($2,1,length($2)-4)}'`
        set trep = `ls raw/$rep".PRM" | awk -F'/' '{print substr($2,4,8)}'`
        set rep = `ls raw/*$trep*.xml | awk -F'/' '{print substr($2,1,length($2)-4)}'`
      endif
      p2p_processing.csh $SAT $ref $rep tmp_conf_$master
      if (-e intf_all/$ref_id"_"$rep_id) rm -rf intf_all/$ref_id"_"$rep_id
      mv intf/$ref_id"_"$rep_id intf_all/$ref_id"_"$rep_id
    end
    rm tmp_conf_$master
    echo ""
    echo "INTERFEROMETRY BATCH - END"
    echo "" 

  endif




