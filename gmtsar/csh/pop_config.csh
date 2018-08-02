#!/bin/csh -f

# create a configure file for p2p_processing.csh
# syntax: pop_config.csh SAT
# SAT can be  ERS, ENVI, ALOS, ALOS_SLC, ALOS2, ALOS2_SCAN
# S1_STRIP, S1_TOPS,, CSK_RAW, CSK_SLC, TSX, RS2

if ($#argv != 1) then
  echo ""
  echo "Usage: pop_config.csh SAT"
  echo ""
  echo "       SAT can be ERS, ENVI, ALOS, ALOS_SLC, ALOS2, ALOS2_SCAN"
  echo "       S1_STRIP, S1_TOPS, CSK_RAW, CSK_SLC, TSX, RS2"
  echo ""
  exit 1
endif

set SAT = `echo $1`

echo "#"
echo "# This is an example configuration file for p2p_processing.csh"
echo "#"
echo "# all the comments or explanations are marked by "\""#"\" 
echo "# The parameters in this configuration file is distinguished by their first word so "
echo "# user should follow the naming of each parameter."
echo "# the parameter name, "\""="\"" sign, parameter value should be separated by space "\"" "\"". "
echo "# leave the parameter value blank if using default value. "
echo "#"
echo "  "
echo "#####################"
echo "# processing stage  #"
echo "#####################"
echo "# 1 - start from preprocess"
echo "# 2 - start from align SLC images"
echo "# 3 - start from make topo_ra "
echo "# 4 - start from make and filter interferograms "
echo "# 5 - start from unwrap phase"
echo "# 6 - start from geocode  "
echo "proc_stage = 1"
echo ""
echo "##################################"
echo "#   parameters for preprocess    #"
echo "#   - pre_proc.csh               #"
echo "##################################"
echo "# num of patches"
echo "num_patches = "
echo ""
echo "# earth radius "
echo "earth_radius ="
echo ""
echo "# near_range"
echo "near_range = "
echo ""
echo "# Doppler centroid "
echo "fd1 = "
echo ""

if ($SAT == "ALOS_SLC") then
  echo "# SLC scale factor to convert float to int "
  echo "SLC_factor = 0.02"
  echo ""
else 
  if ($SAT == "ALOS2" || $SAT == "ALOS2_SCAN") then
    echo "# SLC scale factor to convert float to int"
    echo "SLC_factor = 2.0"
    echo ""
  endif
endif

echo "################################################"
echo "#   parameters for focus and align SLC images  #"
echo "#   - align.csh                                #"
echo "################################################"
echo "#"
echo "#####################################"
echo "#   parameters for make topo_ra     #"
echo "#   - dem2topo_ra.csh               #"
echo "#####################################"
echo "# subtract topo_ra from the phase"
echo "#  (1 -- yes; 0 -- no)"
echo "topo_phase = 1"
echo "# if above parameter = 1 then one should have put dem.grd in topo/"
echo ""
echo "# topo_ra shift (1 -- yes; 0 -- no)"

if ($SAT == "ALOS_SLC" || $SAT == "ALOS" || $SAT == "ERS") then
  echo "shift_topo = 1"
else 
  echo "shift_topo = 0"
endif

echo ""
echo "####################################################"
echo "#   parameters for make and filter interferograms  #"
echo "#   - intf.csh                                     #"
echo "#   - filter.csh                                   #"
echo "####################################################"
echo "# switch the master and slave when doing intf. "
echo "# put "\""1"\"" if assume master as repeat and slave as reference "
echo "# put "\""0"\"" if assume master as reference and slave as repeat [Default]"
echo "# phase = repeat phase - reference phase"
echo "switch_master = 0"
echo ""
echo "# filters "
echo "# look at the filter/ folder to choose other filters"
echo "# for tops processing, to force the decimation factor"
echo "# recommended range decimation to be 8, azimuth decimation to be 2"
if ($SAT == "ALOS2_SCAN") then
  echo "filter_wavelength = 400"
else if ($SAT == "RS2" || $SAT == "TSX") then
  echo "filter_wavelength = 100"
else
  echo "filter_wavelength = 200"
endif
echo ""
echo "# decimation of images "
echo "# decimation control the size of the amplitude and phase images. It is either 1 or 2."
echo "# Set the decimation to be 1 if you want higher resolution images."
echo "# Set the decimation to be 2 if you want images with smaller file size."
echo "# "
if ($SAT == "RS2" || $SAT == "TSX") then
  echo "dec_factor = 1 "
else if ($SAT == "ALOS2_SCAN") then
  echo "dec_factor = 4 "
else
  echo "dec_factor = 2 "
endif
if ($SAT == "S1_TOPS") then
  echo "range_dec = 8"
  echo "azimuth_dec = 2"
endif
echo "#"
echo "# "
echo "#####################################"
echo "#   parameters for unwrap phase     #"
echo "#   - snaphu.csh                    #"
echo "#####################################"
echo "# correlation threshold for snaphu.csh (0~1)"
echo "# set it to be 0 to skip unwrapping."
echo "threshold_snaphu = 0"
echo ""
echo "# interpolate masked or low coherence pixels with their nearest neighbors, 1 means interpolate, "
echo "# others or blank means using original phase, see snaphu.csh and snaphu_interp.csh for details"
echo "# this could be very slow in case a large blank area exist"
echo "near_interp = 1"
echo ""
echo "# region to unwrap in radar coordinates (leave it blank if unwrap the whole region)"
echo "# example 300/5900/0/25000"
echo "region_cut ="
echo ""
echo "# mask the wet region (Lakes/Oceans) before unwrapping (1 -- yes; else -- no)"
echo "mask_water = 1"
echo ""
echo "#"
echo "# Allow phase discontinuity in unrapped phase. This is needed for interferograms having sharp phase jumps."
echo "# defo_max = 0 - used for smooth unwrapped phase such as interseismic deformation"
echo "# defo_max = 65 - will allow a phase jump of 65 cycles or 1.82 m of deformation at C-band"
echo "#"
echo "defomax = 0"
echo ""
echo "#####################################"
echo "#   parameters for geocode          #"
echo "#   - geocode.csh                   #"
echo "#####################################"
echo "# correlation threshold for geocode.csh (0< threshold <=1), set 0 to skip"
echo "threshold_geocode = .10"
echo ""
