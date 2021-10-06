#!/bin/csh -f
#
#  Xiaohua XU, Mar 10, 2017
#
# process Sentinel-1A TOPS data
# Automatically process a single frame of interferogram.
# see instruction.txt for details.
#

alias rm 'rm -f'
unset noclobber
#
  if ($#argv !=  7) then
    echo ""
    echo "Usage: p2p_S1_TOPS_Frame.csh Master.SAFE Master.EOF Aligned.SAFE Aligned.EOF config.s1a.txt polarization parallel"
    echo ""
    echo "Example: p2p_S1_TOPS_Frame.csh S1A_IW_SLC__1SDV_20150607T014936_20150607T015003_006261_00832E_3626.SAFE S1A_OPER_AUX_POEORB_OPOD_20150615T155109_V20150525T225944_20150527T005944.EOF S1A_IW_SLC__1SSV_20150526T014935_20150526T015002_006086_007E23_679A.SAFE S1A_OPER_AUX_POEORB_OPOD_20150627T155155_V20150606T225944_20150608T005944.EOF config.s1a.txt vv 1"
    echo ""
    echo "    Place the .SAFE file in the raw folder, DEM in the topo folder"
    echo "    During processing, F1, F2, F3 and merge folder will be generated"
    echo "    Final results will be placed in the merge folder, with phase"	
    echo "    corr [unwrapped phase]."
    echo "    polarization = vv vh hh or hv "
    echo "    parallel = 0-sequential  1-parallel "
    echo ""
    echo "Reference: Xu, X., Sandwell, D.T., Tymofyeyeva, E., González-Ortega, A. and Tong, X., "
    echo "    2017. Tectonic and Anthropogenic Deformation at the Cerro Prieto Geothermal "
    echo "    Step-Over Revealed by Sentinel-1A InSAR. IEEE Transactions on Geoscience and Remote Sensing."
    echo ""
    exit 1
  endif
# start 
# 
# set polarization
# 
  set pol = $6 
  echo $pol
#
# set processing mode seq
#
  set seq = $7
  echo $seq
#:<<supercalifragilisticexpialidocious
#
# determine file names
#
  set pth = `pwd`
  cd raw/$1
  set f1m = `ls */*iw1*$pol*xml | awk '{print substr($1,12,length($1)-15)}'`
  set f2m = `ls */*iw2*$pol*xml | awk '{print substr($1,12,length($1)-15)}'`
  set f3m = `ls */*iw3*$pol*xml | awk '{print substr($1,12,length($1)-15)}'`
  cd ../$3
  set f1s = `ls */*iw1*$pol*xml | awk '{print substr($1,12,length($1)-15)}'`
  set f2s = `ls */*iw2*$pol*xml | awk '{print substr($1,12,length($1)-15)}'`
  set f3s = `ls */*iw3*$pol*xml | awk '{print substr($1,12,length($1)-15)}'`
  cd $pth
#if (6 == 9) then
#
# organize files
#
  mkdir F1
  mkdir F1/raw F1/topo
  cd F1
  sed "s/.*threshold_geocode.*/threshold_geocode = 0/g" ../$5 | sed "s/.*threshold_snaphu.*/threshold_snaphu = 0/g" | sed "s/.*iono_skip_est.*/iono_skip_est = 1/g"> $5
  echo "Linking files for Subswath 1 ..."
  cd topo
  ln -s ../../topo/dem.grd .
  cd ../raw
  ln -s ../topo/dem.grd .
  ln -s ../../raw/$1/*/$f1m.xml .
  ln -s ../../raw/$1/*/$f1m.tiff .
  ln -s ../../raw/$2 ./$f1m.EOF
  ln -s ../../raw/$3/*/$f1s.xml .
  ln -s ../../raw/$3/*/$f1s.tiff .
  ln -s ../../raw/$4 ./$f1s.EOF
  cd ../..

  mkdir F2
  mkdir F2/raw F2/topo
  cd F2
  sed "s/.*threshold_geocode.*/threshold_geocode = 0/g" ../$5 | sed "s/.*threshold_snaphu.*/threshold_snaphu = 0/g" | sed "s/.*iono_skip_est.*/iono_skip_est = 1/g"> $5
  echo "Linking files for Subswath 2 ..."
  cd topo
  ln -s ../../topo/dem.grd .
  cd ../raw
  ln -s ../topo/dem.grd .
  ln -s ../../raw/$1/*/$f2m.xml .
  ln -s ../../raw/$1/*/$f2m.tiff .
  ln -s ../../raw/$2 ./$f2m.EOF
  ln -s ../../raw/$3/*/$f2s.xml .
  ln -s ../../raw/$3/*/$f2s.tiff .
  ln -s ../../raw/$4 ./$f2s.EOF
  cd ../..

  mkdir F3
  mkdir F3/raw F3/topo
  cd F3
  echo "Linking files for Subswath 3 ..."
  sed "s/.*threshold_geocode.*/threshold_geocode = 0/g" ../$5 | sed "s/.*threshold_snaphu.*/threshold_snaphu = 0/g" | sed "s/.*iono_skip_est.*/iono_skip_est = 1/g"> $5
  cd topo
  ln -s ../../topo/dem.grd .
  cd ../raw
  ln -s ../topo/dem.grd .
  ln -s ../../raw/$1/*/$f3m.xml .
  ln -s ../../raw/$1/*/$f3m.tiff .
  ln -s ../../raw/$2 ./$f3m.EOF
  ln -s ../../raw/$3/*/$f3s.xml .
  ln -s ../../raw/$3/*/$f3s.tiff .
  ln -s ../../raw/$4 ./$f3s.EOF
  cd ../..
# 
# process data
# 
  if ($seq == 0) then
    cd F1
    p2p_processing.csh S1_TOPS $f1m $f1s $5
    cd ../F2
    p2p_processing.csh S1_TOPS $f2m $f2s $5
    cd ../F3
    p2p_processing.csh S1_TOPS $f3m $f3s $5
    cd ..
  else if ($seq == 1) then
    cd F1
    p2p_processing.csh S1_TOPS $f1m $f1s $5 >&log&
    cd ../F2
    p2p_processing.csh S1_TOPS $f2m $f2s $5 >&log&
    cd ../F3
    p2p_processing.csh S1_TOPS $f3m $f3s $5 >&log&
    cd ..
    wait
  else
    echo "Invalid mode"
    exit
  endif
#supercalifragilisticexpialidocious
#
# merge_unwrap_geocode
#
#endif

  mkdir merge
  cd merge
  ln -s ../topo/dem.grd .
  ln -s ../F1/intf/*/gauss* .
  if (-f tmp.filelist) then
    rm tmp.filelist
  endif
  set pth1 = `ls ../F1/intf/*/*PRM | awk NR==1'{print $1}' | awk -F"/" '{for (i=1;i<NF;i++) printf("%s/",$i)}'`
  set prm1m = `echo $f1m | awk '{ print "S1_"substr($1,16,8)"_"substr($1,25,6)"_F"substr($1,7,1)".PRM"}'`
  set prm1s = `echo $f1s | awk '{ print "S1_"substr($1,16,8)"_"substr($1,25,6)"_F"substr($1,7,1)".PRM"}'`

  #set prm1m = `ls ../F1/intf/*/*PRM | awk NR==1'{print $1}' | awk -F"/" '{print $NF}'`
  #set prm1s = `ls ../F1/intf/*/*PRM | awk NR==2'{print $1}' | awk -F"/" '{print $NF}'`
  echo $pth1":"$prm1m":"$prm1s > tmp.filelist
  
  set pth2 = `ls ../F2/intf/*/*PRM | awk NR==1'{print $1}' | awk -F"/" '{for (i=1;i<NF;i++) printf("%s/",$i)}'`
  set prm2m = `echo $f2m | awk '{ print "S1_"substr($1,16,8)"_"substr($1,25,6)"_F"substr($1,7,1)".PRM"}'`
  set prm2s = `echo $f2s | awk '{ print "S1_"substr($1,16,8)"_"substr($1,25,6)"_F"substr($1,7,1)".PRM"}'`

  #set prm2m = `ls ../F2/intf/*/*PRM | awk NR==1'{print $1}' | awk -F"/" '{print $NF}'`
  #set prm2s = `ls ../F2/intf/*/*PRM | awk NR==2'{print $1}' | awk -F"/" '{print $NF}'`
  echo $pth2":"$prm2m":"$prm2s >> tmp.filelist

  set pth3 = `ls ../F3/intf/*/*PRM | awk NR==1'{print $1}' | awk -F"/" '{for (i=1;i<NF;i++) printf("%s/",$i)}'`
  set prm3m = `echo $f3m | awk '{ print "S1_"substr($1,16,8)"_"substr($1,25,6)"_F"substr($1,7,1)".PRM"}'`
  set prm3s = `echo $f3s | awk '{ print "S1_"substr($1,16,8)"_"substr($1,25,6)"_F"substr($1,7,1)".PRM"}'`
  
  #set prm3m = `ls ../F3/intf/*/*PRM | awk NR==1'{print $1}' | awk -F"/" '{print $NF}'`
  #set prm3s = `ls ../F3/intf/*/*PRM | awk NR==2'{print $1}' | awk -F"/" '{print $NF}'`
  echo $pth3":"$prm3m":"$prm3s >> tmp.filelist


  set iono = `grep correct_iono ../$5 | awk '{print $3}'`
  if ($iono != 0) then
    sed "s/.*threshold_geocode.*/threshold_geocode = 0/g" ../$5 | sed "s/.*threshold_snaphu.*/threshold_snaphu = 0/g" | sed "s/.*iono_skip_est.*/iono_skip_est = 1/g" > $5
    merge_unwrap_geocode_tops.csh tmp.filelist $5

    cd ..
    mkdir iono
    cd iono
    mkdir intf_h intf_l intf_o iono_correction
    cd intf_h
    echo "../../F1/iono_phase/intf_h/:"$prm1m":"$prm1s > tmp.filelist
    echo "../../F2/iono_phase/intf_h/:"$prm2m":"$prm2s >> tmp.filelist
    echo "../../F3/iono_phase/intf_h/:"$prm3m":"$prm3s >> tmp.filelist
    sed "s/.*threshold_geocode.*/threshold_geocode = 0/g" ../../$5 | sed "s/.*threshold_snaphu.*/threshold_snaphu = 0.1/g" | sed "s/.*iono_skip_est.*/iono_skip_est = 1/g" > $5
    ln -s ../../topo/dem.grd .
    ln -s ../../merge/trans.dat .
    merge_unwrap_geocode_tops.csh tmp.filelist $5
    cp ../../F1/SLC/params* .
    cd ../intf_l

    echo "../../F1/iono_phase/intf_l/:"$prm1m":"$prm1s > tmp.filelist
    echo "../../F2/iono_phase/intf_l/:"$prm2m":"$prm2s >> tmp.filelist
    echo "../../F3/iono_phase/intf_l/:"$prm3m":"$prm3s >> tmp.filelist
    sed "s/.*threshold_geocode.*/threshold_geocode = 0/g" ../../$5 | sed "s/.*threshold_snaphu.*/threshold_snaphu = 0.1/g" | sed "s/.*iono_skip_est.*/iono_skip_est = 1/g" > $5
    ln -s ../../topo/dem.grd .
    ln -s ../../merge/trans.dat .
    merge_unwrap_geocode_tops.csh tmp.filelist $5
    cp ../../F1/SLC/params* .
    cd ../intf_o

    echo "../../F1/iono_phase/intf_o/:"$prm1m":"$prm1s > tmp.filelist
    echo "../../F2/iono_phase/intf_o/:"$prm2m":"$prm2s >> tmp.filelist
    echo "../../F3/iono_phase/intf_o/:"$prm3m":"$prm3s >> tmp.filelist
    sed "s/.*threshold_geocode.*/threshold_geocode = 0/g" ../../$5 | sed "s/.*threshold_snaphu.*/threshold_snaphu = 0.1/g" | sed "s/.*iono_skip_est.*/iono_skip_est = 1/g" > $5
    ln -s ../../topo/dem.grd .
    ln -s ../../merge/trans.dat .
    merge_unwrap_geocode_tops.csh tmp.filelist $5
    cp ../../F1/SLC/params* .
    cd ../iono_correction

    estimate_ionospheric_phase.csh ../intf_h ../intf_l ../intf_o ../../merge 0.8 0.8 
    cd ../../merge
    ln -s ../iono/iono_correction/ph_iono_orig.grd .
    cp ../$5 .

    merge_unwrap_geocode_tops.csh tmp.filelist $5
   
  else
    cp ../$5 .
    merge_unwrap_geocode_tops.csh tmp.filelist $5

  endif




