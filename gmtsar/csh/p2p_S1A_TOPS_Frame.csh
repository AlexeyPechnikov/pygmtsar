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
  if ($#argv < 5 || $#argv > 6) then
    echo ""
    echo "Usage: p2p_S1A_TOPS.csh Master.SAFE Master.EOF Slave.SAFE Slave.EOF config.s1a.txt mode"
    echo ""
    echo "Example: p2p_S1A_TOPS.csh S1A_IW_SLC__1SDV_20150607T014936_20150607T015003_006261_00832E_3626.SAFE S1A_OPER_AUX_POEORB_OPOD_20150615T155109_V20150525T225944_20150527T005944.EOF S1A_IW_SLC__1SSV_20150526T014935_20150526T015002_006086_007E23_679A.SAFE S1A_OPER_AUX_POEORB_OPOD_20150627T155155_V20150606T225944_20150608T005944.EOF config.s1a.txt "
    echo ""
    echo "	Place the .SAFE file in the raw folder, DEM in the topo folder"
    echo "	During processing, F1, F2, F3 and merge folder will be generated"
    echo "	Final results will be placed in the merge folder, with phase"	
    echo "	, corr and amp [unwrapped phase]."
    echo "	mode = 11/12, process vv, mode = 21/22, process vh, default is vv"
    echo "	mode = 11/21, process sequentially, mode = 12/22, process parallely, default is sequentially"
    echo ""
    exit 1
  endif
# start 
#
# determine mode
#

  if ($#argv == 5) then
    set md = 1
    set seq = 1
    echo "Processing VV data squentially..."
  else if ($6 == 11) then
    set md = 1
    set seq = 1
    echo "Processing VV data squentially..."
  else if ($6 == 21) then
    set md = 2
    set seq = 1
    echo "Processing VH data squentially..."
  else if ($6 == 12) then
    set md = 1
    set seq = 2
    echo "Processing VV data parallelly..."
  else if ($6 == 22) then
    set md = 2
    set seq = 2
    echo "Processing VH data parallelly..."
  else
    echo "Invalid mode ..."
    exit 1
  endif

# 
# set polarization
# 
  if ($md == 1) then
    set pol = vv
  else
    set pol = vh
  endif
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
#
# organize files
#
  mkdir F1
  mkdir F1/raw F1/topo
  cd F1
  sed "s/.*threshold_geocode.*/threshold_geocode = 0/g" ../$5 | sed "s/.*threshold_snaphu.*/threshold_snaphu = 0/g" > $5
  cd topo
  ln -s ../../topo/dem.grd .
  cd ../raw
  ln -s ../topo/dem.grd .
  ln -s ../../raw/$1/*/$f1m.xml .
  ln -s ../../raw/$1/*/$f1m.tiff .
  ln -s ../../raw/$2 .
  ln -s ../../raw/$3/*/$f1s.xml .
  ln -s ../../raw/$3/*/$f1s.tiff .
  ln -s ../../raw/$4 .
  cd ../..

  mkdir F2
  mkdir F2/raw F2/topo
  cd F2
  sed "s/.*threshold_geocode.*/threshold_geocode = 0/g" ../$5 | sed "s/.*threshold_snaphu.*/threshold_snaphu = 0/g" > $5
  cd topo
  ln -s ../../topo/dem.grd .
  cd ../raw
  ln -s ../topo/dem.grd .
  ln -s ../../raw/$1/*/$f2m.xml .
  ln -s ../../raw/$1/*/$f2m.tiff .
  ln -s ../../raw/$2 .
  ln -s ../../raw/$3/*/$f2s.xml .
  ln -s ../../raw/$3/*/$f2s.tiff .
  ln -s ../../raw/$4 .
  cd ../..

  mkdir F3
  mkdir F3/raw F3/topo
  cd F3
  sed "s/.*threshold_geocode.*/threshold_geocode = 0/g" ../$5 | sed "s/.*threshold_snaphu.*/threshold_snaphu = 0/g" > $5
  cd topo
  ln -s ../../topo/dem.grd .
  cd ../raw
  ln -s ../topo/dem.grd .
  ln -s ../../raw/$1/*/$f3m.xml .
  ln -s ../../raw/$1/*/$f3m.tiff .
  ln -s ../../raw/$2 .
  ln -s ../../raw/$3/*/$f3s.xml .
  ln -s ../../raw/$3/*/$f3s.tiff .
  ln -s ../../raw/$4 .
  cd ../..

# 
# process data
# 
  if ($seq == 1) then
    cd F1/raw
    align_tops.csh $f1m $2 $f1s $4 dem.grd
    set mpre1 = `echo $f1m | awk '{ print "S1A"substr($1,16,8)"_"substr($1,25,6)"_F"substr($1,7,1)}'`
    set spre1 = `echo $f1s | awk '{ print "S1A"substr($1,16,8)"_"substr($1,25,6)"_F"substr($1,7,1)}'`
    cd ../../F2/raw
    align_tops.csh $f2m $2 $f2s $4 dem.grd
    set mpre2 = `echo $f2m | awk '{ print "S1A"substr($1,16,8)"_"substr($1,25,6)"_F"substr($1,7,1)}'`
    set spre2 = `echo $f2s | awk '{ print "S1A"substr($1,16,8)"_"substr($1,25,6)"_F"substr($1,7,1)}'`
    cd ../../F3/raw
    align_tops.csh $f3m $2 $f3s $4 dem.grd
    set mpre3 = `echo $f3m | awk '{ print "S1A"substr($1,16,8)"_"substr($1,25,6)"_F"substr($1,7,1)}'`
    set spre3 = `echo $f3s | awk '{ print "S1A"substr($1,16,8)"_"substr($1,25,6)"_F"substr($1,7,1)}'`
    cd ../../F1
    p2p_S1A_TOPS.csh $mpre1 $spre1 $5
    cd ../F2
    p2p_S1A_TOPS.csh $mpre2 $spre2 $5
    cd ../F3
    p2p_S1A_TOPS.csh $mpre3 $spre3 $5
  else if ($seq == 2) then
    cd F1/raw
    align_tops.csh $f1m $2 $f1s $4 dem.grd >& log &
    set mpre1 = `echo $f1m | awk '{ print "S1A"substr($1,16,8)"_"substr($1,25,6)"_F"substr($1,7,1)}'`
    set spre1 = `echo $f1s | awk '{ print "S1A"substr($1,16,8)"_"substr($1,25,6)"_F"substr($1,7,1)}'`
    cd ../../F2/raw
    align_tops.csh $f2m $2 $f2s $4 dem.grd >& log &
    set mpre2 = `echo $f2m | awk '{ print "S1A"substr($1,16,8)"_"substr($1,25,6)"_F"substr($1,7,1)}'`
    set spre2 = `echo $f2s | awk '{ print "S1A"substr($1,16,8)"_"substr($1,25,6)"_F"substr($1,7,1)}'`
    cd ../../F3/raw
    align_tops.csh $f3m $2 $f3s $4 dem.grd >& log &
    set mpre3 = `echo $f3m | awk '{ print "S1A"substr($1,16,8)"_"substr($1,25,6)"_F"substr($1,7,1)}'`
    set spre3 = `echo $f3s | awk '{ print "S1A"substr($1,16,8)"_"substr($1,25,6)"_F"substr($1,7,1)}'`
    wait
    cd ../../F1
    p2p_S1A_TOPS.csh $mpre1 $spre1 $5 >&log&
    cd ../F2
    p2p_S1A_TOPS.csh $mpre2 $spre2 $5 >&log&
    cd ../F3
    p2p_S1A_TOPS.csh $mpre3 $spre3 $5 >&log&
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

  mkdir merge
  cd merge
  ln -s ../$5 .
  ln -s ../topo/dem.grd .
  if (-f tmp.filelist) then
    rm tmp.filelist
  endif
  set pth1 = `ls ../F1/intf/*/*PRM | awk NR==1'{print $1}' | awk -F"/" '{for (i=1;i<NF;i++) printf("%s/",$i)}'`
  set prm1m = `ls ../F1/intf/*/*PRM | awk NR==1'{print $1}' | awk -F"/" '{print $NF}'`
  set prm1s = `ls ../F1/intf/*/*PRM | awk NR==2'{print $1}' | awk -F"/" '{print $NF}'`
  echo $pth1":"$prm1m":"$prm1s > tmp.filelist
  
  set pth2 = `ls ../F2/intf/*/*PRM | awk NR==1'{print $1}' | awk -F"/" '{for (i=1;i<NF;i++) printf("%s/",$i)}'`
  set prm2m = `ls ../F2/intf/*/*PRM | awk NR==1'{print $1}' | awk -F"/" '{print $NF}'`
  set prm2s = `ls ../F2/intf/*/*PRM | awk NR==2'{print $1}' | awk -F"/" '{print $NF}'`
  echo $pth2":"$prm2m":"$prm2s >> tmp.filelist

  set pth3 = `ls ../F3/intf/*/*PRM | awk NR==1'{print $1}' | awk -F"/" '{for (i=1;i<NF;i++) printf("%s/",$i)}'`
  set prm3m = `ls ../F3/intf/*/*PRM | awk NR==1'{print $1}' | awk -F"/" '{print $NF}'`
  set prm3s = `ls ../F3/intf/*/*PRM | awk NR==2'{print $1}' | awk -F"/" '{print $NF}'`
  echo $pth3":"$prm3m":"$prm3s >> tmp.filelist

  merge_unwrap_geocode_tops.csh tmp.filelist $5



    

