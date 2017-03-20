#!/bin/csh -f
#       $Id$
#
# Xiaohua(Eric) Xu, Mar 18, 2017

#
# used for creating a new frame based on two input frames.
#

  if ($#argv != 5 && $#argv != 4) then
    echo ""
    echo "Usage: create_frame_tops.csh ****1.SAFE ****2.SAFE ****.EOF two_pins.llt [mode]"
    echo "       create one frame based on the two input frames, precise/restituted orbit is required"
    echo ""
    echo "  format of two_pins.llt (lon2/lat2 comes later than lon1/lat1 in orbit time):"
    echo "    lon1 lat1"
    echo "    lon2 lat2"
    echo ""
    echo "  outputs:"
    echo "    new.SAFE --> datetime1-datetime2.SAFE"
    echo ""
    echo "  Note:"
    echo "    The two input .SAFE file should be in time order, if the two_pin.llt is at wrong location, the program will output all the bursts."
    echo "    mode = 1, output vv; mode = 2, output vh; default is vv"
    echo ""
    exit 1
  endif

  if ($#argv == 4) then
    set mode = `echo "vv"`
  else 
    if ($5 == 1) then
      set mode = `echo "vv"`
    else
      set mode = `echo "vh"`
    endif
  endif

  echo "Combining $mode data..."

  set pth = `pwd`
  set file1 = $1
  set file2 = $2
  set orb = $3
  set tps = $4

  if (-d new.SAFE) then
    rm -r new.SAFE
  endif

  mkdir new.SAFE
  mkdir new.SAFE/annotation new.SAFE/measurement

# work on the first subswath
  cd $file1/annotation
  set f1 = `ls *iw1*$mode*xml | awk '{print substr($1,1,length($1)-4)}'`
  cd $pth
  cd $file2/annotation
  set f2 = `ls *iw1*$mode*xml | awk '{print substr($1,1,length($1)-4)}'`
  cd $pth

  cd new.SAFE
  ln -s ../$file1/annotation/$f1.xml .
  ln -s ../$file1/measurement/$f1.tiff .
  ln -s ../$file2/annotation/$f2.xml .
  ln -s ../$file2/measurement/$f2.tiff .

  make_s1a_tops $f1.xml $f1.tiff tmp1 0
  ext_orb_s1a tmp1.PRM ../$orb tmp1
  
  set ll1 = `awk NR==1'{print $0}' ../$tps`
  set ll2 = `awk NR==2'{print $0}' ../$tps`

  set azi1 = `echo "$ll1 0" | SAT_llt2rat tmp1.PRM 1 | awk '{print $2}'`
  set azi2 = `echo "$ll2 0" | SAT_llt2rat tmp1.PRM 1 | awk '{print $2}'`

  echo "Working on bursts covering $azi1 - $azi2 ..."
  assemble_tops $azi1 $azi2 $f1 $f2 new
  
  set tail1 = `echo $f1 | awk '{print substr($1,length($1)-16,17)}'`
  set t1 = `grep startTime new.xml | awk -F">" '{print $2}' | awk -F"<" '{print substr($1,1,4)substr($1,6,2)substr($1,9,2)"t"substr($1,12,2)substr($1,15,2)substr($1,18,2)}'`
  set t2 = `grep stopTime new.xml | awk -F">" '{print $2}' | awk -F"<" '{print substr($1,1,4)substr($1,6,2)substr($1,9,2)"t"substr($1,12,2)substr($1,15,2)substr($1,18,2)}'`
  mv new.xml annotation/"s1a-iw1-slc-vv-"$t1"-"$t2"-"$tail1.xml
  mv new.tiff measurement/"s1a-iw1-slc-vv-"$t1"-"$t2"-"$tail1.tiff
  rm *.tiff *.xml tmp1*
  cd ..

# wrok on the second subswath
  cd $file1/annotation
  set f1 = `ls *iw2*$mode*xml | awk '{print substr($1,1,length($1)-4)}'`
  cd $pth
  cd $file2/annotation
  set f2 = `ls *iw2*$mode*xml | awk '{print substr($1,1,length($1)-4)}'`
  cd $pth

  cd new.SAFE
  ln -s ../$file1/annotation/$f1.xml .
  ln -s ../$file1/measurement/$f1.tiff .
  ln -s ../$file2/annotation/$f2.xml .
  ln -s ../$file2/measurement/$f2.tiff .

  make_s1a_tops $f1.xml $f1.tiff tmp1 0
  ext_orb_s1a tmp1.PRM ../$orb tmp1
  
  set ll1 = `awk NR==1'{print $0}' ../$tps`
  set ll2 = `awk NR==2'{print $0}' ../$tps`

  set azi1 = `echo "$ll1 0" | SAT_llt2rat tmp1.PRM 1 | awk '{print $2}'`
  set azi2 = `echo "$ll2 0" | SAT_llt2rat tmp1.PRM 1 | awk '{print $2}'`

  echo "Working on bursts covering $azi1 - $azi2 ..."
  assemble_tops $azi1 $azi2 $f1 $f2 new 
  
  set tail1 = `echo $f1 | awk '{print substr($1,length($1)-16,17)}'`
  set t1 = `grep startTime new.xml | awk -F">" '{print $2}' | awk -F"<" '{print substr($1,1,4)substr($1,6,2)substr($1,9,2)"t"substr($1,12,2)substr($1,15,2)substr($1,18,2)}'`
  set t2 = `grep stopTime new.xml | awk -F">" '{print $2}' | awk -F"<" '{print substr($1,1,4)substr($1,6,2)substr($1,9,2)"t"substr($1,12,2)substr($1,15,2)substr($1,18,2)}'`
  mv new.xml annotation/"s1a-iw2-slc-vv-"$t1"-"$t2"-"$tail1.xml
  mv new.tiff measurement/"s1a-iw2-slc-vv-"$t1"-"$t2"-"$tail1.tiff
  rm *.tiff *.xml tmp1*
  cd ..

# wrok on the third subswath
  cd $file1/annotation
  set f1 = `ls *iw3*$mode*xml | awk '{print substr($1,1,length($1)-4)}'`
  cd $pth
  cd $file2/annotation
  set f2 = `ls *iw3*$mode*xml | awk '{print substr($1,1,length($1)-4)}'`
  cd $pth

  cd new.SAFE
  ln -s ../$file1/annotation/$f1.xml .
  ln -s ../$file1/measurement/$f1.tiff .
  ln -s ../$file2/annotation/$f2.xml .
  ln -s ../$file2/measurement/$f2.tiff .

  make_s1a_tops $f1.xml $f1.tiff tmp1 0
  ext_orb_s1a tmp1.PRM ../$orb tmp1

  set ll1 = `awk NR==1'{print $0}' ../$tps`
  set ll2 = `awk NR==2'{print $0}' ../$tps`

  set azi1 = `echo "$ll1 0" | SAT_llt2rat tmp1.PRM 1 | awk '{print $2}'`
  set azi2 = `echo "$ll2 0" | SAT_llt2rat tmp1.PRM 1 | awk '{print $2}'`

  echo "Working on bursts covering $azi1 - $azi2 ..."
  assemble_tops $azi1 $azi2 $f1 $f2 new

  set tail1 = `echo $f1 | awk '{print substr($1,length($1)-16,17)}'`
  set t1 = `grep startTime new.xml | awk -F">" '{print $2}' | awk -F"<" '{print substr($1,1,4)substr($1,6,2)substr($1,9,2)"t"substr($1,12,2)substr($1,15,2)substr($1,18,2)}'`
  set t2 = `grep stopTime new.xml | awk -F">" '{print $2}' | awk -F"<" '{print substr($1,1,4)substr($1,6,2)substr($1,9,2)"t"substr($1,12,2)substr($1,15,2)substr($1,18,2)}'`
  mv new.xml annotation/"s1a-iw3-slc-vv-"$t1"-"$t2"-"$tail1.xml
  mv new.tiff measurement/"s1a-iw3-slc-vv-"$t1"-"$t2"-"$tail1.tiff
  rm *.tiff *.xml tmp1*
  cd ..

  cp $file1/manifest.safe new.SAFE/

# edit the name of the new .SAFE file 
  set tail2 = `echo $file1 | awk '{print substr($1,length($1)-22,23)}'`
  set t1 = `ls new.SAFE/annotation/*.xml | awk '{print substr($1,length($1)-43,6)}' | gmt gmtinfo -C | awk '{print $1}'`
  set t2 = `ls new.SAFE/annotation/*.xml | awk '{print substr($1,length($1)-27,6)}' | gmt gmtinfo -C | awk '{print $2}'`
  set date1 = `ls new.SAFE/annotation/*.xml | awk '{print substr($1,length($1)-52,8)}' | gmt gmtinfo -C | awk '{print $1}'`

  mv new.SAFE "S1A_IW_SLC__1SSV_"$date1"T"$t1"_"$date1"T"$t2"_"$tail2




