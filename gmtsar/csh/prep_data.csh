#!/bin/csh -f
#
# By Xiaohua XU, Jun, 2015
#
# prepare input data.in for preproc_batch_tops.csh
#
# alias wgetasf to 'wget --http-user=**** --http-password=****' in .cshrc or .tcshrc file
# requires having an account on with ASF
#
# set a local directory that stores S1A and S1B orbits. e.g. orb_dir/S1A and orb_dir/S1B
#
set orb_dir = "/geosat2/InSAR_Processing/Sentinel_Orbits"

rm data.in
ls *.xml > text.dat
set mstem = `awk 'NR==1 {print $0}' text.dat | awk '{print substr($1,16,8)}'`
set mname = `awk 'NR==1 {print $0}' text.dat | awk '{print substr($1,1,64)}'`
set rec = 0

if (! -f orbits.list) then
  set url_root = "https://s1qc.asf.alaska.edu/aux_poeorb/"
  wget $url_root -O orbits.html
  grep EOF orbits.html | awk -F'"' '{print $2}' > orbits.list
  rm orbits.html
endif

echo "Looping over all the lines"
foreach line ( ` awk '{ print $0 }' < text.dat ` )
  set stem = `echo $line | awk '{print substr($1,16,8)}'`
  set name = `echo $line | awk '{print substr($1,1,64)}'`
  if ($stem != $mstem) then
    echo "Writing record $mstem"
    #set n1 = `echo $mstem | awk '{print $1-1}'`
    #set n2 = `echo $mstem | awk '{print $1+1}'`

    # uses date for time manipulation -ben
    set n1 = `date -v-1d -jf "%Y%m%d" $mstem +%Y%m%d`
    set n2 = `date -v+1d -jf "%Y%m%d" $mstem +%Y%m%d`
    set SAT = `echo $mname | awk '{print toupper(substr($1,1,3))}'`

    #cp ../../../../orbit/*$n1*$n2* .
    set orb = `grep $SAT orbits.list | grep $n1 | grep $n2 | tail -1`
    if (! -f $orb) then
      if (-f $orb_dir/$SAT/$orb) then 
        cp $orb_dir/$SAT/$orb .
      else
        wgetasf $url_root/$orb
      endif
    endif
    #if (! -f $orb) wget $url_root"/"$orb
    #set orb = `ls *$n1*$n2*`
    #set orb = `ls *$mstem*EOF`
    set tmp_rec = `echo $rec $orb | awk '{print $1":"$2}'`
    echo $tmp_rec >> data.in
    set rec = `echo $name`
    set mstem = `echo $stem`
    set mname = `echo $name`
  else
    if ($rec == 0) then 
      echo "Setting a new name for rec"
      set rec = `echo $name`
    else
      set tmp_rec = `echo $name $rec | awk '{print $2":"$1}'`
      set rec = `echo $tmp_rec`
    endif
  endif
end

set n1 = `date -v-1d -jf "%Y%m%d" $mstem +%Y%m%d`
set n2 = `date -v+1d -jf "%Y%m%d" $mstem +%Y%m%d`
set SAT = `echo $mname | awk '{print toupper(substr($1,1,3))}'`
echo "Writing record $mstem"
#set orb = `ls *$n1*$n2*`
set orb = `grep $SAT orbits.list | grep $n1 | grep $n2 | tail -1`
if (! -f $orb) then
  if (-f $orb_dir/$SAT/$orb) then
    cp $orb_dir/$SAT/$orb .
  else
    wgetasf $url_root/$orb
  endif
endif
#if (! -f $orb)  wget $url_root"/"$orb
#set orb = `ls *$mstem*EOF`
set tmp_rec = `echo $rec $orb | awk '{print $1":"$2}'`
echo $tmp_rec >> data.in
rm text.dat 
rm orbits.list
