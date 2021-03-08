#!/bin/csh -f
#       $Id$
#  Xiaohua Xu, Dec 17, 2018
#
#  resample a focused SLC image to a certain prf and rng_samp_rate
#
#

if ($#argv != 3) then
  echo ""
  echo "Usage: samp_slc.csh image_stem new_prf new_rng_samp_rate "
  echo ""
  echo "Example: samp_slc.csh IMG-ALPSRP022200660-H1.1__A 2000 3200000"
  echo ""
  echo "Note: A very small change in prf or rng_samp_rate change will be ignored"
  echo "      as they can be accounted during normal processing with resamp."
  echo "      Put 0 in new_prf or new_rng_samp_rate to igore resamping."
  echo ""
  exit 1
endif


set image = $1
set nprf = $2
set nrsr = $3

set oprf = `grep PRF $image.PRM | awk '{print $3}'`
set orsr = `grep rng_samp_rate $image.PRM | awk '{print $3}'`

set trsr = `echo $nrsr $orsr | awk '{printf("%.10f", $1/$2)}'`
set tprf = `echo $nprf $oprf | awk '{printf("%.10f", $1/$2)}'`

cp $image.PRM tmp_master.PRM
cp $image.PRM tmp_aligned.PRM

set ta = `echo $tprf | awk '{if($1>1) print 1;else if($1<1) print 2;else print 0;}'`
set tr = `echo $trsr | awk '{if($1>1) print 1;else if($1<1) print 2;else print 0;}'`

if ($ta == 2 && $nprf != 0) echo "Downsampling along azimuth is not recommended, may cause aliasing ..."
if ($tr == 2 && $nrsr != 0) echo "Downsampling along range   is not recommended, may cause aliasing ..."

if ($nrsr == 0 && $nprf == 0) then
  echo "specify at least one of new_prf/new_rng_samp_rate to be non-zero to run the script" 
  exit 1
endif

if ($nrsr != 0) then
  set r = `echo $trsr | awk '{printf("%.10f", -(1 - 1/$1))}'`
  set tmp = `grep num_rng_bins tmp_master.PRM | awk '{print $3}'`
  set new_num_rng_bins = `echo $tmp $trsr | awk '{printf("%d",$1*$2/4)}' | awk '{printf("%d",$1*4)}'`
  update_PRM tmp_aligned.PRM stretch_r $r
  update_PRM tmp_master.PRM rng_samp_rate $nrsr
  update_PRM tmp_master.PRM num_rng_bins $new_num_rng_bins
  set bytes = `echo $new_num_rng_bins | awk '{printf("%d",$1*4)}'`
  update_PRM tmp_master.PRM bytes_per_line $bytes
  update_PRM tmp_master.PRM good_bytes_per_line $bytes
endif

if ($nprf != 0) then
  set a = `echo $tprf | awk '{printf("%.10f", -(1 - 1/$1))}'`
  set tmp = `grep num_lines tmp_master.PRM | awk '{print $3}'`
  set new_nl = `echo $tmp $tprf | awk '{printf("%d",$1*$2/4)}' | awk '{printf("%d",$1*4)}'`
  update_PRM tmp_aligned.PRM a_stretch_a $a
  update_PRM tmp_master.PRM PRF $nprf
  update_PRM tmp_master.PRM num_lines $new_nl
  update_PRM tmp_master.PRM num_valid_az $new_nl
  update_PRM tmp_master.PRM num_patches 1
  update_PRM tmp_master.PRM nrows $new_nl
endif

resamp tmp_master.PRM tmp_aligned.PRM tmp_aligned.PRMresamp tmp_aligned.SLCresamp 4

mv tmp_master.PRM $image.PRM
mv tmp_aligned.SLCresamp $image.SLC
rm tmp*


