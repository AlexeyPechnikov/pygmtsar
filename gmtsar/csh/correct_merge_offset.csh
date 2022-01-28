#!/bin/csh -f

if ($#argv != 4) then
  echo ""
  echo "Usage: correct_merge_offset.csh merge_list merge_log input.grd output.grd"
  echo ""
  echo "   used to correct offset between subswaths after merging (e.g. caused by range time diff)"
  echo ""
  exit 1
endif

set merge_list = $1
set merge_log = $2
set input = $3
set output = $4
set wid = 5
set spc = 5

set nf = `wc -l $merge_list | awk '{print $1}'`

set grd1 = `head -1 $merge_list | awk -F: '{print $2}'`
if ($nf == 3) then
  set grd2 = `head -2 $merge_list | tail -1 | awk -F: '{print $2}'`
  set grd3 = `tail -1 $merge_list | awk -F: '{print $2}'`
else if ($nf == 2) then
  set grd2 = `tail -1 $merge_list | awk -F: '{print $2}'`
else
  echo "incorrect number of files found in merge_list"
  exit 1
endif
set incx = `gmt grdinfo $input -C | awk '{print $8}'`
set xmin = `gmt grdinfo $input -C | awk '{print $2}'`
set xmax = `gmt grdinfo $input -C | awk '{print $3}'`
set ymin = `gmt grdinfo $input -C | awk '{print $4}'`
set ymax = `gmt grdinfo $input -C | awk '{print $5}'`

#echo $grd1 $grd2 $grd3
set nx1 = `gmt grdinfo $grd1 -C | awk '{print $10}'`
set n1 = `grep n1 merge_log | awk '{print $NF}'`
set ovl12 = `grep ovl merge_log | awk -F: '{print $2}' | awk -F, '{print $1}'`
set stitch_position1 = `echo $nx1 $n1 $ovl12 | awk '{print $1+$2-$3}'`
set position1 = `echo $incx $stitch_position1 | awk '{print $1*$2}'`

if ($nf == 3) then
  set nx2 = `gmt grdinfo $grd2 -C | awk '{print $10}'`
  set n2 = `grep n2 merge_log | awk '{print $NF}'`
  set ovl23 = `grep ovl merge_log | awk -F: '{print $2}' | awk -F, '{print $2}'`
  set stitch_position2 = `echo $nx1 $nx2 $n2 $ovl12 $ovl23 | awk '{print $1+$2-1-$4-$5+$3}'`
  set position2 = `echo $incx $stitch_position2 | awk '{print $1*$2}'`
endif

echo "Stitch positions $stitch_position1 $stitch_position2 ..."

set R1 = `echo $incx $stitch_position1 $ymin $ymax | awk '{print "-R"$1*($2-'$spc'-'$wid')"/"$1*($2-'$spc')"/"$3"/"$4}'`
set R2 = `echo $incx $stitch_position1 $ymin $ymax | awk '{print "-R"$1*($2+1+'$spc')"/"$1*($2+1+'$spc'+'$wid')"/"$3"/"$4}'`
gmt grdcut $input $R1 -G"tmp1_"$output
gmt grdcut $input $R2 -G"tmp2_"$output
gmt grdedit "tmp2_"$output -R"tmp1_"$output -G"tmp2_"$output
gmt grdmath "tmp2_"$output "tmp1_"$output SUB = tmp12_diff.grd
set diff1 = `gmt grdinfo tmp12_diff.grd -L1 -C | awk '{print $12}'`

if ($nf == 3) then
  set R1 = `echo $incx $stitch_position2 $ymin $ymax | awk '{print "-R"$1*($2-'$spc'-'$wid')"/"$1*($2-'$spc')"/"$3"/"$4}'`
  set R2 = `echo $incx $stitch_position2 $ymin $ymax | awk '{print "-R"$1*($2+1+'$spc')"/"$1*($2+1+'$spc'+'$wid')"/"$3"/"$4}'`
  gmt grdcut $input $R1 -G"tmp1_"$output
  gmt grdcut $input $R2 -G"tmp2_"$output
  gmt grdedit "tmp2_"$output -R"tmp1_"$output -G"tmp2_"$output
  gmt grdmath "tmp2_"$output "tmp1_"$output SUB = tmp23_diff.grd
  set diff2 = `gmt grdinfo tmp23_diff.grd -L1 -C | awk '{print $12}'`
endif

if ($nf == 2) then
  gmt grdcut $input -R$xmin/$position1/$ymin/$ymax -G"tmp1_"$output
  gmt grdcut $input -R$position1/$xmax/$ymin/$ymax -G"tmp2_"$output
  gmt grdmath "tmp2_"$output $diff1 SUB = "tmp2_"$output
  gmt grdpaste "tmp1_"$output "tmp2_"$output -G$output
  rm "tmp1_"$output "tmp2_"$output tmp12_diff.grd
else if ($nf == 3) then
  gmt grdcut $input -R$xmin/$position1/$ymin/$ymax -G"tmp1_"$output
  gmt grdcut $input -R$position1/$position2/$ymin/$ymax -G"tmp2_"$output
  gmt grdcut $input -R$position2/$xmax/$ymin/$ymax -G"tmp3_"$output
  gmt grdmath "tmp2_"$output $diff1 SUB = "tmp2_"$output
  gmt grdmath "tmp3_"$output $diff1 SUB $diff2 SUB = "tmp3_"$output
  gmt grdpaste "tmp1_"$output "tmp2_"$output -G"tmp4_"$output
  gmt grdpaste "tmp4_"$output "tmp3_"$output -G$output
  rm "tmp1_"$output "tmp2_"$output "tmp3_"$output "tmp4_"$output tmp12_diff.grd tmp23_diff.grd
endif




