#!/bin/csh -f 
#       $Id$

# compute the mean and standard deviations of grid files

# Creator: Xiaopeng Tong, David Sandwell 
# Date:    July 6th 2013

# input: a list of the gird files 
# output: the mean and standard deviations of the grids 


alias rm 'rm -f'


if ($#argv != 4) then
  echo ""
  echo "Usage: stack.csh  grid.list scale mean.grd std.grd "
  echo ""
  echo "  compute the mean and standard deviations of grid files "
  echo ""
  echo "  grid.list   --  a list of grd file names"
  echo "  scale       --  a scale factor put 1 if not scale"
  echo "  mean.grd    --  output file name of the mean grid "
  echo "  std.grd     --  output file name of the standard deviation grid "
  echo ""
  echo "  note that the grid of the grd files must be consistent" 
  echo ""
  exit 1
endif

if (! -e $1) then
  echo ""
  echo "no input file found: $1"
  echo ""
  exit 1
endif

set list = $1 
set scale = $2
set outmean = $3
set outstd = $4


# compute the mean 
echo "computing the mean of the grids .."
@ num = 1
foreach name (`cat $list`)
  if (! -e $name) then
    echo " Error: file not found: $name "
    echo ""
    exit 1
  endif
  if ($num == 1) then
    gmt grdmath $name = sum.grd 
  else 
    gmt grdmath $name sum.grd ADD = sumtmp.grd
    mv sumtmp.grd sum.grd
  endif
@ num ++
end  
@ num --
gmt grdmath sum.grd $num DIV = $outmean


# compute the standard deviation
echo "compute the standard deviation .. "
@ num = 1
foreach name (`cat $list`)
  if ($num == 1) then
    gmt grdmath $name $outmean SUB SQR = sum2.grd 
  else 
    gmt grdmath $name $outmean SUB SQR sum2.grd ADD = sum2tmp.grd
    mv sum2tmp.grd sum2.grd
  endif
@ num ++
end 
@ num --
gmt grdmath sum2.grd $num DIV SQRT = $outstd

# scale them
gmt grdmath $outmean $scale MUL = tmp.grd 
mv tmp.grd $outmean
gmt grdmath $outstd $scale MUL = tmp.grd 
mv tmp.grd $outstd

# clean up 
rm sum.grd sum2.grd

#
#  plot the results
#
foreach fname ($outmean $outstd)
    if ("$fname" == "$outmean") then
        set label = "Mean of Image Stack"
    else
        set label = "Std. Dev. of Image Stack"
    endif
    set name = `basename $fname .grd`
    gmt grdgradient $name.grd -Nt.9 -A0. -G$name.grad.grd
    set tmp = `gmt grdinfo -C -L2 $name.grd`
    set limitU = `echo $tmp | awk '{printf("%5.1f", $7)}'`
    set limitL = `echo $tmp | awk '{printf("%5.1f", $6)}'`
    set std = `echo $tmp | awk '{printf("%5.1f", $13)}'`
    gmt makecpt -Cseis -I -Z -T"$limitL"/"$limitU"/0.1 -D > $name.cpt
    set boundR = `gmt grdinfo $name.grd -C | awk '{print ($3-$2)/4}'`
    set boundA = `gmt grdinfo $name.grd -C | awk '{print ($5-$4)/4}'`
    gmt grdimage $name.grd -I$name.grad.grd -C$name.cpt -JX6.5i -Bxaf+lRange -Byaf+lAzimuth -BWSen -X1.3i -Y3i -P -K > $name.ps
    gmt psscale -R$name.grd -J -DJTC+w5/0.2+h+e -C$name.cpt -Bxaf+l"$label" -By -O >> $name.ps
    gmt psconvert -Tf -P -A -Z $name.ps
    echo "Mean of stack map: $name.pdf"
    rm $name.cpt $name.grad.grd
end
