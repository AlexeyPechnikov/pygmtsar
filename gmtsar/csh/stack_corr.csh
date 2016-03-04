#!/bin/csh -f 
#       $Id$

# compute the mean correlation

# Creator: Xiaopeng Tong, David Sandwell 
# Date:    July 6th 2013

# input: a list of the correlation gird files 
# output: the mean correlation 

# Reference: Synthetic Aperture Radar Interferometry, Rosen et al., 2000 
if ($#argv != 2) then
  echo ""
  echo "Usage: stack_corr.csh list meancorr.grd"
  echo ""
  echo " compute the mean correlation from a stack of correlation grid files "
  echo ""
  echo "  list            --      a list of corr.grd file names" 
  echo "  meancorr.grd    --      output file name of the mean correlation grid "
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
set out = $2

# compute the mean 
echo "computing the mean correlation of the grids .."

@ num = 1
foreach cor (`cat $list`) 
  if (! -e $cor) then
    echo " Error: file not found: $cor "
    echo ""
    exit 1
  endif
  gmt grdmath $cor SQR = tmp.grd 
  if ($num == 1) then
    gmt grdmath 1 tmp.grd SUB tmp.grd DIV = sum.grd 
  else 
    gmt grdmath 1 tmp.grd SUB tmp.grd DIV sum.grd ADD = tmp2.grd  
    mv tmp2.grd sum.grd 
  endif
  @ num ++
end
@ num --
gmt grdmath 1 sum.grd $num DIV 1 ADD DIV SQRT = $out

# clean up
rm tmp.grd sum.grd

