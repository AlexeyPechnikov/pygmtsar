#!/bin/csh -f 
#       $Id$

# compute the mean and standard deviations of grid files

# input: a list of corr.grd 
# output: 


alias rm 'rm -f'


if ($#argv != 2) then
  echo ""
  echo "Usage: stack_coherence_mask.csh grid_list average_corr_threshold "
  echo ""
  echo "  stack the coherence and create a mask_def.grd "
  echo ""
  echo "  note that the grid of the grd files must be consistent" 
  echo ""
  exit 1
endif

set list = $1 
set threshold = $2

set file0 = `head -1 $list`
cp $file0 ./mask_def.grd
gmt grdmath mask_def.grd 0 MUL = mask_def.grd
set ct = `wc -l $list | awk '{print $1}'`
foreach line (`cat $list`)
  gmt grdmath $line mask_def.grd ADD = mask_def.grd
end

gmt grdmath mask_def.grd $ct DIV = mask_def.grd
gmt grdmath mask_def.grd $threshold GE 0 NAN = mask_def.grd

