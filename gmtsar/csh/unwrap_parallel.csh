#!/bin/csh -f
# 
# By Xiaohua XU, 03/12/2018
#
# Unwrap interferograms parallelly using GNU parallel
#
# IMPORTANT: put a script called unwrap_intf.csh in the current folder or in binary PATH
# e.g. 
#   cd $1
#   snaphu[_interp].csh 0.1 0 
#   cd ..
#

if ($#argv != 3) then
  echo ""
  echo "Usage: unwrap_parallel.csh unwrap_script intflist Ncores"
  echo ""
  echo "    Run unwrapping jobs parallelly. Need to install GNU parallel first."
  echo "    Note, run this in the intf_all folder where all the interferograms are stored. "
  echo ""
  exit
endif

set unwrapscript = $1
set intflist = $2
set ncores = $3
set d1 = `date`

foreach line (`awk '{print $0}' $intflist`)
  echo "$unwrapscript $line > log_$line.txt" >> unwrap.cmd
end
# long option --jobs is not supported on MacOSX
parallel -j $ncores < unwrap.cmd

echo ""
echo "Finished all unwrapping jobs..."
echo ""

set d2 = `date`
# long option --jobs is not supported on MacOSX
#echo "parallel -j $ncores < intf_tops.cmd" | mail -s "Unwrapping finished" "balabala@gmail.com" 
