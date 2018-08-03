#/bin/csh -f
# 
# By Xiaohua XU, 03/12/2018
#
# Unwrap interferograms parallelly using GNU parallel
#

if ($#argv != 1) then
  echo ""
  echo "Usage: unwrap_parallel.csh Ncores"
  echo ""
  echo "    Run unwrapping jobs parallelly. Need to install GNU parallel first."
  echo "    Note, run this in the intf_all folder where all the interferograms are stored. "
  echo ""
  exit
endif

set ncores = $1
set d1 = `date`

foreach line (`awk '{print $0}' $1`)
  echo "unwrap_intf.csh $line > log_$line.txt" >> unwrap.cmd
end

parallel --jobs $ncores < unwrap.cmd

echo ""
echo "Finished all unwrapping jobs..."
echo ""

set d2 = `date`

#echo "parallel --jobs $ncores < intf_tops.cmd" | mail -s "Unwrapping finished" "balabala@gmail.com" 
