#!/bin/csh -f
#
# by Xiaohua XU, 03/06/2018
#

if ($#argv != 3) then
  echo ""
  echo "Usage: intf_tops_parallel.csh intf.in batch.config Ncores"
  echo ""
  echo "    Run intf jobs parallelly. Need to install GNU parallel first."
  echo "    e.g. sudo port install parallel"
  echo ""
  exit
endif


rm -f intf_tops.cmd

set t1 = `date`
set ncores = $3
foreach intf (`awk '{print $0}' $1`)
   set date1 =  `echo $intf |awk -F":" '{print $1}'|cut -c 4-11`
   set date2 =  `echo $intf |awk -F":" '{print $2}'|cut -c 4-11`
   set logfile = "intf_"$date1"_"$date2".log"
   set infile = "intf_"$date1"_"$date2".in"
   echo $intf > $infile
   echo "intf_tops.csh $infile $2 > $logfile" >> intf_tops.cmd
end
# long option --jobs is not supported on MacOSX
parallel -j $ncores < intf_tops.cmd

set t2 = `date`
set dir0 = `pwd`
#echo "Job started on $t1 and finishe don $t2 at $dir0 "|mail -s "TOPS intfs job finished" balabala@gmail.com

