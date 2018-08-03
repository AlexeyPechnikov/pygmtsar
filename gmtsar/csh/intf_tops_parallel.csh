#!/bin/csh -f

if ($#argv != 3) then
  echo ""
  echo "Usage: intf_tops_batch.csh intf.in batch.config Ncores"
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

parallel --jobs $ncores < intf_tops.cmd

set t2 = `date`
set dir0 = `pwd`
echo "Job started on $t1 and finishe don $t2 at $dir0 "|mail -s "TOPS intfs job finished" sddyxxh@gmail.com

