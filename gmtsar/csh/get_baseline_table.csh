#!/bin/csh -f

if ($#argv != 2) then
  echo ""
  echo "Usage: get_baseline_table.csh prmlist master_prm"
  echo ""
  echo "     product: baseline_table.dat"
  echo ""
  exit 1
endif

set list = $1
set master = $2

if (-f baseline_table.dat) rm baseline_table.dat

foreach prm (`cat $list`)
  baseline_table.csh $master $prm >> baseline_table.dat 
end
