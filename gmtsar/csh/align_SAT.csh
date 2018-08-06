#!/bin/tcsh -f
# align images using orbits and DEM
# by Kang Wang in June. 2018
#

if ($#argv != 2) then
   echo ""
   echo "Usage: align_SAT.csh data_type(RAW/SLC) align.in"
   echo ""
   exit 1
endif

if ($1 != RAW && $1 != SLC) then
  echo "data type must be RAW or SLC"
  exit 1
endif

if (! -e dem.grd) then
   echo "DEM file not found!"
   exit 1
endif

 echo "Downsample the DEM data"
 gmt grdfilter dem.grd -D3 -Fg2 -I12s -Ni -Gflt.grd
 gmt grd2xyz --FORMAT_FLOAT_OUT=%lf flt.grd -s > topo.llt

 set master = `awk '{if (NR==1) print $1}' $2`
 awk '{if (NR>1) print $1}' $2 > slave.list

 mv $master.PRM $master.PRM0

 echo "calculating the SAT height for master $master"
 calc_dop_orb $master.PRM0 $master.log 0
 cat $master.PRM0 $master.log > $master.PRM
 echo "fdd1                    = 0" >> $master.PRM
 echo "fddd1                   = 0" >> $master.PRM
 rm -f $master.log
 
 set earth_radius = `grep earth_radius $master.PRM|tail -1 | awk '{print $3}'`
 set fd1 = `grep fd1 $master.PRM |tail -1 |awk '{print $3}'`
 update_PRM $master.PRM earth_radius $earth_radius
 update_PRM $master.PRM fd1 $fd1

 SAT_llt2rat $master".PRM" 1 < topo.llt > master.ratll

  set rmax = `grep num_rng_bins $master".PRM" | awk '{print $3}'`
  set amax = `grep num_lines $master".PRM" | awk '{print $3}'`

if ($1 == RAW) then
 echo "Focusing the master - START"
  sarp.csh $master.PRM
 echo "Focusing the master - END"
endif

set nslave = `cat slave.list |wc -l`
set islave = 0
foreach slave (`cat slave.list`)
    @ islave = $islave + 1
    echo "working on " $slave  " [ $islave / $nslave ]"


    mv $slave.PRM $slave.PRM0
    calc_dop_orb $slave.PRM0 $slave.log $earth_radius

    cat $slave.PRM0 $slave.log > $slave.PRM
    echo "fdd1                    = 0" >> $slave.PRM
    echo "fddd1                   = 0" >> $slave.PRM
    
    update_PRM $slave.PRM earth_radius $earth_radius
    update_PRM $slave.PRM fd1 $fd1

    rm -f $slave.log

    if ($1 == RAW) then
     echo "Focusing the slave - START"
     sarp.csh $slave.PRM
     echo "Focusing the slave - END"
    endif

    SAT_llt2rat $slave".PRM" 1 < topo.llt > slave.ratll
    paste master.ratll slave.ratll | awk '{printf("%.6f %.6f %.6f %.6f %d\n", $1, $6 - $1, $2, $7 - $2, "100")}' > tmp.dat
    awk '{if($1 > 0 && $1 < '$rmax' && $3 > 0 && $3 < '$amax') print $0 }' < tmp.dat > offset.dat
    awk '{ printf("%.6f %.6f %.6f \n",$1,$3,$2) }' < offset.dat > r.xyz
    awk '{ printf("%.6f %.6f %.6f \n",$1,$3,$4) }' < offset.dat > a.xyz

   fitoffset.csh 3 3 offset.dat >> $slave".PRM"
   resamp $master.PRM $slave.PRM $slave.PRMresamp $slave.SLCresamp 4
   mv $slave.SLC $slave.SLC_old
   mv $slave.SLCresamp $slave.SLC
   cp $slave.PRMresamp $slave.PRM
end
