#!/bin/csh -f
#       $Id$

# Select pairs according to the given threshold in time and baseline
# used for time series analysis

# Xiaohua(Eric) Xu, Jan 21 2016
#

  if ($#argv != 3) then
    echo ""
    echo "Usage: select_pairs.csh baseline_table.dat threshold_time threshold_baseline"
    echo "  generate the input file for intf_tops.csh with given threshold of time and baseline"
    echo ""
    echo "  outputs:"
    echo "    intf.in"
    echo ""
    exit 1
  endif

  set file = $1
  set dt = `echo $2 | awk '{print $0}'`
  set db = `echo $3 | awk '{printf $0}'`

# loop over possible pairs
  rm intf.in  

  awk '{print 2014+$3/365.25, $5, $1}' < $1 > text
  set region = `gmt gmtinfo text -C | awk '{print $1-0.5, $2+0.5, $3-50, $4+50}'`
  gmt pstext text -JX8.8i/6.8i -R$region[1]/$region[2]/$region[3]/$region[4] -D0.2/0.2 -X1.5i -Y1i -K -N -F+f8,Helvetica+j5 > baseline.ps  

  foreach line1 (`awk '{print $1":"$2":"$3":"$4":"$5}' < $file`)
    foreach line2 (`awk '{print $1":"$2":"$3":"$4":"$5}' < $file`)
      set t1 = `echo $line1 | awk -F: '{print $3}'`
      set t2 = `echo $line2 | awk -F: '{print $3}'`
      set b1 = `echo $line1 | awk -F: '{print $5}'`
      set b2 = `echo $line2 | awk -F: '{print $5}'`
      set n1 = `echo $line1 | awk -F: '{print $1}'`
      set n2 = `echo $line2 | awk -F: '{print $1}'`
      #if ($t1 < $t2 && $t2 - $t1 < $dt && $db0 < $db) then
      if ($t1 < $t2 & $t2 - $t1 < $dt) then
        set db0 = `echo $b1 $b2 | awk '{printf "%d", sqrt(($1-$2)*($1-$2))}'`
        if ($db0 < $db) then
          echo $n1 $n2 | awk '{print $1":"$2}' >> intf.in
          echo $t1 $b1 | awk '{print $1/365.25+2014, $2}' >> tmp
          echo $t2 $b2 | awk '{print $1/365.25+2014, $2}' >> tmp
          gmt psxy tmp -R -J -K -O >> baseline.ps
          rm tmp
        endif
      endif
    end
  end

  awk '{print $1,$2}' < text > text2
  gmt psxy text2 -Sp0.2c -G0 -R -JX -Ba0.5:"year":/a50g00f25:"baseline (m)":WSen -O >> baseline.ps
  
