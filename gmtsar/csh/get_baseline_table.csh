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

awk '{printf("%.6f %.6f %s\n",$2%1000.0/365.25+int($2/1000.0),$5,$1)}' baseline_table.dat > tmp_text
set region = `gmt gmtinfo tmp_text -C | awk '{print $1-0.5, $2+0.5, $3-50, $4+50}'`
gmt pstext tmp_text -JX8.8i/6.8i -R$region[1]/$region[2]/$region[3]/$region[4] -D0.2/0.2 -X1.5i -Y1i -K -N -F+f8,Helvetica+j5 > baseline.ps
awk '{print $1,$2}' < tmp_text > tmp_text2
gmt psxy tmp_text2 -Sp0.2c -G0 -R -JX -Ba0.5:"year":/a50g00f25:"baseline (m)":WSen -O >> baseline.ps
rm tmp_text tmp_text2
gmt psconvert baseline.ps -Tf -A
rm baseline.ps

