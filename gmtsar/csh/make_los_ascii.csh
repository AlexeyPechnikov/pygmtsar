#!/bin/csh -f
# subsample  and LOS grid and add the look vector 
# need to include the PRM and LED file in the local directory
#
  if ($#argv <1 || $#argv >5) then
	echo ""
	echo "Usage: make_los_ascii.csh los.grd dem.grd -I0.01/0.01 PRM_file Satellite "
	echo ""
	echo "Example: make_los_ascii.csh los_ll.grd dem.grd -I0.01/0.01 IMG-HH-ALOS2046743050-150405-WBDR1.1__D-F1.PRM ALOS"
	echo ""
	echo "output: los.lltnde"
	exit 1
  endif

if ( -f ~/.quiet ) then
    set V = ""
else
	set V = "-V"
endif

gmt  grdsample $2 -Gtmp_topo.grd `grdinfo $1 -I-` `grdinfo $1 -I` -F
gmt  grdmath $1 0 MUL 1 ADD tmp_topo.grd MUL = tmp_topo.grd
gmt  grd2xyz $1 > tmp.xyz
gmt  grd2xyz tmp_topo.grd > tmp_topo.xyz
gmt  blockmedian tmp.xyz `grdinfo $1 -I-` $3 $V > tmp_b.xyz
gmt  blockmedian tmp_topo.xyz `grdinfo $1 -I-` $3 $V > tmp_topo_b.xyz

  set SAT = $5
  if($SAT == ERS) then
 	set SAT = 'ENVI'
  endif

  $SAT"_"look $4 < tmp_topo_b.xyz > tmp_topo_b_n.lltn
  awk '{ a=$3;getline <"tmp_topo_b_n.lltn";print $1,$2,$3,$4,$5,$6,a,-1}' tmp_b.xyz > los.lltnde
  rm tmp*

