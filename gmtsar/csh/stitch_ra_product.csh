#!/bin/csh -f
#       $Id$
#
  if ($#argv != 4 && $#argv != 5 ) then
    echo " "
    echo "Usage: stitch_ra_product.csh PRM1 grid1 PRM2 grid2 [compute_trans]" 
    echo " "
    echo "  Note: set compute_trans to 1 to re-compute projection matric"
    echo "    dem.grd is required and need to cover the two grids."
    echo ""
    echo "  Product: saved in stitch_product.grd"
    echo " "
    exit 1
  endif 

  set buf = 500

  if ($#argv == 5) then
    set compute_trans = $5
  else
    set compute_trans = 0
  endif

  set prm1 = $1
  set grid1 = $2
  set prm2 = $3
  set grid2 = $4

  set fs1 = `grep rng_samp_rate $prm1 | tail -1 | awk '{printf("%d", $3)}'`
  set fs2 = `grep rng_samp_rate $prm2 | tail -1 | awk '{printf("%d", $3)}'`
  set near1 = `grep near_range $prm1 | awk '{printf("%d", $3)}'`
  set near2 = `grep near_range $prm2 | awk '{printf("%d", $3)}'`
  set t_start1 = `grep clock_start $prm1 | grep -v SC_clock_start | awk '{print $3}'`
  set t1 = `echo $t_start1 | awk '{printf("%d",$1*86400.0)}'`
  set t_start2 = `grep clock_start $prm2 | grep -v SC_clock_start | awk '{print $3}'`
  set t2 = `echo $t_start2 | awk '{printf("%d",$1*86400.0)}'`
  set nl1 = `grep num_lines $prm1 | awk '{print $3}'`
  set nl2 = `grep num_lines $prm2 | awk '{print $3}'`
  set xinc1 = `gmt grdinfo -C $grid1 | awk '{print $8}'`
  set xinc2 = `gmt grdinfo -C $grid2 | awk '{print $8}'`
  set nx1 = `gmt grdinfo -C $grid1 | awk '{print $10}'`
  set nx2 = `gmt grdinfo -C $grid2 | awk '{print $10}'`
  set yinc1 = `gmt grdinfo -C $grid1 | awk '{print $9}'`
  set yinc2 = `gmt grdinfo -C $grid2 | awk '{print $9}'`
  set ny1 = `gmt grdinfo -C $grid1 | awk '{print $11}'`
  set ny2 = `gmt grdinfo -C $grid2 | awk '{print $11}'`
  set prf1 = `grep PRF $prm1 | awk '{print $3}'`
  set prf2 = `grep PRF $prm2 | awk '{print $3}'`

  set r1x1 = `gmt grdinfo -C $grid1 | awk '{print $2}'`
  set r1x2 = `gmt grdinfo -C $grid1 | awk '{print $3}'`
  set r1y1 = `gmt grdinfo -C $grid1 | awk '{print $4}'`
  set r1y2 = `gmt grdinfo -C $grid1 | awk '{print $5}'`
  set r2x1 = `gmt grdinfo -C $grid2 | awk '{print $2}'`
  set r2x2 = `gmt grdinfo -C $grid2 | awk '{print $3}'`
  set r2y1 = `gmt grdinfo -C $grid2 | awk '{print $4}'`
  set r2y2 = `gmt grdinfo -C $grid2 | awk '{print $5}'`

echo $fs1
 
  set rng_pixel = `echo $fs1 | awk '{printf("%.12f",299792458.0/$1/2.0)}'` 
  echo "range pixel size is $rng_pixel ..."
  
  if ($fs1 != $fs2) then
    echo "[ERROR]: The two grids have different range sampling rate ..."
    exit 1
  endif

  if ($xinc1 != $xinc2) then
    echo "[ERROR]: The two grids have different range pixel spcaing ..."
    exit 1
  endif
  if ($yinc1 != $yinc2) then
    echo "[ERROR]: The two grids have different azimuth pixel spcaing ..."
    exit 1
  endif

  echo "parameters of grid 1:" $near1 $t_start1 $nl1 $xinc1 $yinc1
  echo "parameters of grid 2:" $near2 $t_start2 $nl2 $xinc2 $yinc2

  if ($t1 >= $t2) then
    echo "[ERROR]: Time of grid2 is ahead of time of grid 1 ..."
    exit 1
  endif

  set dy_pixel = `echo $t_start1 $t_start2 $yinc1 $prf1 | awk '{printf ("%d",($2-$1)*86400.0*$4/$3)}'`
  echo "shift in azimuth is $dy_pixel pixels ..."
  set ny_ovlp = `echo $ny1 $dy_pixel | awk '{printf("%d",$1-$2)}'`

  if ($near1 <= $near2) then
    set dx_pixel = `echo $near1 $near2 $xinc1 $rng_pixel | awk '{printf ("%d",($2-$1)/$3/$4)}'`    
    echo "shift in near range is $dx_pixel pixels ..."
    set rtmp = `echo $xinc2 $dx_pixel | awk '{print $1*$2}'`
    gmt grdcut $grid2 -R0/$rtmp/0/$r2y2 -Gtmp.grd
    gmt grdmath tmp.grd 0 MUL 0 NAN = tmp.grd
    gmt grdedit tmp.grd -R-$rtmp/0/0/$r2y2 -Gtmp2.grd
    gmt grdpaste $grid2 tmp2.grd -Gtmp_grid2.grd
    set rtmp2 = `echo $rtmp $r2x2 | awk '{print $1+$2}'`
    gmt grdedit tmp_grid2.grd -R0/$rtmp2/0/$r2y2
    
    if ($rtmp2 >= $r1x2) then
      set dx_pixel2 = `echo $rtmp2 $r1x2 $xinc1 | awk '{printf ("%d",($1-$2)/$3)}'`
      echo "extending first grid with $dx_pixel2 pixels ..."
      set rtmp = `echo $xinc1 $dx_pixel2 | awk '{print $1*$2}'`
      gmt grdcut $grid1 -R0/$rtmp/0/$r1y2 -Gtmp.grd
      gmt grdmath tmp.grd 0 MUL 0 NAN = tmp.grd
      set rtmp2 = `echo $rtmp $r1x2 | awk '{print $1+$2}'`
      gmt grdedit tmp.grd -R$r1x2/$rtmp2/0/$r1y2 -Gtmp2.grd
      gmt grdpaste $grid1 tmp2.grd -Gtmp_grid1.grd
      set rxx = `echo "0/$rtmp2"`
      set r_end = `echo $rtmp2`
    else
      set dx_pixel2 = `echo $rtmp2 $r1x2 $xinc1 | awk '{printf ("%d",($2-$1)/$3)}'`
      echo "extending second grid with $dx_pixel2 pixels ..."
      set rtmp = `echo $xinc1 $dx_pixel2 | awk '{print $1*$2}'`
      gmt grdcut $grid2 -R0/$rtmp/0/$r2y2 -Gtmp.grd
      gmt grdmath tmp.grd 0 MUL 0 NAN = tmp.grd
      set rtmp3 = `echo $rtmp $rtmp2 | awk '{print $1+$2}'`
      gmt grdedit tmp.grd -R$rtmp2/$rtmp3/0/$r2y2 -Gtmp2.grd
      gmt grdpaste tmp2.grd tmp_grid2.grd -Gtmp3.grd
      mv tmp3.grd tmp_grid2.grd
      set rxx = `echo "0/$rtmp3"`
      set r_end = `echo $rtmp3`
    endif
    cp $prm1 stitch_product.PRM
    set LED = `grep led_file $prm1 | awk '{print $3}'`
    set pth = `echo $prm1 | awk -F"/" '{for (i=1;i<NF;i++) printf("%s/",$i)}'`
    cp $pth/$LED .
  else
    set dx_pixel = `echo $near1 $near2 $xinc1 $rng_pixel | awk '{printf ("%d",($1-$2)/$3/$4)}'`
    echo "shift in near range is $dx_pixel pixels ..."
    set rtmp = `echo $xinc1 $dx_pixel | awk '{print $1*$2}'`
    gmt grdcut $grid1 -R0/$rtmp/0/$r1y2 -Gtmp.grd
    gmt grdmath tmp.grd 0 MUL 0 NAN = tmp.grd
    gmt grdedit tmp.grd -R-$rtmp/0/0/$r1y2 -Gtmp2.grd
    gmt grdpaste $grid1 tmp2.grd -Gtmp_grid1.grd
    set rtmp2 = `echo $rtmp $r1x2 | awk '{print $1+$2}'`
    gmt grdedit tmp_grid1.grd -R0/$rtmp2/0/$r1y2

    if ($rtmp2 >= $r2x2) then
      set dx_pixel2 = `echo $rtmp2 $r2x2 $xinc2 | awk '{printf ("%d",($1-$2)/$3)}'`
      echo "extending second grid with $dx_pixel2 pixels ..."
      set rtmp = `echo $xinc2 $dx_pixel2 | awk '{print $1*$2}'`
      gmt grdcut $grid2 -R0/$rtmp/0/$r2y2 -Gtmp.grd
      gmt grdmath tmp.grd 0 MUL 0 NAN = tmp.grd
      set rtmp2 = `echo $rtmp $r2x2 | awk '{print $1+$2}'`
      gmt grdedit tmp.grd -R$r2x2/$rtmp2/0/$r2y2 -Gtmp2.grd
      gmt grdpaste $grid2 tmp2.grd -Gtmp_grid2.grd
      set rxx = `echo "0/$rtmp2"`
      set r_end = `echo $rtmp2`
    else
      set dx_pixel2 = `echo $rtmp2 $r2x2 $xinc2 | awk '{printf ("%d",($2-$1)/$3)}'`
      echo "extending first grid with $dx_pixel2 pixels ..."
      set rtmp = `echo $xinc1 $dx_pixel2 | awk '{print $1*$2}'`
      gmt grdcut $grid1 -R0/$rtmp/0/$r1y2 -Gtmp.grd
      gmt grdmath tmp.grd 0 MUL 0 NAN = tmp.grd 
      set rtmp3 = `echo $rtmp $rtmp2 | awk '{print $1+$2}'`
      gmt grdedit tmp.grd -R$rtmp2/$rtmp3/0/$r1y2 -Gtmp2.grd
      gmt grdpaste tmp2.grd tmp_grid1.grd -Gtmp3.grd
      mv tmp3.grd tmp_grid1.grd
      set rxx = `echo "0/$rtmp3"`
      set r_end = `echo $rtmp3`
    endif
    cp $prm2 stitch_product.PRM
    set LED = `grep led_file $prm2 | awk '{print $3}'`
    set pth = `echo $prm2 | awk -F"/" '{for (i=1;i<NF;i++) printf("%s/",$i)}'`
    cp $pth/$LED .
  endif

  set y1 = `echo $ny1 $ny_ovlp $yinc1 | awk '{printf("%d",$3*int($1-$2/2))}'`
  gmt grdcut tmp_grid1.grd -R$rxx/0/$y1 -Gtmp1.grd
  echo "Cutting first grid to $rxx/0/$y1 ..."

  set y2 = `echo $ny_ovlp $yinc2| awk '{printf("%d",$2*int($1/2))}'`
  gmt grdcut tmp_grid2.grd -R$rxx/$y2/$r2y2 -Gtmp2.grd

  set y2_start = `gmt grdinfo -C tmp1.grd | awk '{print $5}'`
  set ny2_new = `gmt grdinfo -C tmp2.grd | awk '{print $11}'`
  set y2_end = `echo $y2_start $ny2_new $yinc2 | awk '{printf("%d",($2*$3+$1))}'`
  gmt grdedit tmp2.grd -R$rxx/$y2_start/$y2_end -Gtmp3.grd
  echo "Cutting second grid to $rxx/$y2_start/$y2_end ..."

  gmt grdpaste tmp1.grd tmp3.grd -Gstitch_product.grd

# update PRM file
  update_PRM stitch_product.PRM num_rng_bins $r_end
  update_PRM stitch_product.PRM num_lines $y2_end
  update_PRM stitch_product.PRM nrows $y2_end
  update_PRM stitch_product.PRM num_valid_az $y2_end
  update_PRM stitch_product.PRM clock_start $t_start1
  set t_tmp = `grep SC_clock_start $prm1 | awk '{print $3}'`
  update_PRM stitch_product.PRM SC_clock_start $t_tmp
  set bytes = `echo $r_end | awk '{print $1*4}'`
  update_PRM stitch_product.PRM bytes_per_line $bytes
  update_PRM stitch_product.PRM good_bytes_per_line $bytes

  rm tmp*  

  if ($compute_trans == 1) then
    gmt grd2xyz --FORMAT_FLOAT_OUT=%lf dem.grd -s | SAT_llt2rat stitch_product.PRM 1 -bod  > trans.dat
  endif


