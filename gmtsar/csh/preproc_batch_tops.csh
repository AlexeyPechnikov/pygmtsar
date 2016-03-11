#!/bin/csh -f
#       $Id$

# Preprocess and align tops image stacks
# used for time series analysis

# Xiaohua(Eric) Xu, Jan 19 2016
#

  if ($#argv != 3) then
    echo ""
    echo "Usage: preproc_tops_batch.csh data.in dem.grd mode"
    echo "  preprocess and align a set of tops images in data.in, precise orbits required"
    echo ""
    echo "  format of data.in:"
    echo "    image_name:orbit_name"
    echo ""
    echo "  example of data.in"
    echo "    s1a-iw1-slc-vv-20150602t134416-20150602t134444-006195-00813b-001:S1A_OPER_AUX_POEORB_V20150601_20150603.EOF"
    echo "    s1a-iw1-slc-vv-20150626t134418-20150626t134446-006545-008b4a-001:S1A_OPER_AUX_POEORB_V20150625_20150627.EOF"
    echo ""
    echo "  outputs:"
    echo "    baseline.ps align_table.ra (contains info for precise geomatric alignment)"
    echo "    *.PRM *.LED *.SLC(mode 2)"
    echo ""
    exit 1
  endif

  set mode = $3
  set sl = `echo "1"`
   
  if ($mode == 2) then
    gmt grdfilter $2 -D2 -Fg2 -I12s -Gflt.grd
    gmt grd2xyz --FORMAT_FLOAT_OUT=%lf flt.grd -s > topo.llt
  endif

  # first one is master
  set master = `awk 'NR==1 {print $0}' $1 | awk '{print substr($1,16,8)}'` 
  rm baseline_table.dat *.PRM *.SLC *.LED

  # loop over all the files
  foreach line (`awk '{print $0}' $1`)

    # get names
    set image = `echo $line | awk -F: '{print $1}'`
    set orbit = `echo $line | awk -F: '{print $2}'`
    set stem = `echo $image | awk '{print substr($1,16,8)}'`
    
    # generate prms and leds
    make_s1a_tops $image.xml $image.tiff S1A$stem 0 

    # extract precise orbit
    ext_orb_s1a S1A$stem.PRM $orbit S1A$stem

    # get the range and azimuth table
    if ($mode == 2) then
      SAT_llt2rat S1A$stem.PRM 1 < topo.llt > tmp1.dat

      if($sl == 1) then
        make_s1a_tops $image.xml $image.tiff S1A$stem 1 
        ext_orb_s1a S1A$stem.PRM $orbit S1A$stem
        mv tmp1.dat tmpm.dat
        set sl = `echo "2"`
      else
        paste tmpm.dat tmp1.dat | awk awk '{printf("%.6f %.6f %.6f %.6f %d\n", $6, $6-$1, $7, $7-$2, "100")}' > tmp.dat
        set rmax = `grep num_rng_bins S1A$stem.PRM | awk '{print $3}'`
        set amax = `grep num_lines S1A$stem.PRM | awk '{print $3}'`
        awk '{if($1 > 0 && $1 < '$rmax' && $3 > 0 && $3 < '$amax') print $0 }' < tmp.dat > offset.dat
        awk '{ printf("%.6f %.6f %.6f \n",$1,$3,$2) }' < offset.dat > r.xyz
	awk '{ printf("%.6f %.6f %.6f \n",$1,$3,$4) }' < offset.dat > a.xyz
        gmt blockmedian r.xyz -R0/$rmax/0/$amax -I8/4 -r > rtmp.xyz
        gmt blockmedian a.xyz -R0/$rmax/0/$amax -I8/4 -r > atmp.xyz
        gmt surface rtmp.xyz -R0/$rmax/0/$amax -I8/4 -Grtmp.grd -T0.5 -N1000  -r -V
        gmt surface atmp.xyz -R0/$rmax/0/$amax -I8/4 -Gatmp.grd -T0.5 -N1000  -r -V
        gmt grdmath rtmp.grd FLIPUD = r.grd
        gmt grdmath atmp.grd FLIPUD = a.grd
         
        make_s1a_tops $image.xml $image.tiff S1A$stem 1 r.grd a.grd 
        ext_orb_s1a S1A$stem.PRM $orbit S1A$stem
      
        # update PRM parameters
        
        # resample to match super master
        resamp S1A$master.PRM S1A$stem.PRM S1A$stem.PRMresamp S1A$stem.SLCresamp 1
        mv S1A$stem.SLCresamp S1A$stem.SLC
        cp S1A$stem.PRMresamp S1A$stem.PRM       
        
        fitoffset.csh 3 3 offset.dat >> S1A$stem.PRM
        rm tmp.PRM offset.dat 
      endif
    endif
 
    # get the baseline info
    cp S1A$stem.PRM junk1
    calc_dop_orb junk1 junk2 0 0
    cat junk1 junk2 > S1A$stem.PRM

    baseline_table.csh S1A$master.PRM S1A$stem.PRM >> baseline_table.dat
    baseline_table.csh S1A$master.PRM S1A$stem.PRM GMT >> table.gmt
    
    rm junk1 junk2    
  end
  
  # use data inside the frame
  set rmax = `grep num_rng_bins S1A$stem.PRM | awk '{print $3}'`
  set amax = `grep num_lines S1A$stem.PRM | awk '{print $3}'`
  if ($mode == 2) then
    awk '{if($1 > 0 && $1 < '$rmax' && $2 > 0 && $2 < '$amax') print $0 }' < tmp.dat > align_table.ra 
  endif

  # plot the baseline table  (for S1A)
  awk '{print 2014+$1/365.25,$2,$7}' < table.gmt > text
  set region = `gmt gmtinfo text -C | awk '{print $1-0.5, $2+0.5, $3-50, $4+50}'`
  gmt pstext text -JX8.8i/6.8i -R$region[1]/$region[2]/$region[3]/$region[4] -D0.2/0.2 -X1.5i -Y1i -K -N -F+f8,Helvetica+j5 > baseline.ps
  awk '{print $1,$2}' < text > text2
  gmt psxy text2 -Sp0.2c -G0 -R -JX -Ba0.5:"year":/a50g00f25:"baseline (m)":WSen -O >> baseline.ps

  # clean up junks
  if ($mode == 2) then
    rm tmp* topo.llt flt.grd r.grd a.grd atmp* rtmp* r.xyz a.xyz
  endif
  rm text text2 table.gmt
