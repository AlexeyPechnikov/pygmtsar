#!/bin/csh -f
#       $Id$
#
# Xiaohua(Eric) Xu, Mar 20 2017
#
# set a directory with all downloaded orbits. e.g. orb_dir/S1A and orb_dir/S1B
#
# set a local directory that stores S1A and S1B orbits. e.g. orb_dir/S1A and orb_dir/S1B
#
# For linux users, modify the date command accordingly in order to let the script run.
# Some linux system has issues with the if statement comapring tt, t1, t2. This is not fixed yet
#
# alias wgetasf to 'wget --http-user=**** --http-password=****' in .cshrc or .tcshrc file
# requires having an account on with ASF
#
# modified from organize_files_tops.csh for linux uses. date command and some if statements are changed
#


  if ($#argv != 3) then
    echo ""
    echo "Usage: organize_files_tops_linux.csh filelist pins.ll mode"
    echo "  organize one track of S1A TOPS data, redefine frames, precise/restituted orbit is required"
    echo ""
    echo "filest:"
    echo "    pth_filename1"
    echo "    pth_filename2"
    echo "    ......"
    echo ""
    echo "pins.ll:"
    echo "    lon1 lat1"
    echo "    lon2 lat2"
    echo "    ......"
    echo ""
    echo "Note: "
    echo "    files listed in filelist should be the .SAFE directory with absolute path."
    echo "    mode = 1 will tell how many records are gonna be generated. mode = 2 will do the organizing."
    echo ""
    exit 1
  endif

source ~/.cshrc

  set url_root = "https://s1qc.asf.alaska.edu/aux_poeorb"
  set orb_dir = "/geosat2/InSAR_Processing/Sentinel_Orbits"
  if (! -f orbits.list) then
    wget $url_root -O orbits.html
    grep EOF orbits.html | awk -F'"' '{print $2}' > orbits.list
    rm orbits.html
  endif

  set ii = 0
  set mode = $3

  if (-f tmprecord) rm tmprecord

  # divide the list of files into sets, and create frames based on the given pins
  foreach line (`awk '{print $0}' $1`)
    set file1 = `echo $line | awk -F"," '{print $1}'`
    set date1 = `echo $file1 | awk '{print substr($1,length($1)-54,8)}'`
    set SAT1 = `echo $file1 | awk '{print substr($1,length($1)-71,3)}'`
    
    if ($ii == 0) then
      set file0 = `echo $file1`
      set date0 = `echo $date1`
      set SAT0 = `echo $SAT1`
      echo $file1 > tmprecord
      set ii = 1
    else
      # gather files from the same date
      if ($date1 == $date0 && $SAT1 == $SAT0) then
        echo $file1 >> tmprecord
      else

        echo "" | awk '{printf("%s ","Combing")}' 
        set jj = 1
        set t2 = 9999999999
        # examining whether the frames are consecutive
        foreach line2 (`awk '{print $0}' tmprecord`)
          echo $line2 | awk '{printf("%s ",$1)}'
          set tt = `echo $line2 | awk '{print substr($1,length($1)-54,15)}'`
          set ss2 = `echo $tt|awk '{print substr($1,1,4)"/"substr($1,5,2)"/"substr($1,7,2)" "substr($1,10,2)":"substr($1,12,2)":"substr($1,14,2)}'`
          set t1 = `date --date="$ss2" +%s`
          #set t1 = `date -jf "%Y%m%dT%H%M%S" $tt +%s`
          set test = `echo $t1 $t2 | awk '{if ($1 > $2) print 1; else print 0}'`
          if ($test == 1) set jj = 0
          set tt = `echo $line2 | awk '{print substr($1,length($1)-38,15)}'`
          set ss2 = `echo $tt|awk '{print substr($1,1,4)"/"substr($1,5,2)"/"substr($1,7,2)" "substr($1,10,2)":"substr($1,12,2)":"substr($1,14,2)}'`
          set t2 = `date --date="$ss2" +%s`
          #set t2 = `date -jf "%Y%m%dT%H%M%S" $tt +%s`
        end
        echo "" | awk '{printf("%s\n","...")}'

        # get the orbit file names and download
        #set n1 = `date -v-1d -jf "%Y%m%d" $date0 +%Y%m%d`
        #set n2 = `date -v+1d -jf "%Y%m%d" $date0 +%Y%m%d`
        set n1 = ` date --date="$date0 - 1 day" +%Y%m%d `
        set n2 = ` date --date="$date0 + 1 day" +%Y%m%d `

        set orb = `grep $SAT0 orbits.list | grep $n1 | grep $n2 | tail -1`

echo $n1 $n2 $orb

        if (! -f $orb) then
          if (-f $orb_dir/$SAT0/$orb) then
            ln -s $orb_dir/$SAT0/$orb .
          else
            wgetasf $url_root"/"$orb
          endif
        endif

        # compute azimuth for the start and end 
        set pin1 = `head -1 $2 | awk '{print $1,$2}'` 
        set f1 = `head -1 tmprecord`
        make_s1a_tops $f1/annotation/*iw1*vv*xml $f1/measurement/*iw1*vv*tiff tmp2 0
        ext_orb_s1a tmp2.PRM $orb tmp2
        set tmpazi = `echo $pin1 | awk '{print $1,$2,0}' | SAT_llt2rat tmp2.PRM 1 | awk '{printf("%d",$2+0.5)}'`
        # refinie the calculation in case the pin is far away from the starting frame. (baseline error)
        shift_atime_PRM.csh tmp2.PRM $tmpazi
        set azi1 = `echo $pin1 | awk '{print $1,$2,0}' | SAT_llt2rat tmp2.PRM 1 | awk '{printf("%d",$2+0.5 + '$tmpazi')}'`
        
        set pin2 = `tail -1 $2 | awk '{print $1,$2}'`
        set f2 = `tail -1 tmprecord`
        make_s1a_tops $f2/annotation/*iw1*vv*xml $f2/measurement/*iw1*vv*tiff tmp2 0
        ext_orb_s1a tmp2.PRM $orb tmp2
        set tmpazi = `echo $pin2 | awk '{print $1,$2,0}' | SAT_llt2rat tmp2.PRM 1 | awk '{printf("%d",$2+0.5)}'`
        # refinie the calculation in case the pin is far away from the starting frame.
        shift_atime_PRM.csh tmp2.PRM $tmpazi
        set azi2 = `echo $pin2 | awk '{print $1,$2,0}' | SAT_llt2rat tmp2.PRM 1 | awk '{printf("%d",$2+0.5 + '$tmpazi')}'`

        set nl = `grep num_lines tmp2.PRM | awk '{print $3}'`

        if ($azi1 > 0 && $azi2 < $nl && $jj != 0) then  
          awk '{print $1","$2}' $2 > tmpllt
          set pin0 = `awk NR==1'{print $0}' tmpllt`
          foreach line2 (`awk '{print $0}' tmpllt`)
            if ($line2 != $pin0) then
              echo $pin0 | awk -F"," '{print $1,$2}' > tmp1llt
              echo $line2 | awk -F"," '{print $1,$2}' >> tmp1llt
              if ($mode != 1) then
                create_frame_tops.csh tmprecord $orb tmp1llt 1
                set newfile = `ls -t -d *.SAFE | awk NR==1'{print $0}'`
                set Frame1 = `grep azimuthAnxTime $newfile/annotation/*iw1*vv*xml | head -1 | awk -F">" '{print $2}' | awk -F"<" '{printf("F%.4d", $1+0.5)}'`
                set Frame2 = `grep azimuthAnxTime $newfile/annotation/*iw1*vv*xml | tail -1 | awk -F">" '{print $2}' | awk -F"<" '{printf("F%.4d", $1+0.5)}'` 
                set dirname = `echo $Frame1"_"$Frame2`
                echo "Created Frame $Frame1 - $Frame2 ..."
                echo ""
                if (! -d $dirname) mkdir $dirname
                mv $newfile $dirname
              else
                echo ""
                echo "Frames on date $date0 will be re-organized..."
                echo ""
              endif
              set pin0 = `echo $line2`
            endif
          end 
        else
          if ($jj == 0) then
            echo ""
            echo "SKIP $date0, as it stopped observation in the middle ..."
            echo ""
          else
            echo ""
            echo "SKIP $date0, as it does not have enough scenes ..."
            echo ""
          endif
        endif

        echo $file1 > tmprecord
        set file0 = `echo $file1`
        set date0 = `echo $date1`
        set SAT0 = `echo $SAT1`
      endif

    endif
  end 

  # proces the last set of files
  echo "" | awk '{printf("%s ","Combing")}' 
  set jj = 1
  set t2 = 9999999999
  foreach line2 (`awk '{print $0}' tmprecord`)
    echo $line2 | awk '{printf("%s ",$1)}'
    set tt = `echo $line2 | awk '{print substr($1,length($1)-54,15)}'`
    set ss2 = `echo $tt|awk '{print substr($1,1,4)"/"substr($1,5,2)"/"substr($1,7,2)" "substr($1,10,2)":"substr($1,12,2)":"substr($1,14,2)}'`
    set t1 = `date --date="$ss2" +%s`

    #set t1 = `date -jf "%Y%m%dT%H%M%S" $tt +%s`
    set test = `echo $t1 $t2 | awk '{if ($1 > $2) print 1; else print 0}'`
    if ($test == 1) set jj = 0
    set tt = `echo $line2 | awk '{print substr($1,length($1)-38,15)}'`
    set ss2 = `echo $tt|awk '{print substr($1,1,4)"/"substr($1,5,2)"/"substr($1,7,2)" "substr($1,10,2)":"substr($1,12,2)":"substr($1,14,2)}'`
    set t2 = `date --date="$ss2" +%s`
    #set t2 = `date -jf "%Y%m%dT%H%M%S" $tt +%s`
  end
  echo "" | awk '{printf("%s\n","...")}'

  # get the orbit file names and download
  set n1 = ` date --date="$date0 - 1 day" +%Y%m%d `
  set n2 = ` date --date="$date0 + 1 day" +%Y%m%d `
  #set n1 = `date -v-1d -jf "%Y%m%d" $date0 +%Y%m%d`
  #set n2 = `date -v+1d -jf "%Y%m%d" $date0 +%Y%m%d`
  set orb = `grep $SAT0 orbits.list | grep $n1 | grep $n2 | tail -1`

  if (! -f $orb) then
    if (-f $orb_dir/$SAT0/$orb) then
      ln -s $orb_dir/$SAT0/$orb .
    else
      wgetasf $url_root"/"$orb
    endif
  endif

  # check the start and the end, make sure the start comes later than the first line of the first file and the end comes before the last line of the last file
  set pin1 = `head -1 $2 | awk '{print $1,$2}'` 
  set f1 = `head -1 tmprecord`
  make_s1a_tops $f1/annotation/*iw1*vv*xml $f1/measurement/*iw1*vv*tiff tmp2 0
  ext_orb_s1a tmp2.PRM $orb tmp2
  set tmpazi = `echo $pin1 | awk '{print $1,$2,0}' | SAT_llt2rat tmp2.PRM 1 | awk '{printf("%d",$2+0.5)}'`
  shift_atime_PRM.csh tmp2.PRM $tmpazi
  set azi1 = `echo $pin1 | awk '{print $1,$2,0}' | SAT_llt2rat tmp2.PRM 1 | awk '{printf("%d",$2+0.5 + '$tmpazi')}'`

  set pin2 = `tail -1 $2 | awk '{print $1,$2}'` 
  set f2 = `tail -1 tmprecord`
  make_s1a_tops $f2/annotation/*iw1*vv*xml $f2/measurement/*iw1*vv*tiff tmp2 0
  ext_orb_s1a tmp2.PRM $orb tmp2
  set tmpazi = `echo $pin2 | awk '{print $1,$2,0}' | SAT_llt2rat tmp2.PRM 1 | awk '{printf("%d",$2+0.5)}'`
  shift_atime_PRM.csh tmp2.PRM $tmpazi
  set azi2 = `echo $pin2 | awk '{print $1,$2,0}' | SAT_llt2rat tmp2.PRM 1 | awk '{printf("%d",$2+0.5 + '$tmpazi')}'`
  set nl = `grep num_lines tmp2.PRM | awk '{print $3}'`

  # do the assembling
  if ($azi1 >= 0 && $azi2 < $nl && $jj != 0) then  
    awk '{print $1","$2","$3","$4","$5","$6}' $2 > tmpllt
    set pin0 = `awk NR==1'{print $0}' tmpllt`
    foreach line2 (`awk '{print $0}' tmpllt`)
      if ($line2 != $pin0) then
        echo $pin0 | awk -F"," '{print $1,$2,$3,$4,$5,$6}' > tmp1llt
        echo $line2 | awk -F"," '{print $1,$2,$3,$4,$5,$6}' >> tmp1llt

        if ($mode != 1) then
          create_frame_tops.csh tmprecord $orb tmp1llt 1
          set newfile = `ls -t -d *.SAFE | awk NR==1'{print $0}'`
          set Frame1 = `grep azimuthAnxTime $newfile/annotation/*iw1*vv*xml | head -1 | awk -F">" '{print $2}' | awk -F"<" '{printf("F%.4d", $1+0.5)}'`
          set Frame2 = `grep azimuthAnxTime $newfile/annotation/*iw1*vv*xml | tail -1 | awk -F">" '{print $2}' | awk -F"<" '{printf("F%.4d", $1+0.5)}'` 
          set dirname = `echo $Frame1"_"$Frame2`
          echo "Created Frame $Frame1 - $Frame2 ..."
          echo ""
          if (! -d $dirname) mkdir $dirname
          mv $newfile $dirname
        else
          echo ""
          echo "Frames on date $date0 will be re-organized..."
          echo ""
        endif
        set pin0 = `echo $line2`
      endif
    end
  else 
    if ($jj == 0) then
      echo ""
      echo "SKIP $date0, as it stopped observation in the middle ..."
      echo ""
    else
      echo ""
      echo "SKIP $date0, as it does not have enough scenes ..."
      echo ""
    endif
  endif   

  rm tmp*
  #rm *.EOF
  
