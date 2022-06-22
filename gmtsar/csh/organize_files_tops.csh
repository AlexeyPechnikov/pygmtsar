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
# EDIT 04/2022: Replaced the wgetasf calls with ESA api commands

  if ($#argv != 3 && $#argv != 4) then
    echo ""
    echo "Usage: organize_files_tops.csh filelist pins.ll mode [observation_mode]"
    echo "  organize one track of S1A TOPS data, redefine frames, auto-download precise orbits"
    echo "    if restituted orbits are required, user must download separately"
    echo ""
    echo "filest:"
    echo "    pth_filename1"
    echo "    pth_filename2"
    echo "    ......"
    echo ""
    echo "pins.ll:"
    echo "    lon11 lat11 [lon12 lat12] [lon13 lat13]"
    echo "    lon21 lat21 [lon22 lat22] [lon23 lat23]"
    echo "    ......"
    echo ""
    echo "Note: "
    echo "    files listed in filelist should be the .SAFE directory with absolute path."
    echo "    mode = 1 will tell how many records are gonna be generated. mode = 2 will do the organizing."
    echo ""
    exit 1
  endif

  # set local orbit directory
  set orb_dir = "Sentinel_Orbits"

  set ii = 0
  set mode = $3
  set org_mod = "$org_mod"
  if ($#argv == 4) then
    set org_mod = $4
  endif

  if (-f tmprecord) rm tmprecord

  # divide the list of files into sets, and create frames based on the given pins
  foreach line (`awk '{print $0}' $1`)
    set file1 = `echo $line | awk -F"," '{print $1}'`
    set date1 = `echo $file1 | awk '{print substr($1,length($1)-54,8)}'`
    set SAT1 = `echo $file1 | awk '{print substr($1,length($1)-71,3)}'`
    set orbittype="AUX_POEORB" #assume precise orbits
    
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
          set t1 = `date -jf "%Y%m%dT%H%M%S" $tt +%s`
          if ($t1 > $t2) set jj = 0
          set tt = `echo $line2 | awk '{print substr($1,length($1)-38,15)}'`
          set t2 = `date -jf "%Y%m%dT%H%M%S" $tt +%s`
        end
        echo "" | awk '{printf("%s\n","...")}'

        # get the orbit file names and download
        set n1 = `date -v-1d -jf "%Y%m%d" $date0 +%Y%m%d`
        set n2 = `date -v+1d -jf "%Y%m%d" $date0 +%Y%m%d`
        #set orb = `grep $SAT0 orbits.list | grep $n1 | grep $n2 | tail -1`
        
        echo "Required orbit file dates: ${n1} to  ${n2}..."

        # Format SAFEfile date constraints for ESA database query
        set startorbtime = ` echo $n1 | awk '{printf "%d-%s-%sT00:00:00.000Z",substr($1,1,4),substr($1,5,2),substr($1,7,2)}' `
        set endorbtime = ` echo $n2 | awk '{printf "%d-%s-%sT23:59:59.999Z",substr($1,1,4),substr($1,5,2),substr($1,7,2)}' `

        #echo "Querying ESA POD Hub orbit archive..."
        # Run the query
        wget --no-check-certificate --user={gnssguest} --password={gnssguest} --output-document=orbitquery.txt "https://scihub.copernicus.eu/gnss/search?q=beginPosition:[${startorbtime} TO ${endorbtime}] AND endPosition:[${startorbtime} TO ${endorbtime}] AND platformname:Sentinel-1 AND filename:${SAT0}_* AND producttype:${orbittype}"

        echo "Checking query for existing orbit file..."

        set orbit = ` grep "title" orbitquery.txt | tail -1 | awk '{printf "%s.EOF",substr($1,8,73)}' `
        set esaID = ` grep "uuid" orbitquery.txt | awk '{print substr($2,13,36)}' `           
        if (! -f $orbit) then
          # IF esaID is empty that means no uuid was found corresponding to the date window 
          if (${esaID} == "") then
            echo "Query Failed -- possible issues:"
            echo " - an orbit file for those dates may not exist yet"
            echo "     --> check the resistited orbit files (AUX_RESORB) "
            echo ""
            echo " SKIP $date0, as precise orbit file may not exist ..."
            echo ""
            echo $file1 > tmprecord
            set file0 = `echo $file1`
            set date0 = `echo $date1`
            set SAT0 = `echo $SAT1`
            continue
          else
            if (-f $orb_dir/$SAT0/$orbit) then
              ln -s $orb_dir/$SAT0/$orbit
            else
              echo "Query successful -- downloading orbit file..."
              wget --content-disposition --continue --user={gnssguest} --password={gnssguest} "https://scihub.copernicus.eu/gnss/odata/v1/Products('${esaID}')/"`echo '$'`"value"
              echo "...orbit file ${orbit} downloaded"
            endif
          endif  
        else  
          echo "...orbit file already exists..."
          echo " "
        endif

        # compute azimuth for the start and end 
        set pin1 = `head -1 $2 | awk '{print $1,$2}'` 
        set f1 = `head -1 tmprecord`
        make_s1a_tops $f1/annotation/*iw1*"$org_mod"*xml $f1/measurement/*iw1*"$org_mod"*tiff tmp2 0
        ext_orb_s1a tmp2.PRM $orbit tmp2
        set tmpazi = `echo $pin1 | awk '{print $1,$2,0}' | SAT_llt2rat tmp2.PRM 1 | awk '{printf("%d",$2+0.5)}'`
        # refine the calculation in case the pin is far away from the starting frame. (baseline error)
        shift_atime_PRM.csh tmp2.PRM $tmpazi
        set azi1 = `echo $pin1 | awk '{print $1,$2,0}' | SAT_llt2rat tmp2.PRM 1 | awk '{printf("%d",$2+0.5 + '$tmpazi')}'`
        
        set pin2 = `tail -1 $2 | awk '{print $1,$2}'`
        set f2 = `tail -1 tmprecord`
        make_s1a_tops $f2/annotation/*iw1*"$org_mod"*xml $f2/measurement/*iw1*"$org_mod"*tiff tmp2 0
        ext_orb_s1a tmp2.PRM $orbit tmp2
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
                create_frame_tops.csh tmprecord $orbit tmp1llt $org_mod
                set newfile = `ls -t -d *.SAFE | awk NR==1'{print $0}'`
                set Frame1 = `grep azimuthAnxTime $newfile/annotation/*iw1*"$org_mod"*xml | head -1 | awk -F">" '{print $2}' | awk -F"<" '{printf("F%.4d", $1+0.5)}'`
                set Frame2 = `grep azimuthAnxTime $newfile/annotation/*iw1*"$org_mod"*xml | tail -1 | awk -F">" '{print $2}' | awk -F"<" '{printf("F%.4d", $1+0.5)}'` 
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

  # process the last set of files
  echo "" | awk '{printf("%s ","Combing")}' 
  set jj = 1
  set t2 = 9999999999
  foreach line2 (`awk '{print $0}' tmprecord`)
    echo $line2 | awk '{printf("%s ",$1)}'
    set tt = `echo $line2 | awk '{print substr($1,length($1)-54,15)}'`
    set t1 = `date -jf "%Y%m%dT%H%M%S" $tt +%s`
    if ($t1 > $t2) set jj = 0
    set tt = `echo $line2 | awk '{print substr($1,length($1)-38,15)}'`
    set t2 = `date -jf "%Y%m%dT%H%M%S" $tt +%s`
  end
  echo "" | awk '{printf("%s\n","...")}'

  # get the orbit file names and download
  set n1 = `date -v-1d -jf "%Y%m%d" $date0 +%Y%m%d`
  set n2 = `date -v+1d -jf "%Y%m%d" $date0 +%Y%m%d`
  #set orb = `grep $SAT0 orbits.list | grep $n1 | grep $n2 | tail -1`

  echo "Required orbit file dates: ${n1} to  ${n2}..."

  # Format SAFEfile date constraints for ESA database query
  set startorbtime = ` echo $n1 | awk '{printf "%d-%s-%sT00:00:00.000Z",substr($1,1,4),substr($1,5,2),substr($1,7,2)}' `
  set endorbtime = ` echo $n2 | awk '{printf "%d-%s-%sT23:59:59.999Z",substr($1,1,4),substr($1,5,2),substr($1,7,2)}' `

  # Run the query
  wget --no-check-certificate --user={gnssguest} --password={gnssguest} --output-document=orbitquery.txt "https://scihub.copernicus.eu/gnss/search?q=beginPosition:[${startorbtime} TO ${endorbtime}] AND endPosition:[${startorbtime} TO ${endorbtime}] AND platformname:Sentinel-1 AND filename:${SAT0}_* AND producttype:${orbittype}"

  echo "Checking query for existing orbit file..."
      
  set orbit = ` grep "title" orbitquery.txt | tail -1 | awk '{printf "%s.EOF",substr($1,8,73)}' `
  set esaID = ` grep "uuid" orbitquery.txt | awk '{print substr($2,13,36)}' `           
  if (! -f $orbit) then
    # IF esaID is empty that means no uuid was found corresponding to the date window 
    if (${esaID} == "") then
      echo "Query Failed -- possible issues:"
      echo " - an orbit file for those dates may not exist yet"
      echo "     --> check the resistited orbit files (AUX_RESORB)"
      echo ""
      echo " SKIP $date0, as precise orbit file may not exist ..."
      echo ""
      exit 1  
    else
      if (-f $orb_dir/$SAT0/$orbit) then
        ln -s $orb_dir/$SAT0/$orbit
      else
        echo "Query successful -- downloading orbit file..."
        wget --content-disposition --continue --user={gnssguest} --password={gnssguest} "https://scihub.copernicus.eu/gnss/odata/v1/Products('${esaID}')/"`echo '$'`"value"
        echo "...orbit file ${orbit} downloaded"
      endif
    endif  
  else  
    echo "...orbit file already exists..."
    echo " "
  endif


  # check the start and the end, make sure the start comes later than the first line of the first file and the end comes before the last line of the last file
  set pin1 = `head -1 $2 | awk '{print $1,$2}'` 
  set f1 = `head -1 tmprecord`
  make_s1a_tops $f1/annotation/*iw1*"$org_mod"*xml $f1/measurement/*iw1*"$org_mod"*tiff tmp2 0
  ext_orb_s1a tmp2.PRM $orbit tmp2
  set tmpazi = `echo $pin1 | awk '{print $1,$2,0}' | SAT_llt2rat tmp2.PRM 1 | awk '{printf("%d",$2+0.5)}'`
  shift_atime_PRM.csh tmp2.PRM $tmpazi
  set azi1 = `echo $pin1 | awk '{print $1,$2,0}' | SAT_llt2rat tmp2.PRM 1 | awk '{printf("%d",$2+0.5 + '$tmpazi')}'`

  set pin2 = `tail -1 $2 | awk '{print $1,$2}'` 
  set f2 = `tail -1 tmprecord`
  make_s1a_tops $f2/annotation/*iw1*"$org_mod"*xml $f2/measurement/*iw1*"$org_mod"*tiff tmp2 0
  ext_orb_s1a tmp2.PRM $orbit tmp2
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
          create_frame_tops.csh tmprecord $orbit tmp1llt $org_mod
          set newfile = `ls -t -d *.SAFE | awk NR==1'{print $0}'`
          set Frame1 = `grep azimuthAnxTime $newfile/annotation/*iw1*"$org_mod"*xml | head -1 | awk -F">" '{print $2}' | awk -F"<" '{printf("F%.4d", $1+0.5)}'`
          set Frame2 = `grep azimuthAnxTime $newfile/annotation/*iw1*"$org_mod"*xml | tail -1 | awk -F">" '{print $2}' | awk -F"<" '{printf("F%.4d", $1+0.5)}'` 
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
  
