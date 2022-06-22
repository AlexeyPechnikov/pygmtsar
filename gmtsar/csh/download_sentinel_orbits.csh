#!/bin/csh -f

# Written 04/05/2022 by Katherine Guns with aid from code snippets by Xiaohua Xu 
# and from ESA's website pages:
#   https://scihub.copernicus.eu/userguide/BatchScripting
#   https://scihub.copernicus.eu/twiki/do/view/SciHubUserGuide/ODataAPI#URI_Components
#   https://scihub.copernicus.eu/gnss/#/home

if ($#argv != 2) then
    echo ""
    echo "Usage: download_sentinel_orbits.csh safefilelist mode"
    echo "  Downloads precise or restituted orbits for specific Sentinel-1 *.SAFE data files  "
    echo ""
    echo "safefilelist:"
    echo "    absolutepathto/filename1.SAFE"
    echo "    absolutepathto/filename2.SAFE"
    echo "    ......"
    echo "mode:"
    echo "    mode 1 = precise orbits (POEORB)"
    echo "            (most users should choose precise orbits)"
    echo "    mode 2 = temporary (restituted) orbits (RESORB)"
    echo "            (only recent data (~last couple weeks) requires restituted"
    echo "            orbits, because precise orbits are not yet finalized)"
    echo ""
    echo "Example: download_sentinel_orbits.csh SAFEfile.list 1"
    echo ""
    echo "Note: "
    echo "  (1) Files listed in safefilelist should be the .SAFE directory with absolute path."
    echo ""
    exit 1
endif


#-------------------------
# PRECISE ORBITS (POEORB)
#-------------------------

if ($2 == 1) then
    echo " Downloading Precise Orbits (POEORB)..."
    #start working with SAFE file list
    foreach line (` awk -F"/" '{print $(NF)}' $1`)   #pull the name of the SAFE file from end of path
      set orbittype="AUX_POEORB"
      echo " "
      echo "--------------------- "
      echo "Finding orbits for ${line}..."
      set date1 = `echo $line | awk -F"_" '{print substr($6,1,8)}' `                
      set SAT1 = ` echo $line | awk -F"_" '{print $1}' `                 

      # get the orbit file names 
      set n1 = `date -v-1d -jf "%Y%m%d" $date1 +%Y%m%d`
      set n2 = `date -v+1d -jf "%Y%m%d" $date1 +%Y%m%d`

      echo "Required orbit file dates: ${n1} to  ${n2}..." 
  
      # Format SAFEfile date constraints for ESA database query
      set startorbtime = ` echo $n1 | awk '{printf "%d-%s-%sT00:00:00.000Z",substr($1,1,4),substr($1,5,2),substr($1,7,2)}' `
      set endorbtime = ` echo $n2 | awk '{printf "%d-%s-%sT23:59:59.999Z",substr($1,1,4),substr($1,5,2),substr($1,7,2)}' `

      echo "Querying ESA POD Hub archive..."
      # Run the query
      wget --no-check-certificate --user={gnssguest} --password={gnssguest} --output-document=orbitquery.txt "https://scihub.copernicus.eu/gnss/search?q=beginPosition:[${startorbtime} TO ${endorbtime}] AND endPosition:[${startorbtime} TO ${endorbtime}] AND platformname:Sentinel-1 AND filename:${SAT1}_* AND producttype:${orbittype}"

      echo "Checking query for existing orbit file..."
      
      set orbit = ` grep "title" orbitquery.txt | tail -1 | awk '{printf "%s.EOF",substr($1,8,73)}' `
      set esaID = ` grep "uuid" orbitquery.txt | awk '{print substr($2,13,36)}' `           
      if (! -f $orbit) then
        if (${esaID} == "") then
          echo "Query Failed -- possible issues:"
          echo " (1) couldn't connect to ESA POD hub (check manually)"
          echo " (2) an orbit file for those dates does not exist yet"
          echo "     -- check the resistited orbit files             "
          exit 1
        else
          echo "Query successful -- downloading orbit file..."
          wget --content-disposition --continue --user={gnssguest} --password={gnssguest} "https://scihub.copernicus.eu/gnss/odata/v1/Products('${esaID}')/"`echo '$'`"value"
          echo "...orbit file ${orbit} downloaded"
        endif  
      else  
        echo "...orbit file already exists"
        echo " "
      endif
    end
    # clean up
    rm orbitquery.txt
endif

#----------------------------
# RESTITUTED ORBITS (RESORB)
#----------------------------

if ($2 == 2) then
    echo " Downloading temporary Restituted Orbits (RESORB)..."
    #start working with SAFE file list
    foreach line (` awk -F"/" '{print $(NF)}' $1`)   #pull the name of the SAFE file from end of path
      set orbittype="AUX_RESORB"
      echo " "
      echo "--------------------- "
      echo "Finding orbits for ${line}..."
      set date1 = `echo $line | awk -F"_" '{print substr($6,1,8)}' `                
      set datetime1 = `echo $line | awk -F"_" '{printf "%s-%s-%s:%s:%s:%s",substr($6,1,4),substr($6,5,2),substr($6,7,2),substr($6,10,2),substr($6,12,2),substr($6,14,2)}' `                
      set datetime2 = `echo $line | awk -F"_" '{printf "%s-%s-%s:%s:%s:%s",substr($7,1,4),substr($7,5,2),substr($7,7,2),substr($7,10,2),substr($7,12,2),substr($7,14,2)}' `                
      set SAT1 = ` echo $line | awk -F"_" '{print $1}' `                 

      # get the orbit file names 
      set n1 = `date -v-3H -jf "%Y-%m-%d:%H:%M:%S" $datetime1 +"%Y-%m-%d %H:%M:%S" `
      set n2 = `date -v+3H -jf "%Y-%m-%d:%H:%M:%S" $datetime2 +"%Y-%m-%d %H:%M:%S" `

      echo "Required orbit file dates: ${n1} to  ${n2}..." 
  
      # Format SAFEfile date constraints for ESA database query
      set startorbtime = ` echo $n1 | awk '{printf "%sT%s.000Z",$1,$2}' `
      set endorbtime = ` echo $n2 | awk '{printf "%sT%s.000Z",$1,$2}' ` 
      
      echo "Querying ESA POD Hub archive..."
      # Run the query
      wget --no-check-certificate --user={gnssguest} --password={gnssguest} --output-document=orbitquery.txt "https://scihub.copernicus.eu/gnss/search?q=beginPosition:[${startorbtime} TO ${endorbtime}] AND endPosition:[${startorbtime} TO ${endorbtime}] AND platformname:Sentinel-1 AND filename:${SAT1}_* AND producttype:${orbittype}"

      echo "Checking query for existing restituted orbit file..."
      
      set orbit = ` grep "title" orbitquery.txt | tail -1 | awk '{printf "%s.EOF",substr($1,8,73)}' `
      set esaID = ` grep "uuid" orbitquery.txt | awk '{print substr($2,13,36)}' `           
      if (! -f $orbit) then
        if (${esaID} == "") then
          echo "Query Failed -- possible issues:"
          echo " (1) couldn't connect to ESA POD hub (check manually)"
          echo " (2) an orbit file for those dates does not exist yet"
          echo " "
          exit 1
        else
          echo "Query successful -- downloading restituted orbit file..."
          wget --content-disposition --continue --user={gnssguest} --password={gnssguest} "https://scihub.copernicus.eu/gnss/odata/v1/Products('${esaID}')/"`echo '$'`"value"
          echo "...restituted orbit file ${orbit} downloaded"
        endif  
      else  
        echo "...restituted orbit file already exists"
        echo " "
      endif
    end
    # clean up
    rm orbitquery.txt
endif
