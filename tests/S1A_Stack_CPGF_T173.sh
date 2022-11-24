#!/bin/sh

if [ "$#" -eq 0 ]
then
    echo "Default launch saves downloaded dataset for future usage"
    force=0
elif [ "$#" -eq 1 ]
then
    echo "Force launch removes downloaded dataset to save disk space"
    force=1
else
    echo "Incorrect number of arguments"
    exit 1
fi

wget -q --show-progress --progress=bar:force:noscroll -c http://topex.ucsd.edu/gmtsar/tar/S1A_Stack_CPGF_T173.tar.gz
rm -rf raw_orig topo
tar xvzf S1A_Stack_CPGF_T173.tar.gz -C .
# cleanup
if [ "$force" -eq 1 ]
then
    echo "Removing downloaded datasets to free disk space following 'force' command line argument'..."
    rm *.tar.gz
fi
echo "Checking disk space and downloaded files"
df -h
ls -lh
echo "Running Python test script..."
rm -fr *.jpg
python3 ./S1A_Stack_CPGF_T173.py && echo SUCCESS
