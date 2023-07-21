#!/bin/sh

#wget -q --show-progress --progress=bar:force:noscroll -c http://topex.ucsd.edu/gmtsar/tar/S1A_Stack_CPGF_T173.tar.gz
./icloud_download.sh "https://www.icloud.com/iclouddrive/0b2SKaQ9FokMxD75yu6FF5lyg#S1A_Stack_CPGF_T173.tar.gz"
# unpack the dataset
rm -rf raw_orig topo && tar xvzf S1A_Stack_CPGF_T173.tar.gz -C . && rm -f raw_orig/._*
# cleanup
if [ "$1" != "" ]
then
    echo "Removing downloaded datasets to free disk space following 'force' command line argument..."
    rm *.tar.gz
fi

echo "Checking disk space"
df -h

echo "Running Python test script..."
rm -fr *.jpg
# prevent error
# Exception: "OSError(24, 'Too many open files')"
# default limit on MacOS is small (256)
ulimit -n 10000
python3 ./S1A_Stack_CPGF_T173.py && echo SUCCESS
