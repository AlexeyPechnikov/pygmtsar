#!/bin/sh

wget -c http://topex.ucsd.edu/gmtsar/tar/S1A_Stack_CPGF_T173.tar.gz
tar xvzf S1A_Stack_CPGF_T173.tar.gz -C .
python3 ./S1A_Stack_CPGF_T173.py && echo SUCCESS
