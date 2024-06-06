#!/bin/sh

rm -f SAT_baseline
gcc \
    -I/usr/local/GMTSAR/gmtsar -L/usr/local/GMTSAR/gmtsar -lgmtsar \
    -I/opt/homebrew/Cellar/gmt/6.5.0_1/include/gmt -L/opt/homebrew/Cellar/gmt/6.5.0_1/include/gmt \
    -o SAT_baseline SAT_baseline.c
./SAT_baseline S1_20201222_ALL_F2.PRM S1_20210103_ALL_F2.PRM
