#!/bin/sh
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
# Replacement for prep_data_linux.csh and prep_data.csh
# suppose required and only required orbit files already downloaded by "eof" tool
set -e

ls *.EOF  | sort -t '_' -k7 > orbits.list
ls *.tiff | sort | xargs -I {} -n 1 basename {} .tiff > scenes.list
paste -d ':' scenes.list orbits.list > data.in
