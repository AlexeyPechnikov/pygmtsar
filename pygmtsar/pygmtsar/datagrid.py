#!/usr/bin/env python3
# Alexey Pechnikov, Oct, 2022, https://github.com/mobigroup/gmtsar

class datagrid:

    # NetCDF options
    chunksize = 512
    compression = dict(zlib=True, complevel=3, chunksizes=(chunksize, chunksize))
    engine = 'h5netcdf'
