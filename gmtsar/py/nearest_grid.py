#!/usr/bin/env python3
# Alexey Pechnikov, Aug, 2021, https://github.com/mobigroup/gmtsar
# python3 -m pip install install xarray numpy scipy --upgrade
import pytest

def nearest_grid(in_grid_filename,  search_radius=300):
    from scipy.spatial import cKDTree
    import xarray as xr
    import numpy as np

    in_grid = xr.open_dataarray(in_grid_filename)
    ys, xs = np.meshgrid(range(in_grid.shape[1]), range(in_grid.shape[0]))
    ys = ys.reshape(-1)
    xs = xs.reshape(-1)
    zs = in_grid.values.reshape(-1)
    mask = np.where(~np.isnan(zs))
    # on regular source grid some options should be redefined for better performance
    tree = cKDTree(np.column_stack((ys[mask],xs[mask])), compact_nodes=False, balanced_tree=False)
    # use distance_limit
    d, inds = tree.query(np.column_stack((ys,xs)), k = 1, distance_upper_bound=search_radius, workers=8)
    # replace not available indexes by zero (see distance_upper_bound)
    fakeinds = np.where(~np.isinf(d), inds, 0)
    # produce the same output array as dataset to be able to add global attributes
    out_grid = xr.zeros_like(in_grid).to_dataset()
    out_grid['z'].values = np.where(~np.isinf(d), zs[mask][fakeinds], np.nan).reshape(in_grid.shape)
    # magic: add GMT attribute to prevent coordinates shift for 1/2 pixel
    out_grid.attrs['node_offset'] = 1
    return out_grid

def main():
    import sys

    if not len(sys.argv) in [3,4]:
        print (f"Usage: {sys.argv[0]} input.grd output.grd [search_radius_pixels]")
        exit(0)

    in_grid_filename = sys.argv[1]
    out_grid_filename = sys.argv[2]
    if len(sys.argv) == 4:
        search_radius = float(sys.argv[3])
        out_grid = nearest_grid(in_grid_filename, search_radius)
    else:
        out_grid = nearest_grid(in_grid_filename)
    #print ("in_grid", in_grid, "out_grid", out_grid, "search_radius", search_radius)
    # save to NetCDF file
    compression = dict(zlib=True, complevel=3, chunksizes=[128,128])
    out_grid.to_netcdf(out_grid_filename, encoding={'z': compression})

def test_main():
    import subprocess
    import os

    in_grid_filename = '../testdata/nearest_grid.grd'
    out_grid_filename = '../testdata/_pytest.nearest_grid.grd'

    output = """\
: Title: z
: Command: 
: Remark: 
: Pixel node registration used [Geographic grid]
: Grid file format: nf = GMT netCDF format (32-bit float), CF-1.7
: x_min: 240.866666667 x_max: 244.033333333 x_inc: 0.00166666666667 (6 sec) name: longitude n_columns: 1900
: y_min: 32.3194444444 y_max: 36.1111111111 y_inc: 0.00138888888889 (5 sec) name: latitude n_rows: 2730
: v_min: -276.653930664 v_max: 26.5470123291 name: z
: scale_factor: 1 add_offset: 0
: 219936 nodes (4.2%) set to NaN
: mean: -66.259027874 stdev: 30.5838879722 rms: 72.976933393
: format: netCDF-4 chunk_size: 128,128 shuffle: on deflation_level: 3
"""
    assert os.path.exists(in_grid_filename)
    out_grid = nearest_grid(in_grid_filename)
    # save to NetCDF file
    compression = dict(zlib=True, complevel=3, chunksizes=[128,128])
    out_grid.to_netcdf(out_grid_filename, encoding={'z': compression})
    # check the NetCDF file
    result = subprocess.run(['gmt', 'grdinfo', '-L2', out_grid_filename], stdout=subprocess.PIPE)
    result = result.stdout.decode('utf-8').replace(out_grid_filename,'')
    assert result == output
    os.remove(out_grid_filename)

if __name__ == "__main__":
    # execute only if run as a script
    main()
