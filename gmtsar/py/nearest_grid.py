#!/usr/bin/env python3
# python3 -m pip install install xarray numpy scipy --upgrade
import pytest

def nearest_grid(in_grid, out_grid, search_radius=300):
    from scipy.spatial import cKDTree
    import xarray as xr
    import numpy as np

    in0 = xr.open_dataarray(in_grid)
    ys, xs = np.meshgrid(range(in0.shape[1]), range(in0.shape[0]))
    ys = ys.reshape(-1)
    xs = xs.reshape(-1)
    zs = in0.values.reshape(-1)
    mask = np.where(~np.isnan(zs))
    # on regular source grid some options should be redefined for better performance
    tree = cKDTree(np.column_stack((ys[mask],xs[mask])), compact_nodes=False, balanced_tree=False)
    # use distance_limit
    d, inds = tree.query(np.column_stack((ys,xs)), k = 1, distance_upper_bound=search_radius, workers=8)
    # replace not available indexes by zero (see distance_upper_bound)
    fakeinds = np.where(~np.isinf(d), inds, 0)
    # produce the same output array as dataset to be able to add global attributes
    out1 = xr.zeros_like(in0).to_dataset()
    out1['z'].values = np.where(~np.isinf(d), zs[mask][fakeinds], np.nan).reshape(in0.shape)
    compression = dict(zlib=True, complevel=3, chunksizes=[128,128])
    # magic: add GMT attribute to prevent coordinates shift for 1/2 pixel
    out1.attrs['node_offset'] = np.int32(1)
    out1.to_netcdf(out_grid, encoding={'z': compression})

def main():
    import sys

    if not len(sys.argv) in [3,4]:
        print (f"Usage: {sys.argv[0]} input.grd output.grd [search_radius_pixels]")
        exit(0)

    in_grid = sys.argv[1]
    out_grid = sys.argv[2]
    if len(sys.argv) == 4:
        search_radius = float(sys.argv[3])
        nearest_grid(in_grid, out_grid, search_radius)
    else:
        nearest_grid(in_grid, out_grid)
    #print ("in_grid", in_grid, "out_grid", out_grid, "search_radius", search_radius)

def test_main():
    import subprocess
    import os

    in_grid = '../testdata/nearest_grid.grd'
    out_grid = '../testdata/_pytest.nearest_grid.grd'

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
    assert os.path.exists(in_grid)
    nearest_grid(in_grid, out_grid)
    result = subprocess.run(['gmt', 'grdinfo', '-L2', out_grid], stdout=subprocess.PIPE)
    result = result.stdout.decode('utf-8').replace(out_grid,'')
    assert result == output
    os.remove(out_grid)

if __name__ == "__main__":
    # execute only if run as a script
    main()
