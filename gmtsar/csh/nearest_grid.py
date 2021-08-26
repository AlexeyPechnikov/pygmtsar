#!/usr/bin/env python3
# python3 -m pip install install xarray numpy scipy --upgrade

from scipy.spatial import cKDTree
import xarray as xr
import numpy as np
import sys

if not len(sys.argv) in [3,4]:
    print (f"Usage: {sys.argv[0]} input.grd output.grd [search_radius_pixels]")
    exit(0)

in_grid = sys.argv[1]
out_grid = sys.argv[2]
if len(sys.argv) == 4:
    search_radius = float(sys.argv[3])
else:
    search_radius = 300
#print ("in_grid", in_grid, "out_grid", out_grid, "search_radius", search_radius)

in0 = xr.open_dataarray(in_grid)
ys, xs = np.meshgrid(range(in0.shape[1]), range(in0.shape[0]))
ys = ys.reshape(-1)
xs = xs.reshape(-1)
zs = in0.values.reshape(-1)
mask = np.where(~np.isnan(zs))
# on regular source grid some options should be redefined
tree = cKDTree(np.column_stack((ys[mask],xs[mask])), compact_nodes=False, balanced_tree=False)
# use distance_limit
d, inds = tree.query(np.column_stack((ys,xs)), k = 1, distance_upper_bound=search_radius, workers=8)
# replace not available indexes by zero (see distance_upper_bound)
fakeinds = np.where(~np.isinf(d), inds, 0)
# produce the same output array
out1 = xr.zeros_like(in0)
out1.values = np.where(~np.isinf(d), zs[mask][fakeinds], np.nan).reshape(in0.shape)
compression = dict(zlib=True, complevel=3)
# magic: add GMT attribute to prevent coordinates shift for 1/2 pixel
out1.attrs['node_offset'] = 1
out1.to_netcdf(out_grid, encoding={out1.name: compression})
