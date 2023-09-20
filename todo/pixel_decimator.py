# use weighted decimation
# amplify by weights and multi-look and divide by weights multi-look
#decimator = lambda da: da.coarsen({'y': 2, 'x': 2}, boundary='trim').mean()
def pixel_decimator(self, resolution_meters=60, grid=(1, 4), debug=False):
    """
    Return function for pixel decimation to the specified output resolution.

    Parameters
    ----------
    resolution_meters : int, optional
        DEM grid resolution in meters. The same grid is used for geocoded results output.
    grid : tuple, optional
        Grid size for pixel decimation in the format (vertical, horizontal).
    debug : bool, optional
        Boolean flag to print debug information.

    Returns
    -------
    callable
        Post-processing function for SBAS.ints() and SBAS.intf_parallel().

    Examples
    --------
    Decimate computed interferograms to default DEM resolution 60 meters:
    decimator = sbas.pixel_decimator()
    sbas.intf_parallel(pairs, func=decimator)
    """
    import numpy as np
    import dask

    # special cases: scale factor should be 4*N or 2*N to prevent rounding issues
    # grid can be defined as [] or () or xarray dataarray
    grid_as_coeffs = isinstance(grid,(list, tuple))
    if   (grid_as_coeffs and grid==(1, 1)) or (not grid_as_coeffs and grid.x.diff('x')[0].item()) == 1:
        xscale0 = 4
    elif (grid_as_coeffs and grid==(1, 2)) or (not grid_as_coeffs and grid.x.diff('x')[0].item()) == 2:
        xscale0 = 2
    else:
        xscale0 = 1
    if debug:
        print (f'DEBUG: scale to square grid: xscale0={xscale0}')

    dy, dx = self.pixel_size(grid)
    yscale, xscale = int(np.round(resolution_meters/dy)), int(np.round(resolution_meters/dx/xscale0))
    if debug:
        print (f'DEBUG: average per subswaths ground pixel size in meters: y={dy}, x={dx}')
    if yscale <= 1 and xscale <= 1:
        # decimation is not possible
        if debug:
            print (f'DEBUG: decimator = lambda da: da')
        return lambda da: da

    # weighted decimate function
    def decimator_weighted(da, weight=None):
        # workaround for Google Colab when we cannot save grids with x,y coordinate names
        # also supports geographic coordinates
        yname = [varname for varname in ['y', 'lat', 'a'] if varname in da.dims][0]
        xname = [varname for varname in ['x', 'lon', 'r'] if varname in da.dims][0]
        #print ('da', da)
        #if debug:
        #    print (f"Decimate y variable '{yname}' for scale 1/{yscale} and x variable '{xname}' for scale 1/{xscale}")
        # avoid creating the large chunks
        with dask.config.set(**{'array.slicing.split_large_chunks': True}):
            # weight dimensions should be (y,x), NAN values in raster ignored in weighting
            if weight is not None:
                assert da.shape == weight.shape, 'Array and weight shapes are different'
                # unify grid coordinate names
                # process in fake coordinates to prevent square grid issue
                if yname != 'a':
                    da = da.rename({yname: 'a', xname: 'r'})
                print ('da', da)
                wyname = [varname for varname in ['y', 'lat', 'a'] if varname in weight.dims][0]
                wxname = [varname for varname in ['x', 'lon', 'r'] if varname in weight.dims][0]
                # process in fake coordinates — to prevent square grid issue
                if wyname != 'a':
                    weight = weight.rename({wyname: 'a', wxname: 'r'})
                # rechunk if needed
                #weight = weight.chunk(da.chunks)
                print ('weight', weight)
                return (weight*da).coarsen({'a': yscale, 'r': xscale0*xscale}, boundary='trim').sum()/\
                       weight.where(~np.isnan(da)).coarsen({'a': yscale, 'r': xscale0*xscale}, boundary='trim').sum()
            else:
                #print ('yname', yname, 'xname', xname)
                #print ('da', da)
                # do not rename dimensions — to prevent square grid issue, it is opposite to the case above
                return da.coarsen({yname: yscale, xname: xscale0*xscale}, boundary='trim').mean()

    # return callback function
    if debug:
        print (f"""DEBUG: decimator_weighted = lambda da, weight: (weight*da).coarsen({{'y': {yscale}, 'x': {xscale0*xscale}}}, boundary='trim').sum()/
               weight.where(~np.isnan(da)).coarsen({{'y': {yscale}, 'x': {xscale0*xscale}}}, boundary='trim').sum()""")
        print (f"DEBUG: decimator = lambda da: da.coarsen({{'y': {yscale}, 'x': {xscale0*xscale}}}, boundary='trim').sum()")
    return lambda da, weight=None: decimator_weighted(da, weight)
    # TEST
    #return lambda da, weight=None: da

SBAS.pixel_decimator = pixel_decimator

#decimator = sbas.pixel_decimator(resolution_meters=60, debug=True)
#sbas.intf(1, pairs.head(1), wavelength=200, psize=None, func=decimator)
#sbas.intf(1, pairs.head(1), wavelength=100, psize=None, weight=1/adi, func=decimator, coarsen=None, debug=True)
