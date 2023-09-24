# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2021, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_topo import Stack_topo
from .tqdm_dask import tqdm_dask

class Stack_phasediff(Stack_topo):

    def phasediff(self, pair, topo='auto', scale=2.5e-07, method='cubic', chunksize=None):
        import dask
        import xarray as xr
        import numpy as np

        # TODO
        date1, date2 = pair

        if topo == 'auto':
            topo = self.get_topo()

        if chunksize is None:
            # SLC chunksize is unlimited for x dimension and chunked for y
            chunksize = self.chunksize // 8

        # calculate the combined earth curvature and topography correction
        def calc_drho(rho, topo, earth_radius, height, b, alpha, Bx):
            sina = np.sin(alpha)
            cosa = np.cos(alpha)
            c = earth_radius + height
            # compute the look angle using equation (C26) in Appendix C
            # GMTSAR uses long double here
            #ret = earth_radius + topo.astype(np.longdouble)
            ret = earth_radius + topo
            cost = ((rho**2 + c**2 - ret**2) / (2. * rho * c))
            #if (cost >= 1.)
            #    die("calc_drho", "cost >= 0");
            sint = np.sqrt(1. - cost**2)
            # Compute the offset effect from non-parallel orbit
            term1 = rho**2 + b**2 - 2 * rho * b * (sint * cosa - cost * sina) - Bx**2
            drho = -rho + np.sqrt(term1)
            del term1, sint, cost, ret, c, cosa, sina
            return drho

        def block_phasediff(ys, xs):
            # use outer variables topo, data1, data2, prm1, prm2, scale
            # build topo block
            dy, dx = topo.y.diff('y').item(0), topo.x.diff('x').item(0)
            # convert indices 0.5, 1.5,... to 0,1,... for easy calculations
            # fill NaNs by zero because typically DEM is missed outside of land areas
            block_topo = topo.sel(y=slice(ys[0]-2*dy, ys[-1]+2*dy), x=slice(xs[0]-2*dx, xs[-1]+2*dx))\
                        .compute(n_workers=1)\
                        .fillna(0)\
                        .interp({'y': ys, 'x': xs}, method=method)\
                        .assign_coords(y=ys.astype(int), x=xs.x.astype(int))
            # set dimensions
            xdim = prm1.get('num_rng_bins')
            ydim = prm1.get('num_patches') * prm1.get('num_valid_az')

            # set heights
            htc = prm1.get('SC_height')
            ht0 = prm1.get('SC_height_start')
            htf = prm1.get('SC_height_end')

            # compute the time span and the time spacing
            tspan = 86400 * abs(prm2.get('SC_clock_stop') - prm2.get('SC_clock_start'))
            assert (tspan >= 0.01) and (prm2.get('PRF') >= 0.01), 'Check sc_clock_start, sc_clock_end, or PRF'

            #from scipy import constants
            # setup the default parameters
            # constant from GMTSAR code for consistency
            SOL = 299792456.0
            #drange = constants.speed_of_light / (2 * prm2.get('rng_samp_rate'))
            drange = SOL / (2 * prm2.get('rng_samp_rate'))
            alpha = prm2.get('alpha_start') * np.pi / 180
            cnst = -4 * np.pi / prm2.get('radar_wavelength')

            # calculate initial baselines
            Bh0 = prm2.get('baseline_start') * np.cos(prm2.get('alpha_start') * np.pi / 180)
            Bv0 = prm2.get('baseline_start') * np.sin(prm2.get('alpha_start') * np.pi / 180)
            Bhf = prm2.get('baseline_end')   * np.cos(prm2.get('alpha_end')   * np.pi / 180)
            Bvf = prm2.get('baseline_end')   * np.sin(prm2.get('alpha_end')   * np.pi / 180)
            Bx0 = prm2.get('B_offset_start')
            Bxf = prm2.get('B_offset_end')

            # first case is quadratic baseline model, second case is default linear model
            if prm2.get('baseline_center') != 0 or prm2.get('alpha_center') != 0 or prm2.get('B_offset_center') != 0:
                Bhc = prm2.get('baseline_center') * np.cos(prm2.get('alpha_center') * np.pi / 180)
                Bvc = prm2.get('baseline_center') * np.sin(prm2.get('alpha_center') * np.pi / 180)
                Bxc = prm2.get('B_offset_center')

                dBh = (-3 * Bh0 + 4 * Bhc - Bhf) / tspan
                dBv = (-3 * Bv0 + 4 * Bvc - Bvf) / tspan
                ddBh = (2 * Bh0 - 4 * Bhc + 2 * Bhf) / (tspan * tspan)
                ddBv = (2 * Bv0 - 4 * Bvc + 2 * Bvf) / (tspan * tspan)

                dBx = (-3 * Bx0 + 4 * Bxc - Bxf) / tspan
                ddBx = (2 * Bx0 - 4 * Bxc + 2 * Bxf) / (tspan * tspan)
            else:
                dBh = (Bhf - Bh0) / tspan
                dBv = (Bvf - Bv0) / tspan
                dBx = (Bxf - Bx0) / tspan
                ddBh = ddBv = ddBx = 0

            # calculate height increment
            dht = (-3 * ht0 + 4 * htc - htf) / tspan
            ddht = (2 * ht0 - 4 * htc + 2 * htf) / (tspan * tspan)

            # multiply xr.ones_like(topo) for correct broadcasting
            near_range = xr.ones_like(block_topo)*(prm1.get('near_range') + \
                block_topo.x * (1 + prm1.get('stretch_r')) * drange) + \
                xr.ones_like(block_topo)*(block_topo.y * prm1.get('a_stretch_r') * drange)

            # calculate the change in baseline and height along the frame if topoflag is on
            time = block_topo.y * tspan / (ydim - 1)        
            Bh = Bh0 + dBh * time + ddBh * time**2
            Bv = Bv0 + dBv * time + ddBv * time**2
            Bx = Bx0 + dBx * time + ddBx * time**2
            B = np.sqrt(Bh * Bh + Bv * Bv)
            alpha = np.arctan2(Bv, Bh)
            height = ht0 + dht * time + ddht * time**2

            # calculate the combined earth curvature and topography correction
            drho = calc_drho(near_range, block_topo, prm1.get('earth_radius'), height, B, alpha, Bx)

            # make topographic and model phase corrections
            # GMTSAR uses float32 complex operations with precision lost
            #phase_shift = np.exp(1j * (cnst * drho).astype(np.float32))
            phase_shift = np.exp(1j * (cnst * drho))
            del block_topo, near_range, drho, height, B, alpha, Bx, Bv, Bh, time

            # calculate phase difference
            # convert indices 0.5, 1.5,... to 0,1,... for easy calculations
            block_data1 = data1.sel(y=ys, x=xs)\
                .compute(n_workers=1)\
                .assign_coords(y=ys.astype(int), x=xs.x.astype(int))
            block_data2 = data2.sel(y=ys, x=xs)\
                .compute(n_workers=1)\
                .assign_coords(y=ys.astype(int), x=xs.x.astype(int))
        
            intf = block_data1 * phase_shift * np.conj(block_data2)
            del block_data1, block_data2, phase_shift
            return intf.where(abs(intf)>scale**2).astype(np.complex64)
    
        prm1 = self.PRM(date1)
        prm2 = self.PRM(date2)

        # update and add required parameters
        prm2.set(prm1.SAT_baseline(prm2, tail=9)).fix_aligned()
        prm1.set(prm1.SAT_baseline(prm1).sel('SC_height','SC_height_start','SC_height_end')).fix_aligned()

        # check the grids
        assert prm1.get('num_rng_bins') == prm2.get('num_rng_bins'), 'The dimensions of range do not match'
        assert prm1.get('num_patches') * prm1.get('num_valid_az') == prm2.get('num_patches') * prm2.get('num_valid_az'), \
            'The dimensions of azimuth do not match'

        # read datafiles and convert to complex values
        slc1 = prm1.read_SLC_int(scale=scale, chunksize=chunksize)
        data1 = (slc1.re + 1j * slc1.im)
        slc2 = prm2.read_SLC_int(scale=scale, chunksize=chunksize)
        data2 = (slc2.re + 1j * slc2.im)
        del slc1, slc2
    
        if topo is None:
            # calculation is straightforward and does not require delayed wrappers
            return (data1 * np.conj(data2)).where(abs(data1)>scale**2).rename('phasediff')
    
        ychunks, xchunks = data1.chunks
        ychunksize = ychunks[0]
        xchunksize = xchunks[0]
        #print (ychunksize, xchunksize)
    
        # split to equal chunks and rest
        ys_blocks = np.array_split(data1.y, np.arange(0,data1.y.size, ychunksize)[1:])
        xs_blocks = np.array_split(data1.x, np.arange(0,data1.x.size, xchunksize)[1:])
        print ('ys_blocks.size', len(ys_blocks), 'xs_blocks.size', len(xs_blocks))
    
        blocks_total = []
        for ys_block in ys_blocks:
            blocks = []
            for xs_block in xs_blocks:
                block = dask.array.from_delayed(dask.delayed(block_phasediff)(ys_block, xs_block),
                                                shape=(ys_block.size, xs_block.size), dtype=np.complex64)
                blocks.append(block)
                del block
            blocks_total.append(blocks)
            del blocks
        intf = xr.DataArray(dask.array.block(blocks_total), coords=data1.coords).rename('phasediff')
        del blocks_total
        return intf

    def stack_phasediff(self, pairs, topo_fromfile=None, chunksize=None, n_jobs=-1):
        """
        Build phase difference for all pairs.

        Parameters
        ----------
        pairs : list
            List of date pairs (baseline pairs).
        topo_fromfile : str, optional
            Define radar coordinate topography. Default is topography created by geocode().
        n_jobs : int, optional
            Number of parallel processing jobs. n_jobs=-1 means all the processor cores used.
        Returns
        -------
        None
        """
        from tqdm.auto import tqdm
        import joblib

        # convert pairs (list, array, dataframe) to 2D numpy array
        pairs = self.get_pairs(pairs)[['ref', 'rep']].astype(str).values

        with self.tqdm_joblib(tqdm(desc='Phase Difference', total=len(pairs))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(self.phasediff)(pair, topo_fromfile, chunksize) \
                for pair in pairs)
