# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2023, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_topo import Stack_topo
from .tqdm_dask import tqdm_dask

class Stack_phasediff(Stack_topo):

    @staticmethod
    def interferogram(phase):
        import numpy as np
        import warnings
        # suppress Dask warning "RuntimeWarning: invalid value encountered in divide"
        warnings.filterwarnings('ignore', module='dask')
        warnings.filterwarnings('ignore', module='dask.core')
        return np.arctan2(phase.imag, phase.real)

#     @staticmethod
#     def correlation(I1, I2, amp):
#         import xarray as xr
#         import numpy as np
#         # constant from GMTSAR code
#         thresh = 5.e-21
#         i = I1 * I2
#         corr = xr.where(i > 0, amp / np.sqrt(i), 0)
#         corr = xr.where(corr < 0, 0, corr)
#         corr = xr.where(corr > 1, 1, corr)
#         # mask too low amplitude areas as invalid
#         # amp1 and amp2 chunks are high for SLC, amp has normal chunks for NetCDF
#         return xr.where(i >= thresh, corr, np.nan).chunk(a.chunksizes).rename('phase')

    def correlation(self, phase, data):
        """
        Example:
        data_200m = stack.multilooking(np.abs(sbas.open_data()), wavelength=200, coarsen=(4,16))
        intf2_200m = stack.multilooking(intf2, wavelength=200, coarsen=(4,16))
        stack.correlation(intf2_200m, data_200m)

        Note:
        Multiple interferograms require the same data grids, allowing us to speed up the calculation
        by saving filtered data to a disk file.
        """
        import pandas as pd
        import dask
        import xarray as xr
        import numpy as np
        import warnings
        # suppress Dask warning "RuntimeWarning: invalid value encountered in divide"
        warnings.filterwarnings('ignore', module='dask')
        warnings.filterwarnings('ignore', module='dask.core')
        # constant from GMTSAR code
        # apply square root because we compare multiplication of amplitudes instead of intensities
        thresh = np.sqrt(5.e-21)

        # convert pairs (list, array, dataframe) to 2D numpy array
        pairs, dates = self.get_pairs(phase, dates=True)
        pairs = pairs[['ref', 'rep']].astype(str).values

        # check correctness for user-defined data argument
        assert not np.issubdtype(data.dtype, np.complexfloating), 'Use np.abs(data) to convert complex values to amplitudes before multilooking'

        stack = []
        for stack_idx, pair in enumerate(pairs):
            date1, date2 = pair
            # calculate correlation
            amps = data.sel(date=date1) * data.sel(date=date2)
            # mask too low amplitude areas as invalid
            corr = xr.where(amps >= thresh, np.abs(phase.sel(pair=' '.join(pair))) / amps, np.nan)
            corr = xr.where(corr < 0, 0, corr)
            corr = xr.where(corr > 1, 1, corr)
            del amps
            # add to stack
            stack.append(corr)
            del corr
        # cleanup
        del data

        return xr.concat(stack, dim='pair').rename('correlation')

    def phasediff(self, pairs, data='auto', topo='auto', method='cubic'):
        import pandas as pd
        import dask
        import xarray as xr
        import numpy as np

        # convert pairs (list, array, dataframe) to 2D numpy array
        pairs, dates = self.get_pairs(pairs, dates=True)
        pairs = pairs[['ref', 'rep']].astype(str).values

        if topo == 'auto':
            topo = self.get_topo()

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

        def block_phasediff(stack_idx, ylim, xlim):
            # disable "distributed.utils_perf - WARNING - full garbage collections ..."
            from dask.distributed import utils_perf
            utils_perf.disable_gc_diagnosis()

            # unpack input stacks
            prm1,  prm2  = stack_prm[stack_idx]
            data1, data2 = stack_data[stack_idx]

            # convert indices 0.5, 1.5,... to 0,1,... for easy calculations
            block_data1 = data1.isel(y=slice(ylim[0], ylim[1]), x=slice(xlim[0], xlim[1])).compute(n_workers=1)
            block_data2 = data2.isel(y=slice(ylim[0], ylim[1]), x=slice(xlim[0], xlim[1])).compute(n_workers=1)
            del data1, data2

            if abs(block_data1).sum() == 0:
                intf = np.nan * xr.zeros_like(block_data1)
                del block_data1, block_data2
                return intf

            ys = block_data1.y.astype(int)
            xs = block_data1.x.astype(int)

            block_data1 = block_data1.assign_coords(y=ys, x=xs)
            block_data2 = block_data2.assign_coords(y=ys, x=xs)

            # use outer variables topo, data1, data2, prm1, prm2
            # build topo block
            dy, dx = topo.y.diff('y').item(0), topo.x.diff('x').item(0)
            if dy == 1 and dx == 1:
                # topography is already in the original resolution
                block_topo = topo.isel(y=slice(ylim[0], ylim[1]), x=slice(xlim[0], xlim[1]))\
                            .compute(n_workers=1)\
                            .fillna(0)\
                            .assign_coords(y=ys, x=xs)
            else:
                # topography resolution is different, interpolation with extrapolation required
                # convert indices 0.5, 1.5,... to 0,1,... for easy calculations
                # fill NaNs by zero because typically DEM is missed outside of land areas
                block_topo = topo.sel(y=slice(ys[0]-3*dy, ys[-1]+3*dy), x=slice(xs[0]-3*dx, xs[-1]+3*dx))\
                            .compute(n_workers=1)\
                            .fillna(0)\
                            .interp({'y': block_data1.y, 'x': block_data1.x}, method=method, kwargs={'fill_value': 'extrapolate'})\
                            .assign_coords(y=ys, x=xs)
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

            from scipy import constants
            # setup the default parameters
            # constant from GMTSAR code for consistency
            #SOL = 299792456.0
            drange = constants.speed_of_light / (2 * prm2.get('rng_samp_rate'))
            #drange = SOL / (2 * prm2.get('rng_samp_rate'))
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
            # GMTSAR uses float32 complex operations with precision loss
            #phase_shift = np.exp(1j * (cnst * drho).astype(np.float32))
            phase_shift = np.exp(1j * (cnst * drho))
            del block_topo, near_range, drho, height, B, alpha, Bx, Bv, Bh, time

            # calculate phase difference
            intf = block_data1 * phase_shift * np.conj(block_data2)
            del block_data1, block_data2, phase_shift
            return intf.astype(np.complex64)

        if data == 'auto':
            # open datafiles required for all the pairs
            data = self.open_data(dates)
        # define blocks
        chunks = data.chunks
        ychunks,xchunks = chunks[1], chunks[2]
        ychunks = np.concatenate([[0], np.cumsum(ychunks)])
        xchunks = np.concatenate([[0], np.cumsum(xchunks)])
        ylims = [(y1, y2) for y1, y2 in zip(ychunks, ychunks[1:])]
        xlims = [(x1, x2) for x1, x2 in zip(xchunks, xchunks[1:])]
        #print ('ylims', ylims)
        #print ('xlims', xlims)

        stack_prm  = []
        stack_data = []
        stack = []
        for stack_idx, pair in enumerate(pairs):
            date1, date2 = pair

            data1 = data.sel(date=date1)
            data2 = data.sel(date=date2)

            # prepare for delayed stack processing
            stack_data.append((data1, data2))

            # prepare for delayed stack processing
            prm1 = self.PRM(date1)
            prm2 = self.PRM(date2)
            # it does not work because attributes are the same for all the grids
            #prm1 = PRM.from_str(data1.prm)
            #prm2 = PRM.from_str(data2.prm)
            # directory and filename required for SAT_... tools to locate LED file
            #prm1.filename = os.path.join(self.basedir, prm1.get('led_file'))
            #prm2.filename = os.path.join(self.basedir, prm2.get('led_file'))
            #print ('prm1.filename', prm1.filename)

            # update and add required parameters
            prm2.set(prm1.SAT_baseline(prm2, tail=9)).fix_aligned()
            # TEST - remove PRMs cross-linking
            #prm2.set(prm2.SAT_baseline(prm2, tail=9)).fix_aligned()
            prm1.set(prm1.SAT_baseline(prm1).sel('SC_height','SC_height_start','SC_height_end')).fix_aligned()
            stack_prm.append((prm1, prm2))

            # check the grids
            #assert prm1.get('num_rng_bins') == prm2.get('num_rng_bins'), 'The dimensions of range do not match'
            #assert prm1.get('num_patches') * prm1.get('num_valid_az') == prm2.get('num_patches') * prm2.get('num_valid_az'), \
            #    'The dimensions of azimuth do not match'

            if topo is None:
                # calculation is straightforward and does not require delayed wrappers
                intf = (data1 * np.conj(data2))
            else:
                # split to equal chunks and rest
                #ys_blocks = np.array_split(data[0].y, np.arange(0,data.y.size, self.chunksize)[1:])
                #xs_blocks = np.array_split(data[0].x, np.arange(0,data.x.size, self.chunksize)[1:])
                #print ('ys_blocks.size', len(ys_blocks), 'xs_blocks.size', len(xs_blocks))

                blocks_total = []
                for ylim in ylims:
                    blocks = []
                    for xlim in xlims:
                        block = dask.array.from_delayed(dask.delayed(block_phasediff)(stack_idx, ylim, xlim),
                                                        shape=((ylim[1]-ylim[0]), (xlim[1]-xlim[0])), dtype=np.complex64)
                        blocks.append(block)
                        del block
                    blocks_total.append(blocks)
                    del blocks
                intf = xr.DataArray(dask.array.block(blocks_total), coords={'y': data.y, 'x': data.x})
                del blocks_total

            # add to stack
            stack.append(intf)
            # cleanup
            del intf, data1, data2, prm1, prm2

        # cleanup
        del data

        coord_pair = [' '.join(pair) for pair in pairs]
        coord_ref = xr.DataArray(pd.to_datetime(pairs[:,0]), coords={'pair': coord_pair})
        coord_rep = xr.DataArray(pd.to_datetime(pairs[:,1]), coords={'pair': coord_pair})

        return xr.concat(stack, dim='pair').assign_coords(ref=coord_ref, rep=coord_rep, pair=coord_pair).rename('phasediff')

    @staticmethod
    def goldstein(phase, corr, psize=32):
        import xarray as xr
        import numpy as np
        import dask
        import warnings
        # suppress Dask warning "RuntimeWarning: invalid value encountered in divide"
        warnings.filterwarnings('ignore', module='dask')
        warnings.filterwarnings('ignore', module='dask.core')

        def apply_pspec(data, alpha):
            # NaN is allowed value
            assert not(alpha < 0), f'Invalid parameter value {alpha} < 0'
            wgt = np.power(np.abs(data)**2, alpha / 2)
            data = wgt * data
            return data

        def make_wgt(nxp, nyp):
            # Create arrays of horizontal and vertical weights
            wx = 1.0 - np.abs(np.arange(nxp // 2) - (nxp / 2.0 - 1.0)) / (nxp / 2.0 - 1.0)
            wy = 1.0 - np.abs(np.arange(nyp // 2) - (nyp / 2.0 - 1.0)) / (nyp / 2.0 - 1.0)
            # Compute the outer product of wx and wy to create the top-left quadrant of the weight matrix
            quadrant = np.outer(wy, wx)
            # Create a full weight matrix by mirroring the quadrant along both axes
            wgt = np.block([[quadrant, np.flip(quadrant, axis=1)],
                            [np.flip(quadrant, axis=0), np.flip(np.flip(quadrant, axis=0), axis=1)]])
            return wgt

        def patch_goldstein_filter(data, corr, wgt, psize):
            """
            Apply the Goldstein adaptive filter to the given data.

            Args:
                data: 2D numpy array of complex values representing the data to be filtered.
                corr: 2D numpy array of correlation values. Must have the same shape as `data`.

            Returns:
                2D numpy array of filtered data.
            """
            # Calculate alpha
            alpha = 1 - (wgt * corr).sum() / wgt.sum()
            data = np.fft.fft2(data, s=(psize,psize))
            data = apply_pspec(data, alpha)
            data = np.fft.ifft2(data, s=(psize,psize))
            return wgt * data

        def apply_goldstein_filter(data, corr, psize):
            # Create an empty array for the output
            out = np.zeros(data.shape, dtype=np.complex64)
            # ignore processing for empty chunks 
            if np.all(np.isnan(data)):
                return out
            # Create the weight matrix
            wgt_matrix = make_wgt(psize, psize)
            # Iterate over windows of the data
            for i in range(0, data.shape[0] - psize, psize // 2):
                for j in range(0, data.shape[1] - psize, psize // 2):
                    # Create proocessing windows
                    data_window = data[i:i+psize, j:j+psize]
                    corr_window = corr[i:i+psize, j:j+psize]
                    wgt_window = wgt_matrix[:data_window.shape[0],:data_window.shape[1]]
                    # Apply the filter to the window
                    filtered_window = patch_goldstein_filter(data_window, corr_window, wgt_window, psize)
                    # Add the result to the output array
                    slice_i = slice(i, min(i + psize, out.shape[0]))
                    slice_j = slice(j, min(j + psize, out.shape[1]))
                    out[slice_i, slice_j] += filtered_window[:slice_i.stop - slice_i.start, :slice_j.stop - slice_j.start]
            return out

        assert phase.shape == corr.shape, 'ERROR: phase and correlation variables have different shape'
        if len(phase.dims) == 2:
            stackvar = None
        else:
            stackvar = phase.dims[0]
    
        stack =[]
        for ind in range(len(phase) if stackvar is not None else 1):
            # Apply function with overlap; psize//2 overlap is not enough (some empty lines produced)
            # use complex data and real correlation
            block = dask.array.map_overlap(lambda phase, corr: apply_goldstein_filter(phase, corr, psize),
                                           (phase[ind] if stackvar is not None else phase).data,
                                           (corr[ind]  if stackvar is not None else corr).data,
                                           depth=psize // 2 + 2,
                                           dtype=np.complex64, 
                                           meta=np.array(()))
            # Calculate the phase
            stack.append(block)
            del block

        if stackvar is not None:
            ds = xr.DataArray(dask.array.stack(stack), coords=phase.coords)
        else:
            ds = xr.DataArray(stack[0], coords=phase.coords)
        del stack
        # replace zeros produces in NODATA areas
        return ds.where(ds).rename('phase')
