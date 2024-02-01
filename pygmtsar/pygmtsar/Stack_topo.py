# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2023, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_trans_inv import Stack_trans_inv

class Stack_topo(Stack_trans_inv):

    def get_topo(self):
        """
        Get the radar topography grid.

        Returns
        -------
        xarray.DataArray
            The 'topo' variable data is a xarray.DataArray.

        Examples
        --------
        Get DEM:
        topo = stack.get_topo()

        Notes
        -----
        This method returns 'topo' variable from inverse radar transform grid.
        """
        return self.get_trans_inv()['ele'].rename('topo')

    def plot_topo(self, topo='auto', caption='Topography on WGS84 ellipsoid, [m]',
                  quantile=None, vmin=None, vmax=None, cmap='gray', aspect=None, **kwargs):
        import numpy as np
        import matplotlib.pyplot as plt

        if isinstance(topo, str) and topo == 'auto':
            topo = self.get_topo()
            
        if quantile is not None:
            assert vmin is None and vmax is None, "ERROR: arguments 'quantile' and 'vmin', 'vmax' cannot be used together"

        if quantile is not None:
            vmin, vmax = np.nanquantile(topo, quantile)

        plt.figure()
        topo.plot.imshow(cmap=cmap, vmin=vmin, vmax=vmax)
        self.plot_AOI(**kwargs)
        self.plot_POI(**kwargs)
        if aspect is not None:
            plt.gca().set_aspect(aspect)
        plt.xlabel('Range')
        plt.ylabel('Azimuth')
        plt.title(caption)

    def topo_phase(self, pairs, topo='auto', debug=False):
        """
        decimator = sbas.decimator(resolution=15, grid=(1,1))
        topophase = sbas.topo_phase(pairs, topo=decimator(sbas.get_topo()))
        """
        import pandas as pd
        import dask
        import dask.dataframe
        import xarray as xr
        import numpy as np
        #from tqdm.auto import tqdm
        import joblib
        import warnings
        # suppress Dask warning "RuntimeWarning: invalid value encountered in divide"
        warnings.filterwarnings('ignore')
        warnings.filterwarnings('ignore', module='dask')
        warnings.filterwarnings('ignore', module='dask.core')

        if debug:
            print ('DEBUG: topo_phase')

        # convert pairs (list, array, dataframe) to 2D numpy array
        pairs, dates = self.get_pairs(pairs, dates=True)
        pairs = pairs[['ref', 'rep']].astype(str).values

        if isinstance(topo, str) and topo == 'auto':
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

        def block_phase(prm1, prm2, ylim, xlim):
            # use outer variables date, stack_prm
            # disable "distributed.utils_perf - WARNING - full garbage collections ..."
            from dask.distributed import utils_perf
            utils_perf.disable_gc_diagnosis()
            import warnings
            # suppress Dask warning "RuntimeWarning: invalid value encountered in divide"
            warnings.filterwarnings('ignore')
            warnings.filterwarnings('ignore', module='dask')
            warnings.filterwarnings('ignore', module='dask.core')

            # for lazy Dask ddataframes
            #prm1 = PRM(prm1)
            #prm2 = PRM(prm2)
            #prm1,  prm2  = stack_prm[stack_idx]
            #data1, data2 = stack_data[stack_idx]

            # use outer variables topo, data1, data2, prm1, prm2
            # build topo block
            block_topo = topo.isel(y=slice(ylim[0], ylim[1]), x=slice(xlim[0], xlim[1]))\
                        .compute(n_workers=1)\
                        .fillna(0)

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

            return phase_shift.astype(np.complex64)

    #     # prepare lazy PRM
    #     # this is the "true way" but processing is ~40% slower due to additional Dask tasks
    #     def block_prms(date1, date2):
    #         prm1 = self.PRM(date1)
    #         prm2 = self.PRM(date2)
    #         prm2.set(prm1.SAT_baseline(prm2, tail=9)).fix_aligned()
    #         prm1.set(prm1.SAT_baseline(prm1).sel('SC_height','SC_height_start','SC_height_end')).fix_aligned()
    #         return (prm1.df, prm2.df)
    #     # Define metadata explicitly to match PRM dataframe
    #     prm_meta = pd.DataFrame(columns=['name', 'value']).astype({'name': 'str', 'value': 'object'}).set_index('name')

        # immediately prepare PRM
        # here is some delay on the function call but the actual processing is faster
        def prepare_prms(pair):
            date1, date2 = pair
            prm1 = self.PRM(date1)
            prm2 = self.PRM(date2)
            prm2.set(prm1.SAT_baseline(prm2, tail=9)).fix_aligned()
            prm1.set(prm1.SAT_baseline(prm1).sel('SC_height','SC_height_start','SC_height_end')).fix_aligned()
            return {(date1, date2): (prm1, prm2)}

        #with self.tqdm_joblib(tqdm(desc=f'Pre-Processing PRM', total=len(pairs))) as progress_bar:
        prms = joblib.Parallel(n_jobs=-1)(joblib.delayed(prepare_prms)(pair) for pair in pairs)
        # convert the list of dicts to a single dict
        prms = {k: v for d in prms for k, v in d.items()}

        # define blocks
        chunks = topo.chunks
        ychunks, xchunks = chunks[0], chunks[1]
        ychunks = np.concatenate([[0], np.cumsum(ychunks)])
        xchunks = np.concatenate([[0], np.cumsum(xchunks)])
        ylims = [(y1, y2) for y1, y2 in zip(ychunks, ychunks[1:])]
        xlims = [(x1, x2) for x1, x2 in zip(xchunks, xchunks[1:])]
        #print ('ylims', ylims)
        #print ('xlims', xlims)

        stack = []
        for stack_idx, pair in enumerate(pairs):
            date1, date2 = pair

            # Create a Dask DataFrame with provided metadata for each Dask block
            #prms = dask.delayed(block_prms)(date1, date2)
            #prm1 = dask.dataframe.from_delayed(dask.delayed(prms[0]), meta=prm_meta)
            #prm2 = dask.dataframe.from_delayed(dask.delayed(prms[1]), meta=prm_meta)
            prm1, prm2 = prms[(date1, date2)]

            blocks_total = []
            for ylim in ylims:
                blocks = []
                for xlim in xlims:
                    delayed = dask.delayed(block_phase)(prm1, prm2, ylim, xlim)
                    block = dask.array.from_delayed(delayed,
                                                    shape=((ylim[1]-ylim[0]), (xlim[1]-xlim[0])),
                                                    dtype=np.complex64)
                    blocks.append(block)
                    del block, delayed
                blocks_total.append(blocks)
                del blocks
            intf = xr.DataArray(dask.array.block(blocks_total), coords={'y': topo.y, 'x': topo.x})
            del blocks_total

            # add to stack
            stack.append(intf)
            # cleanup
            del intf, prm1, prm2
        del prms

        coord_pair = [' '.join(pair) for pair in pairs]
        coord_ref = xr.DataArray(pd.to_datetime(pairs[:,0]), coords={'pair': coord_pair})
        coord_rep = xr.DataArray(pd.to_datetime(pairs[:,1]), coords={'pair': coord_pair})

        return xr.concat(stack, dim='pair').assign_coords(ref=coord_ref, rep=coord_rep, pair=coord_pair).where(np.isfinite(topo)).rename('phase')

    def topo_slope(self, topo='auto', edge_order=1):
        import xarray as xr
        import numpy as np

        if isinstance(topo, str) and topo == 'auto':
            topo = self.get_topo()

        dy, dx = self.get_spacing(topo)
        # The cell sizes (dx, dy) are passed to np.gradient to account for non-uniform spacing
        gradient_y, gradient_x = np.gradient(topo, dy, dx, edge_order=edge_order)
        slope_radians = np.arctan(np.sqrt(gradient_x**2 + gradient_y**2))
        slope_degrees = np.degrees(slope_radians)
        slope = xr.DataArray(slope_degrees, coords=topo.coords, dims=topo.dims, name="slope_degrees")
        slope.attrs["units"] = "degrees"
        slope.attrs["description"] = "Slope calculated from DEM, accounting for non-uniform grid spacing"

        return slope

    def topo_variation(self):
        """
        Topography related range cell variation to detect layovers and shadows.
        """
        drange = self.get_trans_inv().ll.diff('x')
        return drange/drange.mean()
