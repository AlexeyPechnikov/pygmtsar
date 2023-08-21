# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2021, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .SBAS_topo import SBAS_topo
from .tqdm_dask import tqdm_dask

class SBAS_intf(SBAS_topo):

    def intf(self, subswath, pair, **kwargs):
        """
        Generate an interferogram for a given pair of dates representing synthetic aperture radar (SAR) images.

        Parameters:
        ----------
        subswath : int
            The subswath number to use.

        pair : tuple
            A tuple of strings representing the dates for the pair of images in the form 'YYYYMMDD'.

        kwargs : dict, optional
            Additional keyword arguments for the PRM.intf() method. This can be used to customize the interferogram generation process.

        Raises:
        ------
        OSError
            If the required input files or directories cannot be accessed.

        Notes:
        ------
        This method generates an interferogram by first creating PRM objects for the reference (ref) and repeat (rep) images. It then calls the intf() method of the PRM object for the reference image, passing the PRM object for the repeat image and any additional keyword arguments. The output is an interferogram stored as a grid file.
        """
        import pandas as pd
        import numpy as np

        # convert to 2D single-element array
        if isinstance(pair, pd.DataFrame):
            assert len(pair) == 1, 'Only single-record DataFrame or 1 D or 2D pair of dates allowed'
            pair = self.get_pairs(pair)[['ref','rep']].astype(str).values
        else:
            pair = [pair] if np.asarray(pair).ndim == 1 else pair

        # extract dates from pair
        date1, date2 = pair[0]

        prm_ref = self.PRM(subswath, date1)
        prm_rep = self.PRM(subswath, date2)

        if 'topo_file' in kwargs:
            topo_file = kwargs['topo_file']
            del kwargs['topo_file']
        else:
            topo_file = self.get_filename('topo', subswath)
        #print ('SBAS intf kwargs', kwargs)
        prm_ref.intf(prm_rep,
                     basedir=self.basedir,
                     topo_fromfile = topo_file,
                     **kwargs)

    def intf_parallel(self, pairs, weight=None, n_jobs=-1, chunksize=None, **kwargs):
        """
        Build interferograms for all the subswaths in parallel.

        Parameters
        ----------
        pairs : list
            List of date pairs (baseline pairs).
        n_jobs : int, optional
            Number of parallel processing jobs. n_jobs=-1 means all the processor cores used.
        wavelength : float, optional
            Filtering wavelength in meters.
        psize : int, optional
            Patch size for modified Goldstein adaptive filter (power of two).
        func : function, optional
            Post-processing function usually used for decimation.

        Returns
        -------
        None

        Examples
        --------
        For default 60m DEM resolution and other default parameters use command below:
        pairs = [sbas.to_dataframe().index.unique()]
        decimator = lambda dataarray: dataarray.coarsen({'y': 4, 'x': 4}, boundary='trim').mean()
        sbas.intf_parallel(pairs, func=decimator)
        """
        import xarray as xr
        import pandas as pd
        import numpy as np
        import dask
        from tqdm.auto import tqdm
        import joblib
        import os

        # convert pairs (list, array, dataframe) to 2D numpy array
        pairs = self.get_pairs(pairs)[['ref', 'rep']].astype(str).values

        subswaths = self.get_subswaths()

        # for now (Python 3.10.10 on MacOS) joblib loads the code from disk instead of copying it
        kwargs['chunksize'] = chunksize

        # materialize lazy weights
        if weight is not None and isinstance(weight, xr.DataArray):
            if len(subswaths) == 1:
                weights = [weight]
            else:
                raise ValueError(f"Argument weight should be a list or a tuple of DataArray corresponding to subswaths")
        elif weight is not None and isinstance(weight, (list, tuple)):
            weights = weight
        else:
            # form list of None
            weights = [None for swath in subswaths]
        
    #     if weight is not None and isinstance(weight, (list, tuple)):
    #         for idx, subswath in enumerate(subswaths):
    #             weight_filename = self.get_filename('intfweight', subswath)
    #             if os.path.exists(weight_filename):
    #                 os.remove(weight_filename)
    #             # workaround to save NetCDF file correct
    #             handler = weight[idx].rename('weight').rename({'y':'a','x':'r'}).\
    #                 to_netcdf(weight_filename,
    #                           encoding={'weight': self.compression(weight[idx].shape, chunksize=chunksize)},
    #                           engine=self.engine,
    #                           compute=False)
    #             tqdm_dask(dask.persist(handler), desc=f'Materialize weight sw{subswath}')
    #         # set base name for all subswaths
    #         kwargs['weight'] = 'intfweight'

        # this way does not work properly for long interferogram series on MacOS
        # see https://github.com/mobigroup/gmtsar/commit/3eea6a52ddc608639e5e06306bce2f973a184fd6
        with self.tqdm_joblib(tqdm(desc='Interferograms', total=len(pairs)*len(subswaths))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(self.intf)(subswath, pair, weight=weight, **kwargs) \
                for (subswath, weight) in zip(subswaths, weights) for pair in pairs)

#         # workaround: start a set of jobs together but not more than available cpu cores at once
#         from joblib.externals import loky
#         if n_jobs == -1:
#             n_jobs = joblib.cpu_count()
#         # create list of arrays [subswath, date1, date2] where all the items are strings
#         subpairs = [[subswath, pair[0], pair[1]] for subswath in subswaths for pair in pairs]
#         n_chunks = int(np.ceil(len(subpairs)/n_jobs))
#         chunks = np.array_split(subpairs, n_chunks)
#         #print ('n_jobs', n_jobs, 'n_chunks', n_chunks, 'chunks', [len(chunk) for chunk in chunks])
#         with tqdm(desc='Interferograms', total=len(subpairs)) as pbar:
#             for chunk in chunks:
#                 loky.get_reusable_executor(kill_workers=True).shutdown(wait=True)
#                 with joblib.parallel_backend('loky', n_jobs=n_jobs):
#                     # convert string subswath to integer value
#                     joblib.Parallel()(joblib.delayed(self.intf)(int(subswath), [date1, date2], **kwargs) \
#                         for (subswath,date1,date2) in chunk)
#                     pbar.update(len(chunk))
