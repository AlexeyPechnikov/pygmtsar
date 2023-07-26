# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2021, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .tqdm_joblib import tqdm_joblib
from .tqdm_dask import tqdm_dask
from .datagrid import datagrid

class SBAS_base(tqdm_joblib, datagrid):

    def __repr__(self):
        return 'Object %s %d items\n%r' % (self.__class__.__name__, len(self.df), self.df)

    def to_dataframe(self):
        """
        Return a Pandas DataFrame for all SBAS scenes.

        Returns
        -------
        pandas.DataFrame
            The DataFrame containing SBAS scenes.

        Examples
        --------
        df = sbas.to_dataframe()
        """
        return self.df

    def dump(self, to_path=None):
        """
        Dump SBAS object state to a pickle file (SBAS.pickle in the processing directory by default).

        Parameters
        ----------
        to_path : str, optional
            Path to the output dump file. If not provided, the default dump file in the processing directory is used.

        Returns
        -------
        None

        Examples
        --------
        Dump the current state to the default dump file in the processing directory:
        sbas.dump()

        Notes
        -----
        This method serializes the state of the SBAS object and saves it to a pickle file. The pickle file can be used to
        restore the SBAS object with its processed data and configuration. By default, the dump file is named "SBAS.pickle"
        and is saved in the processing directory. An alternative file path can be provided using the `to_path` parameter.
        """
        import pickle
        import os

        if to_path is None:
            sbas_pickle = os.path.join(self.basedir, 'SBAS.pickle')
        else:
            if os.path.isdir(to_path):
                sbas_pickle = os.path.join(to_path, 'SBAS.pickle')
            else:
                sbas_pickle = to_path
    
        print (f'NOTE: save state to file {sbas_pickle}')
        pickle.dump(self, open(sbas_pickle, 'wb'))

        return

    @staticmethod
    def restore(from_path):
        """
        Restore SBAS object state from a pickle file (SBAS.pickle in the processing directory by default).

        Parameters
        ----------
        from_path : str
            Path to the input dump file.

        Returns
        -------
        SBAS
            The restored SBAS object.

        Examples
        --------
        Restore the current state from the default dump file in the processing directory:
        sbas.restore()

        Notes
        -----
        This static method restores the state of an SBAS object from a pickle file. The pickle file should contain the
        serialized state of the SBAS object, including its processed data and configuration. By default, the method assumes
        the input file is named "SBAS.pickle" and is located in the processing directory. An alternative file path can be
        provided using the `from_path` parameter. The method returns the restored SBAS object.
        """
        import pickle
        import os

        if os.path.isdir(from_path):
            sbas_pickle = os.path.join(from_path, 'SBAS.pickle')
        else:
            sbas_pickle = from_path

        print (f'NOTE: load state from file {sbas_pickle}')
        return pickle.load(open(sbas_pickle, 'rb'))


    def backup(self, backup_dir, copy=False, debug=False):
        """
        Backup framed SBAS scenes, orbits, DEM, and landmask files to build a minimal reproducible dataset.

        Parameters
        ----------
        backup_dir : str
            The backup directory where the files will be copied.
        copy : bool, optional
            Flag indicating whether to make a copy of scene and orbit files and remove the originals in the work directory.
            If False, the files will be moved to the backup directory. Default is False.
        debug : bool, optional
            Flag indicating whether to print debug information. Default is False.

        Returns
        -------
        None

        Examples
        --------
        Backup the files to the specified directory:
        sbas.backup('backup')

        Open the backup for the reproducible run by defining it as a new data directory:
        sbas = SBAS('backup', 'backup/DEM_WGS84.nc', 'raw')

        Notes
        -----
        This method backs up the framed SBAS scenes, orbits, DEM, and landmask files to a specified backup directory.
        It provides a way to create a minimal reproducible dataset by preserving the necessary files for processing.
        The method creates the backup directory if it does not exist. By default, the method moves the scene and orbit files
        to the backup directory, effectively removing them from the work directory. The DEM and landmask files are always
        copied to the backup directory. If the `copy` parameter is set to True, the scene and orbit files will be copied
        instead of moved. Use caution when setting `copy` to True as it can result in duplicated files and consume
        additional storage space. The method also updates the SBAS object's dataframe to mark the removed files as empty.
        """
        import os
        import shutil

        os.makedirs(backup_dir, exist_ok=True)

        # this optional file is dumped state, copy it if exists
        # auto-generated file can't be a symlink but user-defined symlink target should be copied
        filename = os.path.join(self.basedir, 'SBAS.pickle')
        if os.path.exists(filename):
            if debug:
                print ('DEBUG: copy', filename, backup_dir)
            shutil.copy2(filename, backup_dir, follow_symlinks=True)

        # these files required to continue the processing, do not remove and copy only
        filenames = [self.dem_filename, self.landmask_filename]
        for filename in filenames:
            # DEM and landmask can be not defined
            if filename is None:
                continue
            if debug:
                print ('DEBUG: copy', filename, backup_dir)
            shutil.copy2(filename, backup_dir, follow_symlinks=True)

        # these files are large and are not required to continue the processing
        filenames = []
        for record in self.df.itertuples():
            for filename in [record.datapath, record.metapath, record.orbitpath]:
                filenames += filename if isinstance(filename, list) else [filename]
        for filename in filenames:
            # orbit files can be not defined
            if filename is None:
                continue
            # copy and delete the original later to prevent cross-device links issues
            if debug:
                print ('DEBUG: copy', filename, backup_dir)
            shutil.copy2(filename, backup_dir, follow_symlinks=True)
            if not copy and self.basedir == os.path.dirname(filename):
                # when copy is not needed then delete files in work directory only
                if debug:
                    print ('DEBUG: remove', filename)
                os.remove(filename)

        if not copy:
            # mark as empty all the removed files
            for col in ['datapath','metapath','orbitpath']:
                self.df[col] = None
        return

    def multistem_stem(self, subswath, dt=None):
        """
        Define stem and multistem using datetime    
        """
        from datetime import datetime

        # use master datetime if not defined
        if dt is None:
            dt = self.df.loc[self.master, 'datetime']

        stem = f'S1_{dt.strftime("%Y%m%d_%H%M%S")}_F{subswath}'
        multistem = f'S1_{dt.strftime("%Y%m%d")}_ALL_F{subswath}'
        return (multistem, stem)

    def set_master(self, master):
        """
        Define master scene for SBAS object.

        Parameters
        ----------
        master : str
            Date string representing the master scene.

        Returns
        -------
        SBAS
            Modified instance of the SBAS class.

        Examples
        --------
        Set the master scene to '2022-01-20':
        sbas.set_master('2022-01-20')
        """
        if not master in self.df.index:
            raise Exception('Master image not found')
        self.master = master
        return self

    def get_master(self, subswath=None):
        """
        Return dataframe master record(s) for all or only selected subswath.

        Parameters
        ----------
        subswath : str, optional
            The subswath to select. If None, all subswaths are considered. Default is None.

        Returns
        -------
        pd.DataFrame
            The DataFrame containing master records for the specified subswath.
        """
        df = self.df.loc[[self.master]]
        if not subswath is None:
            df = df[df.subswath == subswath]
        assert len(df) > 0, f'Master record for subswath {subswath} not found'
        return df

    def get_aligned(self, subswath=None, date=None):
        """
        Return dataframe aligned records (excluding master) for selected subswath.

        Parameters
        ----------
        subswath : str, optional
            The subswath to select. If None, all subswaths are considered. Default is None.
        date : datetime, optional
            The date for which to return aligned records. If None, all dates are considered. Default is None.

        Returns
        -------
        pd.DataFrame
            The DataFrame containing aligned records for the specified subswath and date.
        """
        if date is None:
            idx = self.df.index.difference([self.master])
        else:
            idx = [date]
        df = self.df.loc[idx]
        if subswath is None:
            return df
        return df[df.subswath == subswath]

    # enlist all the subswaths
    def get_subswaths(self):
        """
        Enlist all the subswaths.

        Returns
        -------
        np.array
            An array containing all unique subswaths.
        """
        import numpy as np
        # note: df.unique() returns unsorted values so it would be 21 instead of expected 12
        return np.unique(self.df.subswath)
    
    def get_subswath(self, subswath=None):
        """
        Check and return subswath or return an unique subswath to functions which work with a single subswath only.

        Parameters
        ----------
        subswath : str, optional
            The subswath to check. If None, an unique subswath is returned. Default is None.

        Returns
        -------
        str
            The checked or unique subswath.
        """
        # detect all the subswaths
        subswaths = self.get_subswaths()
        assert subswath is None or subswath in subswaths, f'ERROR: subswath {subswath} not found'
        if subswath is not None:
            return subswath
        assert len(subswaths)==1, f'ERROR: multiple subswaths {subswaths} found, merge them first using SBAS.merge_parallel()'
        # define subswath
        return subswaths[0]

    # function is obsolete
    def find_pairs(self, name='phasefilt'):
        """
        Find pairs. This function is obsolete. Use SBAS.pairs() function instead.

        Parameters
        ----------
        name : str, optional
            The name of the phase filter. Default is 'phasefilt'.

        Returns
        -------
        np.ndarray
            An array of pairs.
        """
        print ('NOTE: use SBAS.pairs() wrapper function to get pairs as DataFrame and optionally dates array')        
        pairs = self.pairs(name=name)
        return pairs

    # function is obsolete
    def find_dates(self, pairs=None):
        """
        Find dates. This function is obsolete. Use SBAS.pairs() function instead.

        Parameters
        ----------
        pairs : np.ndarray, optional
            An array of pairs. If None, all pairs are considered. Default is None.

        Returns
        -------
        np.ndarray
            An array of dates.
        """
        print ('NOTE: use SBAS.pairs() wrapper function to get pairs as DataFrame and optionally dates array')
        dates = self.pairs(pairs, dates=True)[1]
        return dates

    def pairs(self, pairs=None, dates=False, name='phasefilt'):
        """
        Get pairs as DataFrame and optionally dates array.

        Parameters
        ----------
        pairs : np.ndarray, optional
            An array of pairs. If None, all pairs are considered. Default is None.
        dates : bool, optional
            Whether to return dates array. Default is False.
        name : str, optional
            The name of the phase filter. Default is 'phasefilt'.

        Returns
        -------
        pd.DataFrame or tuple
            A DataFrame of pairs. If dates is True, also returns an array of dates.
        """
        import pandas as pd
        import numpy as np
        from glob import glob

        if pairs is None:
            # find all the named grids
            pattern = self.get_filenames(None, None, f'????????_????????_{name}')
            filenames = glob(pattern, recursive=False)
            pairs = [filename.split('_')[-3:-1] for filename in sorted(filenames)]
            # return as numpy array for compatibility reasons
            pairs = np.asarray(pairs)
            # check that all the pairs produced from the SBAS scenes
            #dates = list(map(lambda x: x.replace('-',''), self.df.index))
            #invalid = [pair for pair in pairs.flatten() if pair not in dates]
            #assert len(invalid) == 0, 'ERROR: found grids for pairs not in the SBAS scenes. Define valid pairs manually.'

        if not isinstance(pairs, pd.DataFrame):
            # Convert numpy array to DataFrame
            pairs = pd.DataFrame(pairs, columns=['ref', 'rep'])
            # Convert ref and rep columns to datetime format
            pairs['ref'] = pd.to_datetime(pairs['ref'])
            pairs['rep'] = pd.to_datetime(pairs['rep'])
            # Calculate the duration in days and add it as a new column
            pairs['duration'] = (pairs['rep'] - pairs['ref']).dt.days
        else:
            # workaround for baseline_pairs() output
            pairs = pairs.rename(columns={'ref_date': 'ref', 'rep_date': 'rep'})

        if dates:
            # pairs is DataFrame
            dates = np.unique(pairs[['ref', 'rep']].astype(str).values.flatten())
            return (pairs, dates)
        return pairs

    # use the function for open_grids() and save_grids()
    def get_filenames(self, subswath, pairs, name, add_subswath=True):
        """
        Get the filenames of the data grids. The filenames are determined by the subswath, pairs, and name parameters.

        Parameters
        ----------
        subswath : int or None
            The subswath number. If None, the function will determine it based on the dataset.
        pairs : np.ndarray or pd.DataFrame or None
            An array or DataFrame of pairs. If None, the function will open a single grid or a set of subswath grids.
        name : str
            The name of the grid to be opened.
        add_subswath : bool, optional
            Whether to add subswath to the prefix of the filename. Default is True.

        Returns
        -------
        str or list of str
            The filename or a list of filenames of the grids.
        """
        import pandas as pd
        import numpy as np
        import os

        # define subswath when subswath=None and check otherwise
        subswath = self.get_subswath(subswath)
    
        if pairs is None:
            # special case to open a single grid {name}.grd or a set of subswath grids Fn_{name}.grd
            pass
        elif isinstance(pairs, pd.DataFrame):
            # convert to standalone DataFrame first
            pairs = self.pairs(pairs)[['ref', 'rep']].astype(str).values
        else:
            pairs = np.asarray(pairs)

        if add_subswath == True:
            prefix = f'F{subswath}_'
        else:
            prefix = ''

        filenames = []
        if pairs is None:
            # special case to open a single grid {name}.grd or a set of subswath grids Fn_{name}.grd
            filename = os.path.join(self.basedir, f'{prefix}{name}.grd')
            return filename
        elif len(pairs.shape) == 1:
            # read all the grids from files
            for date in sorted(pairs):
                filename = os.path.join(self.basedir, f'{prefix}{name}_{date}.grd'.replace('-',''))
                filenames.append(filename)
        elif len(pairs.shape) == 2:
            # read all the grids from files
            for pair in pairs:
                filename = os.path.join(self.basedir, f'{prefix}{pair[0]}_{pair[1]}_{name}.grd'.replace('-',''))
                filenames.append(filename)
        return filenames

    def save_grids(self, grids, name, func=None, add_subswath=True, chunksize=None, n_jobs=1, interactive=True, **kwargs):
        """
        Save the input grids with the given name. The grids are preprocessed by a given function before saving.

        Parameters
        ----------
        grids : xr.DataArray
            The input data grids to be saved.
        name : str
            The name of the grid to be saved.
        func : callable, list of callable or None, optional
            Function or list of functions to be applied to each grid before saving. 
            If None, no function is applied. Default is None.
        add_subswath : bool, optional
            Whether to add subswath to the prefix of the filename. Default is True.
        chunksize : int or None, optional
            The chunk size to be used when saving the grid. If None, the chunk size is determined automatically. Default is None.
        n_jobs : int, optional
            The number of jobs to run in parallel for saving the grids. Default is 1.
        interactive : bool, optional
            Whether to display a progress bar during the saving process. Default is True.
        **kwargs
            Additional keyword arguments to be passed to the saving function.

        Returns
        -------
        None
        """
        import xarray as xr
        import numpy as np
        from tqdm.auto import tqdm
        import joblib
        import os

        subswaths = self.get_subswaths()
        # pipeline before merging subswaths is well defined and we do not need to use this function
        assert len(subswaths) == 1, 'ERROR: use SBAS.merge() to merge multiple subswaths before'

        if isinstance(grids, xr.DataArray):
            if len(grids.dims) == 3:
                grids = [da for da in grids]
            elif len(grids.dims) == 2:
                # special case for a single 2D grid
                filename = self.get_filenames(None, None, name, add_subswath=add_subswath)
                grids.astype(np.float32).rename(name)\
                    .to_netcdf(filename, encoding={name: self.compression(grids.shape, chunksize=chunksize)}, engine=self.engine)
                return
            else:
                assert 0, 'ERROR: supported 2D and 3D arrays only'

        def preprocess(subswath, da):
            if func is not None:
                if isinstance(func, list):
                    for f in func:
                        da = f(da)
                else:
                    da = func(da)
            # check 3rd dimension
            if 'date' in da[0].coords:
                # "date" for SBAS results for an example
                pair = da.date.item()
            elif 'pair' in grids[0].coords:
                # "pair" for interferograms for an example
                pair = da.pair.item().split(' ')
            # the function returns a set of filenames
            filename = self.get_filenames(subswath, [pair], name, add_subswath=add_subswath)[0]
            #print ('filename', filename)
            if os.path.exists(filename):
                os.remove(filename)
            da.astype(np.float32).rename(name)\
                .to_netcdf(filename, encoding={name: self.compression(da.shape, chunksize=chunksize)}, engine=self.engine)
            return

        # process all the grids
        if interactive:
            with self.tqdm_joblib(tqdm(desc='Saving', total=len(grids)*len(subswaths))) as progress_bar:
                joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(preprocess)(None, grid) for grid in grids)
        else:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(preprocess)(None, grid) for grid in grids)

    # returns all grids in basedir by mask or grids by dates and name
    # Backward-compatible open_grids() returns list of grids fot the name or a single grid for a single subswath
    def open_grids(self, pairs, name, geocode=False, inverse_geocode=False,  mask=None, func=None,
                   crop_valid=False, add_subswath=True, chunksize=None, n_jobs=-1, interactive=True):
        """
        Lazy open PyGMTSAR produced NetCDF grids as 3D data cube and apply a set of functions on it.

        Parameters
        ----------
        pairs : list, pandas.DataFrame, or numpy.ndarray
            List of date pairs (baseline pairs) or dates.
        name : str
            Grid name: 'phasefilt', 'corr', 'unwrap', 'detrend', 'disp'.
        geocode : bool, optional
            Whether to geocode the grid to geographic coordinates. Default is False.
        inverse_geocode : bool, optional
            Whether to inverse geocode the grid to radar coordinates. Default is False.
        mask : xarray.DataArray, optional
            Mask to exclude invalid areas. Default is None.
        func : function or list of functions, optional
            Function(s) to apply to each input grid. Default is None.
        crop_valid : bool, optional
            Whether to crop NaN values in the grids. Default is False.
        add_subswath : bool, optional
            Whether to add subswath name as 'Fn' to grid names. Default is True.
        chunksize : int, optional
            Size of the chunks for processing the data instead of the default value. Default is None.
        n_jobs : int, optional
            Number of CPU cores to use for parallel processing. Default is -1, which means using all available cores.
        interactive : bool, optional
            If True, show the progress indicator. Default is True.

        Returns
        -------
        3D Xarray DataArray or list of 3D Xarray DataArray
            The opened and processed grids as a 3D Xarray DataArray or a list of 3D Xarray DataArray.

        Examples
        --------
        Open unwrapped phase grids in radar coordinates:
        unwraps_ra = sbas.open_grids(pairs, 'unwrap')

        Open PyGMTSAR SBAS grids in radar coordinates and calculate LOS displacement and fill NODATA areas:
        dates = np.unique(pairs.values.flatten() if isinstance(pairs, pd.DataFrame) else pairs.flatten())
        disps = sbas.open_grids(dates, 'disp', func=[sbas.los_displacement_mm, sbas.nearest_grid])

        Open GMTSAR SBAS grids in radar coordinates using compatibility option add_subswath=False:
        sbas.open_grids(sbas.df.index, 'disp', func=sbas.nearest_grid, add_subswath=False)

        Calculate LOS displacement for unwrapped phase grids in radar coordinates:
        unwraps_ra = sbas.open_grids(pairs, 'unwrap')
        los_disp_ra = sbas.los_displacement_mm(unwraps_ra)
        # or the same code in one line
        los_disp_ra = sbas.open_grids(pairs, 'unwrap', func=sbas.los_displacement_mm)
        # Note: here "func" argument for open_grids() function reduces the code to a single command.

        Calculate LOS displacement for detrended unwrapped phase grids in geographic coordinates:
        detrend_ll = sbas.open_grids(pairs, 'detrend', geocode=True)
        los_disp_ll = sbas.los_displacement_mm(detrend_ll)
        # or the same code in one line
        los_disp_ll = sbas.open_grids(pairs, 'detrend', geocode=True, func=sbas.los_displacement_mm)
        # Note: here "func" argument for open_grids() function reduces the code to a single command.

        Open GMTSAR SBAS velocity grid in radar coordinates:
        vel = sbas.open_grids(None, 'vel', add_subswath=False)
        Open GMTSAR SBAS velocity grid in geographic coordinates:
        vel = sbas.open_grids(None, 'vel', geocode=True, add_subswath=False)

        Notes
        -----
        This method lazily opens PyGMTSAR-produced NetCDF grids as a 3D data cube and applies a set of functions on it.
        It can open multiple subswaths of grids and process them in parallel using Dask and Joblib.
        The grids can be geocoded to geographic coordinates or inverse geocoded to radar coordinates.
        Various processing options such as applying functions, cropping NaN values, and adding subswath names are available.
        """
        import pandas as pd
        import xarray as xr
        import numpy as np
        from tqdm.auto import tqdm
        import joblib

        assert not(geocode and inverse_geocode), 'ERROR: Only single geocoding option can be applied'

        if chunksize is None:
            chunksize = self.chunksize

        # iterate all the subswaths
        subswaths = self.get_subswaths()

        # for backward compatibility
        if isinstance(mask, bool):
            print ('NOTE: mask argument changed from boolean to dataarray for SBAS.open_grid() function call')
            mask = None
        if mask is not None:
            assert len(subswaths) == 1, 'ERROR: mask can be applied to a single subswath only'
            nanmask = xr.where((mask == 0)|(np.isnan(mask)), np.nan, 1)

        if pairs is None:
            # special case to open a single grid {name}.grd or a set of subswath grids Fn_{name}.grd
            pass
        elif isinstance(pairs, pd.DataFrame):
            # convert to standalone DataFrame first
            pairs = self.pairs(pairs)[['ref', 'rep']].astype(str).values
        else:
            pairs = np.asarray(pairs)

        def postprocess(da, subswath):
            masked = False
            if mask is not None and self.is_same(mask, da):
                # apply mask when the mask defined in the same coordinates only
                if geocode or inverse_geocode: 
                    # apply the mask later when the coordinates are not transformed
                    da = nanmask * da
                    # apply the mask just once
                    masked = True
            if geocode:
                assert self.is_ra(da), 'ERROR: geocode option requires radar coordinates grid'
                da = self.intf_ra2ll(da, chunksize=chunksize)
            elif inverse_geocode:
                assert self.is_geo(da), 'ERROR: inverse_geocode option requires geographic coordinates grid'
                da = self.intf_ll2ra(da, chunksize=chunksize)
            # apply user-defined post-processing function(s)
            if func is not None:
                if isinstance(func, list):
                    for f in func:
                        da = f(da)
                else:
                    da = func(da)
            # second masking can be used for some cases like to fill low-coherence areas and crop valid area only
            if mask is not None and self.is_same(mask, da):
                # apply mask when the mask defined in the same coordinates only
                da = nanmask * da
                masked = True
            if not masked and mask is not None:
                print ('NOTE: the mask is not applied because it does not correspond to the grid coordinates (radar or geographic)')
            return da

        # use chunksize variable
        def open_grid(filename):
            assert chunksize is not None, 'open_grids() chunksize is not defined'
            da = xr.open_dataarray(filename, engine=self.engine, chunks=chunksize)
            # workaround for Google Colab when we cannot save grids with x,y coordinate names
            if 'a' in da.dims:
                da = da.rename({'a': 'y'})
            if 'r' in da.dims:
                da = da.rename({'r': 'x'})
            return da

        dass = []
        for subswath in subswaths:
            filenames = self.get_filenames(subswath, pairs=pairs, name=name, add_subswath=add_subswath)

            das = []
            if pairs is None:
                # special case to open a single grid {name}.grd or a set of subswath grids Fn_{name}.grd
                #print ('filename', filename)
                da = open_grid(filenames)
                das  = postprocess(da, subswath)
            elif len(pairs.shape) == 1:
                # read all the grids from files
                for filename in filenames:
                    #print (date, filename)
                    da = open_grid(filename)
                    das.append(da)

                # post-processing on a set of 2D rasters
                if interactive:
                    with self.tqdm_joblib(tqdm(desc='Loading', total=len(das))) as progress_bar:
                        das = joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(postprocess)(da, subswath) for da in das)
                else:
                    das = joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(postprocess)(da, subswath) for da in das)

                # prepare to stacking
                das = [da.expand_dims('date') for da in das]

                # allow stack to be extended to largest 1st dimension size
                # to be sure all code work well for this case
                # so user is able to load grids by his own way
                das = xr.concat(das, dim='date')
                das['date'] = sorted(pairs)
            elif len(pairs.shape) == 2:
                # read all the grids from files
                for filename in filenames:
                    #print (filename)
                    da = open_grid(filename)
                    das.append(da)

                # post-processing on a set of 2D rasters
                if interactive:
                    with self.tqdm_joblib(tqdm(desc='Loading', total=len(das))) as progress_bar:
                        das = joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(postprocess)(da, subswath) for da in das)
                else:
                    das = joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(postprocess)(da, subswath) for da in das)

                # prepare to stacking
                das = [da.expand_dims('pair') for da in das]

                # allow stack to be extended to largest 1st dimension size
                # to be sure all code work well for this case
                # so user is able to load grids by his own way
                das = xr.concat(das, dim='pair')
                das['pair'] = [f'{pair[0]} {pair[1]}' for pair in pairs]
                das['ref']  = xr.DataArray([pair[0] for pair in pairs], dims='pair')
                das['rep']  = xr.DataArray([pair[1] for pair in pairs], dims='pair')
            else:
                raise Exception('Use single or two columns Pandas dataset or array as "pairs" argument')

            # crop NaNs
            if crop_valid:
                das = self.cropna(das)
            dass.append(das)

        return dass[0] if len(dass) == 1 else dass

    def open_model(self, name, chunksize=None):
        """
        Opens an xarray 3D Dataset from a NetCDF file and re-chunks it based on the specified chunksize.

        This function takes the name of the model to be opened, reads the NetCDF file, and re-chunks
        the dataset according to the provided chunksize or the default value from the 'sbas' object.
        The 'date' dimension is always chunked with a size of 1.

        Parameters
        ----------
        name : str
            The name of the model file to be opened.
        chunksize : int, optional
            The chunk size to be used for dimensions other than 'date'. If not provided, the default
            chunk size from the 'sbas' object will be used.

        Returns
        -------
        xarray.Dataset
            The re-chunked xarray Dataset read from the specified NetCDF file.

        """
        import xarray as xr
        import pandas as pd
        import os

        if chunksize is None:
            chunksize = self.chunksize

        model_filename = self.get_filenames(None, None, name, add_subswath=False)
        assert os.path.exists(model_filename), f'ERROR: The NetCDF file is missed: {model_filename}'

        # Open the dataset without chunking
        model = xr.open_dataset(model_filename, engine=self.engine)
        chunks = {dim: 1 if dim == 'date' else chunksize for dim in model.dims}
        # Re-chunk the dataset using the chunks dictionary
        model = model.chunk(chunks)

        # convert string dates to dates
        if 'date' in model.dims:
            model['date'] = pd.to_datetime(model['date'].values)
    
        if 'yy' in model.dims and 'xx' in model.dims:
            model = model.rename({'yy': 'lat', 'xx': 'lon'})

        if 'a' in model.dims and 'r' in model.dims:
            model = model.rename({'a': 'y', 'r': 'x'})

        # in case of a single variable return DataArray
        if len(model.data_vars) == 1:
            return model[list(model.data_vars)[0]]
        return model

    def save_model(self, model, name=None, caption='Saving 3D datacube', chunksize=None, debug=False):
        """
        Save an xarray 3D Dataset to a NetCDF file and re-chunks it based on the specified chunksize.

        The 'date' dimension is always chunked with a size of 1.

        Parameters
        ----------
        model : xarray.Dataset
            The model to be saved.
        chunksize : int, optional
            The chunk size to be used for dimensions other than 'date'. If not provided, the default
            chunk size from the 'sbas' object will be used.
        caption: str
            The text caption for the saving progress bar.

        Returns
        -------
        None
        """
        import xarray as xr
        import dask
        import os

        if chunksize is None:
            chunksize = self.chunksize

        if name is None and isinstance(model, xr.DataArray):
            name = model.name
        elif name is None:
            raise ValueError('Specify name for the output NetCDF file')
        
        # save to NetCDF file
        model_filename = self.get_filenames(None, None, name, add_subswath=False)
        if os.path.exists(model_filename):
            os.remove(model_filename)
        if isinstance(model, xr.DataArray):
            if debug:
                print ('DEBUG: DataArray')
            netcdf_compression = self.compression(model.shape, chunksize=(1, chunksize, chunksize))
            encoding = {model.name: netcdf_compression}
        elif isinstance(model, xr.Dataset):
            if debug:
                print ('DEBUG: Dataset')
            if len(model.dims) == 3:
                encoding = {varname: self.compression(model[varname].shape, chunksize=(1, chunksize, chunksize)) for varname in model.data_vars}
            else:
                encoding = {varname: self.compression(model[varname].shape, chunksize=chunksize) for varname in model.data_vars}
        delayed = model.to_netcdf(model_filename,
                         engine=self.engine,
                         encoding=encoding,
                         compute=False)
        tqdm_dask(dask.persist(delayed), desc=caption)

        # cleanup - sometimes writing NetCDF handlers are not closed immediately and block reading access
        import gc; gc.collect()
