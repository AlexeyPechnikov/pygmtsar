# ----------------------------------------------------------------------------
# PyGMTSAR
#.
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
#.
# Copyright (c) 2023, Alexey Pechnikov
#.
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .datagrid import datagrid
from pygmtsar import tqdm_dask

class IO(datagrid):

    # processing directory
    basedir = '.'

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
        SBAS.restore()

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

    def get_filename(self, name, subswath=None, add_subswath=True):
        import os

        if subswath is None:
            # define subswath when subswath=None and check otherwise
            subswaths = self.get_subswaths()
        else:
            subswaths = [subswath]

        filenames = []
        for swath in subswaths:
            prefix = f'F{swath}_' if add_subswath == True else ''
            filename = os.path.join(self.basedir, f'{prefix}{name}.grd')
            filenames.append(filename)

        # only one filename when only one subswath
        return filenames[0] if subswath is not None or add_subswath is False else filenames

    def get_filenames(self, pairs, name, subswath=None, add_subswath=True):
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

        if subswath is None:
            # define subswath when subswath=None and check otherwise
            subswaths = self.get_subswaths()
        else:
            subswaths = [subswath]

        if isinstance(pairs, pd.DataFrame):
            # convert to standalone DataFrame first
            pairs = self.get_pairs(pairs)[['ref', 'rep']].astype(str).values
        else:
            pairs = np.asarray(pairs)

        filenames_total = []
        for swath in subswaths:
            prefix = f'F{swath}_' if add_subswath == True else ''
            filenames = []
            if len(pairs.shape) == 1:
                # read all the grids from files
                for date in sorted(pairs):
                    filename = os.path.join(self.basedir, f'{prefix}{name}_{date}.grd'.replace('-',''))
                    filenames.append(filename)
            elif len(pairs.shape) == 2:
                # read all the grids from files
                for pair in pairs:
                    filename = os.path.join(self.basedir, f'{prefix}{pair[0]}_{pair[1]}_{name}.grd'.replace('-',''))
                    filenames.append(filename)
            filenames_total.append(filenames)

        return filenames_total if subswath is None else filenames_total[0]

    def load_pairs(self, name='phase', subswath=None):
        import numpy as np
        from glob import glob

        # find all the named grids
        patterns = self.get_filename(f'????????_????????_{name}', subswath)
        if subswath is not None:
            patterns = [patterns]
    
        pairs_total = []
        for pattern in patterns:
            filenames = glob(pattern, recursive=False)
            pairs = [filename.split('_')[-3:-1] for filename in sorted(filenames)]
            # return as numpy array for compatibility reasons
            pairs_total.append(np.asarray(pairs))
    
        return pairs_total if subswath is None else pairs_total[0]

    def open_grid(self, name, subswath=None, add_subswath=True, chunksize=None):
        """
        sbas.open_grid('intf_ll2ra')
        sbas.open_grid('intf_ra2ll')
        sbas.open_grid('intfweight')
        """
        import xarray as xr

        if chunksize is None:
            chunksize = self.chunksize

        if subswath is None:
            # iterate all the subswaths
            subswaths = self.get_subswaths()
        else:
            subswaths = [subswath]

        das = []
        for swath in subswaths:
            filename = self.get_filename(name, swath, add_subswath=add_subswath)
            ds = xr.open_dataset(filename, engine=self.engine, chunks=chunksize)
            if 'a' in ds.dims and 'r' in ds.dims:
                ds = ds.rename({'a': 'y', 'r': 'x'})
            if len(ds.data_vars) == 1:
                das.append(ds[list(ds.data_vars)[0]])
            else:
                das.append(ds)

        return das if subswath is None else das[0]

    def save_grid(self, data, name, subswath=None, caption='Saving 2D grid', chunksize=None):
        import xarray as xr
        import dask
        import os

        if subswath is None:
            subswaths = self.get_subswaths()
        else:
            subswaths = [subswath]
    
        if not isinstance(data, (list, tuple)):
            grids = [data]
        else:
            grids = data
    
        assert len(subswaths) == len(grids), 'ERROR: mismatch between data and subswaths'
    
        if chunksize is None:
            chunksize = self.chunksize

        # save to NetCDF files
        delayeds = []
        for subswath, grid in zip(subswaths, grids):
            filename = self.get_filename(name, subswath)
            if os.path.exists(filename):
                os.remove(filename)

            if isinstance(grid, xr.Dataset):
                encoding = {varname: self.compression(grid[varname].shape, chunksize=chunksize) for varname in grid.data_vars}
            elif isinstance(grid, xr.DataArray):
                encoding = {grid.name: self.compression(grid.shape, chunksize=chunksize)}
            else:
                raise Exception('Argument grid is not xr.Dataset or xr.DataArray object')
            delayed = grid.to_netcdf(filename,
                                  encoding=encoding,
                                  engine=self.engine,
                                  compute=False)
            delayeds.append(delayed)
    
        tqdm_dask(dask.persist(delayeds), desc=caption)
        # cleanup - sometimes writing NetCDF handlers are not closed immediately and block reading access
        import gc; gc.collect()

    def open_stack(self, pairs, name, subswath=None, add_subswath=True, chunksize=None):
        """
        sbas.open_stack(baseline_pairs,'phasefilt')
        """
        import xarray as xr
        import pandas as pd
        import numpy as np
        import os

        if chunksize is None:
            chunksize = self.chunksize

        if subswath is None:
            # iterate all the subswaths
            subswaths = self.get_subswaths() if add_subswath else [None]
        else:
            subswaths = [subswath]

        if isinstance(pairs, pd.DataFrame):
            # convert to standalone DataFrame first
            pairs = self.get_pairs(pairs)[['ref', 'rep']].astype(str).values
        else:
            pairs = np.asarray(pairs)

        def set_pair(ds):
            """
            Extract the pair name from the filename and set it as a coordinate.
            """
            filename = os.path.basename(ds.encoding['source'])
            pair = keys[filename]
            ds = ds.assign_coords(pair=' '.join(pair), ref=pair[0], rep=pair[1])
            if 'a' in ds.dims and 'r' in ds.dims:
                return ds.rename({'a': 'y', 'r': 'x'})
            return ds

        das = []
        for swath in subswaths:
            filenames = self.get_filenames(pairs, name, swath, add_subswath=add_subswath)
            keys = {os.path.basename(filename): pair for filename, pair in zip(filenames, pairs)}

            ds = xr.open_mfdataset(
                filenames,
                engine=self.engine,
                chunks=chunksize,
                parallel=True,
                concat_dim='pair',
                combine='nested',
                preprocess=set_pair
            )

            if len(ds.data_vars) == 1:
                das.append(ds[list(ds.data_vars)[0]])
            else:
                das.append(ds)

        return das if subswath is None else das[0]

    def open_stack_slc(self, dates=None, subswath=None, intensity=False, dfact=2.5e-07, chunksize=None):
        import xarray as xr
        import pandas as pd
        import numpy as np
        import dask

        if chunksize is None:
            chunksize = self.chunksize

        if subswath is None:
            subswaths = self.get_subswaths()
        else:
            subswaths = [subswath]

        stacks = []
        for swath in subswaths:
            if dates is None:
                dates = self.df[self.df['subswath']==swath].index.values
            #print ('dates', dates)
            # select radar coordinates extent
            rng_max = self.PRM(swath).get('num_rng_bins')
            #print ('azi_max', azi_max, 'rng_max', rng_max)
            # use SLC-related chunks for faster processing
            minichunksize = int(np.round(chunksize**2/rng_max))
            #print (minichunksize)
            slcs = [self.PRM(swath, date).read_SLC_int(intensity=intensity, dfact=dfact, chunksize=minichunksize) for date in dates]
            # build stack
            slcs = xr.concat(slcs, dim='date')
            slcs['date'] = pd.to_datetime(dates)
            stacks.append(slcs)

        return stacks if subswath is None else stacks[0]

    def open_stack_geotif(self, dates=None, subswath=None, intensity=False, chunksize=None):
        """
        tiffs = sbas.open_stack_geotif(['2022-06-16', '2022-06-28'], intensity=True)
        """
        import xarray as xr
        import rioxarray as rio
        import pandas as pd
        import numpy as np
        # from GMTSAR code
        DFACT = 2.5e-07
    
        if chunksize is None:
            chunksize = self.chunksize

        if subswath is None:
            subswaths = self.get_subswaths()
        else:
            subswaths = [subswath]

        stack = []
        for swath in self.get_subswaths():
            if dates is None:
                dates = self.df[self.df['subswath']==swath].index.values
            #print ('dates', dates)
            tiffs = self.df[(self.df['subswath']==swath)&(self.df.index.isin(dates))].datapath.values
            tiffs = [rio.open_rasterio(tif, chunks=chunksize)[0] for tif in tiffs]
            # build stack
            tiffs = xr.concat(tiffs, dim='date')
            tiffs['date'] = pd.to_datetime(dates)
            # 2 and 4 multipliers to have the same values as in SLC
            if intensity:
                stack.append((2*DFACT*np.abs(tiffs))**2)
            else:
                stack.append(2*DFACT*tiffs)

        return stacks if subswath is None else stacks[0]

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

        model_filename = self.get_filename(name, add_subswath=False)
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
            assert model.name is not None, 'Define the grid name or use name argument for the NetCDF filename'
            name = model.name
        elif name is None:
            raise ValueError('Specify name for the output NetCDF file')
        
        # save to NetCDF file
        model_filename = self.get_filename(name, add_subswath=False)
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
