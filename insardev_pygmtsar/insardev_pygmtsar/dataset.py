# ----------------------------------------------------------------------------
# InSAR.dev
# 
# This file is part of the InSAR.dev project: https://InSAR.dev
# 
# Copyright (c) 2025, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from insardev_toolkit import datagrid, tqdm_dask

class dataset(datagrid):

    # redefine
    netcdf_complevel = 1
    
    # work directory
    basedir = '.'

    def _glob_re(self, name, basedir='auto'):
        import os
        import re

        if isinstance(basedir, str) and basedir == 'auto':
            basedir = self.basedir

        filenames = filter(re.compile(name).match, os.listdir(basedir))
        return sorted([os.path.join(basedir, filename) for filename in filenames])

    def get_filename(self, name, basedir='auto'):
        import os

        if isinstance(basedir, str) and basedir == 'auto':
            basedir = self.basedir

        filename = os.path.join(basedir, f'{name}.nc')
        return filename


    def get_filenames(self, pairs, name, basedir='auto'):
        """
        Get the filenames of the data grids. The filenames are determined by the pairs and name parameters.

        Parameters
        ----------
        pairs : np.ndarray or pd.DataFrame or None
            An array or DataFrame of pairs. If None, the function will open a single grid or a set of subswath grids.
        name : str
            The name of the grid to be opened.
        
        Returns
        -------
        str or list of str
            The filename or a list of filenames of the grids.
        """
        import pandas as pd
        import numpy as np
        import os

        if isinstance(basedir, str) and basedir == 'auto':
            basedir = self.basedir

        if isinstance(pairs, pd.DataFrame):
            # convert to standalone DataFrame first
            pairs = self.get_pairs(pairs)[['ref', 'rep']].astype(str).values
        else:
            pairs = np.asarray(pairs)

        if name == '' or name is None:
            name = ''
        else:
            name = name + '_'
 
        filenames = []
        if len(pairs.shape) == 1:
            # read all the grids from files
            for date in sorted(pairs):
                filename = os.path.join(basedir, f'{name}{date}.nc'.replace('-',''))
                filenames.append(filename)
        elif len(pairs.shape) == 2:
            # read all the grids from files
            for pair in pairs:
                filename = os.path.join(basedir, f'{name}{pair[0]}_{pair[1]}.nc'.replace('-',''))
                filenames.append(filename)
        return filenames

    def open_cube(self, name, basedir='auto'):
        """
        Opens an xarray 2D/3D Dataset or dataArray from a NetCDF file.

        This function takes the name of the model to be opened, reads the NetCDF file, and re-chunks
        the dataset according to the provided chunksize or the default value from the 'stack' object.
        The 'date' dimension is always chunked with a size of 1.

        Parameters
        ----------
        name : str
            The name of the model file to be opened.

        Returns
        -------
        xarray.Dataset
            Xarray Dataset read from the specified NetCDF file.

        """
        import xarray as xr
        import pandas as pd
        import numpy as np
        import os

        filename = self.get_filename(name, basedir=basedir)
        assert os.path.exists(filename), f'ERROR: The NetCDF file is missed: {filename}'

        # Workaround: open the dataset without chunking
        data = xr.open_dataset(filename,
                               engine=self.netcdf_engine_read,
                               format=self.netcdf_format)
        
        if 'stack' in data.dims:
            if 'y' in data.coords and 'x' in data.coords:
                multi_index_names = ['y', 'x']
            elif 'lat' in data.coords and 'lon' in data.coords:
                multi_index_names = ['lat', 'lon']
            multi_index = pd.MultiIndex.from_arrays([data.y.values, data.x.values], names=multi_index_names)
            data = data.assign_coords(stack=multi_index).set_index({'stack': ['y', 'x']})
            chunksize = self.chunksize1d
        else:
            chunksize = self.chunksize

        # set the proper chunk sizes
        chunks = {dim: 1 if dim in ['pair', 'date'] else chunksize for dim in data.dims}
        data = data.chunk(chunks)

        # attributes are empty when dataarray is prezented as dataset
        # revert dataarray converted to dataset
        data_vars = list(data.data_vars)
        if len(data_vars) == 1 and 'dataarray' in data.attrs:
            assert data.attrs['dataarray'] == data_vars[0]
            data = data[data_vars[0]]

        # convert string dates to dates
        for dim in ['date', 'ref', 'rep']:
            if dim in data.dims:
                data[dim] = pd.to_datetime(data[dim])

        return data

    def sync_cube(self, data, name=None, spatial_ref=None, caption='Syncing NetCDF 2D/3D Dataset', basedir='auto'):
        import xarray as xr
        if name is None and isinstance(data, xr.DataArray):
            assert data.name is not None, 'Define data name or use "name" argument for the NetCDF filename'
            name = data.name
        elif name is None:
            raise ValueError('Specify name for the output NetCDF file')
        self.save_cube(data, name, spatial_ref, caption, basedir)
        return self.open_cube(name, basedir)

    def save_cube(self, data, name=None, spatial_ref=None, caption='Saving NetCDF 2D/3D Dataset', basedir='auto'):
        """
        Save a lazy and not lazy 2D/3D xarray Dataset or DataArray to a NetCDF file.

        The 'date' or 'pair' dimension is always chunked with a size of 1.

        Parameters
        ----------
        data : xarray.Dataset or xarray.DataArray
            The model to be saved.
        name : str
            The text name for the output NetCDF file.
        caption: str
            The text caption for the saving progress bar.

        Returns
        -------
        None

        Examples
        -------
        stack.save_cube(intf90m, 'intf90m')                              # save lazy 3d dataset
        stack.save_cube(intf90m.phase, 'intf90m')                        # save lazy 3d dataarray
        stack.save_cube(intf90m.isel(pair=0), 'intf90m')                 # save lazy 2d dataset
        stack.save_cube(intf90m.isel(pair=0).phase, 'intf90m')           # save lazy 2d dataarray
        stack.save_cube(intf90m.compute(), 'intf90m')                    # save 3d dataset     
        stack.save_cube(intf90m.phase.compute(), 'intf90m')              # save 3d dataarray
        stack.save_cube(intf90m.isel(pair=0).compute(), 'intf90m')       # save 2d dataset
        stack.save_cube(intf90m.isel(pair=0).phase.compute(), 'intf90m') # save 2d dataarray
        """
        import xarray as xr
        import pandas as pd
        import dask
        import os
        import warnings
        # suppress Dask warning "RuntimeWarning: invalid value encountered in divide"
        warnings.filterwarnings('ignore')
        warnings.filterwarnings('ignore', module='dask')
        warnings.filterwarnings('ignore', module='dask.core')
        import logging
        # prevent warnings "RuntimeWarning: All-NaN slice encountered"
        logging.getLogger('distributed.nanny').setLevel(logging.ERROR)
        # disable "distributed.utils_perf - WARNING - full garbage collections ..."
        try:
            from dask.distributed import utils_perf
            utils_perf.disable_gc_diagnosis()
        except ImportError:
            from distributed.gc import disable_gc_diagnosis
            disable_gc_diagnosis()

        if name is None and isinstance(data, xr.DataArray):
            assert data.name is not None, 'Define data name or use "name" argument for the NetCDF filename'
            name = data.name
        elif name is None:
            raise ValueError('Specify name for the output NetCDF file')

        chunksize = None
        if 'stack' in data.dims and isinstance(data.coords['stack'].to_index(), pd.MultiIndex):
            # replace multiindex by sequential numbers 0,1,...
            data = data.reset_index('stack')
            # single-dimensional data compression required
            chunksize = self.netcdf_chunksize1d

        if isinstance(data, xr.DataArray):
            if data.name is None:
                data = data.rename(name)
            data = data.to_dataset().assign_attrs({'dataarray': data.name})

        is_dask = isinstance(data[list(data.data_vars)[0]].data, dask.array.Array)
        encoding = {varname: self._compression(data[varname].shape, chunksize=chunksize) for varname in data.data_vars}
        #print ('save_cube encoding', encoding)
        #print ('is_dask', is_dask, 'encoding', encoding)

        # save to NetCDF file
        filename = self.get_filename(name, basedir=basedir)
        if os.path.exists(filename):
            os.remove(filename)
        delayed = self.spatial_ref(data, spatial_ref).to_netcdf(filename,
                                 engine=self.netcdf_engine_write,
                                 format=self.netcdf_format,
                                 encoding=encoding,
                                 compute=not is_dask)
        if is_dask:
            tqdm_dask(result := dask.persist(delayed), desc=caption)
            # cleanup - sometimes writing NetCDF handlers are not closed immediately and block reading access
            del delayed, result
            import gc; gc.collect()

    def delete_cube(self, name, basedir='auto'):
        import os

        filename = self.get_filename(name, basedir=basedir)
        #print ('filename', filename)
        if os.path.exists(filename):
            os.remove(filename)

    def sync_stack(self, data, name=None, spatial_ref=None, caption='Saving 2D Stack', basedir='auto', queue=None, timeout=300):
        import xarray as xr
        if name is None and isinstance(data, xr.DataArray):
            assert data.name is not None, 'Define data name or use "name" argument for the NetCDF filenames'
            name = data.name
        elif name is None:
            raise ValueError('Specify name for the output NetCDF files')
        self.delete_stack(name, basedir=basedir)
        self.save_stack(data, name, spatial_ref, caption, basedir, queue, timeout)
        return self.open_stack(name, basedir)

    def open_stack(self, name, stack=None, basedir='auto'):
        """
        Examples:
        stack.open_stack('data')
        stack.open_stack('data', ['2018-03-23'])
        stack.open_stack('data', ['2018-03-23', '2018-03-11'])
        stack.open_stack('phase15m')
        stack.open_stack('intf90m',[['2018-02-21','2018-03-11']])
        stack.open_stack('intf90m', stack.get_pairs([['2018-02-21','2018-03-11']]))
        """
        import xarray as xr
        import pandas as pd
        import numpy as np
        import glob
        import os
    
        if name == '' or name is None:
            name = ''
        else:
            name = name + '_'

        if stack is None:
            # look for all stack files
            #filenames = self.get_filenames(['*'], name)[0]
            #filenames = self.get_filename(f'{name}_????????_????????')
            # like data_20180323.nc or intf60m_20230114_20230219.nc
            filenames = self._glob_re(name + '[0-9]{8}(_[0-9]{8})*.nc', basedir=basedir)
        elif isinstance(stack, (list, tuple, np.ndarray)) and len(np.asarray(stack).shape) == 1:
            # dates
            filenames = self.get_filenames(np.asarray(stack), name, basedir=basedir)
        else:
            # pairs
            filenames = self.get_filenames(stack, name, basedir=basedir)
        #print ('filenames', filenames)

        data = xr.open_mfdataset(
            filenames,
            engine=self.netcdf_engine_read,
            format=self.netcdf_format,
            parallel=True,
            concat_dim='stackvar',
            chunks={"stackvar": 1},
            combine='nested'
        )
        
        if 'stack' in data.dims:
            if 'y' in data.coords and 'x' in data.coords:
                multi_index_names = ['y', 'x']
            elif 'lat' in data.coords and 'lon' in data.coords:
                multi_index_names = ['lat', 'lon']
            multi_index = pd.MultiIndex.from_arrays([data.y.values, data.x.values], names=multi_index_names)
            data = data.assign_coords(stack=multi_index).set_index({'stack': ['y', 'x']}).chunk({'stack': self.chunksize1d})
        else:
            dims = list(data.dims)
            data = data.chunk({dims[0]: 1, dims[1]: self.chunksize, dims[2]: self.chunksize})

        # revert dataarray converted to dataset
        data_vars = list(data.data_vars)
        if 'dataarray' in data.attrs:
            data = data[data.attrs['dataarray']]

        for dim in ['pair', 'date']:
            #if dim in data.coords:
            if dim in (data.data_vars if isinstance(data, xr.Dataset) else data.coords):
                if data[dim].shape == () or 'stack' in data.dims:
                    if data[dim].shape == ():
                        data = data.assign_coords(pair=('stackvar', [data[dim].values]))
                    data = data.rename({'stackvar': dim}).set_index({dim: dim})
                else:
                    data = data.swap_dims({'stackvar': dim})
    
        # convert string (or already timestamp) dates to dates
        for dim in ['date', 'ref', 'rep']:
            if dim in data.dims:
                if not data[dim].shape == ():
                    data[dim] = pd.to_datetime(data[dim])
                else:
                    data[dim].values = pd.to_datetime(data['date'].values)
    
        return data

    # use save_mfdataset
    def save_stack(self, data, name, spatial_ref=None, caption='Saving 2D Stack', basedir='auto', queue=None, timeout=None):
        import numpy as np
        import xarray as xr
        import pandas as pd
        import dask
        import os
        from dask.distributed import get_client
        import warnings
        # suppress Dask warning "RuntimeWarning: invalid value encountered in divide"
        warnings.filterwarnings('ignore')
        warnings.filterwarnings('ignore', module='dask')
        warnings.filterwarnings('ignore', module='dask.core')
        # Filter out Dask "Restarting worker" warnings
        warnings.filterwarnings("ignore", module="distributed.nanny")
        import logging
        # Suppress Dask "Restarting worker" warnings
        logging.getLogger('distributed.nanny').setLevel(logging.ERROR)
        # disable "distributed.utils_perf - WARNING - full garbage collections ..."
        try:
            from dask.distributed import utils_perf
            utils_perf.disable_gc_diagnosis()
        except ImportError:
            from distributed.gc import disable_gc_diagnosis
            disable_gc_diagnosis()
    
        # Dask cluster client
        client = get_client()
        
        if isinstance(data, xr.Dataset):
            stackvar = data[list(data.data_vars)[0]].dims[0]
            is_dask = isinstance(data[list(data.data_vars)[0]].data, dask.array.Array)
        elif isinstance(data, xr.DataArray):
            stackvar = data.dims[0]
            is_dask = isinstance(data.data, dask.array.Array)
        else:
            raise Exception('Argument grid is not xr.Dataset or xr.DataArray object')
        #print ('is_dask', is_dask, 'stackvar', stackvar)
        stacksize = data[stackvar].size
    
        if queue is None:
            queue = self.netcdf_queue
        if queue is None:
            # process all the stack items in a single operation
            queue = stacksize
    
        if 'stack' in data.dims and isinstance(data.coords['stack'].to_index(), pd.MultiIndex):
            # replace multiindex by sequential numbers 0,1,...
            data = data.reset_index('stack')

        if isinstance(data, xr.DataArray):
            data = data.to_dataset().assign_attrs({'dataarray': data.name})
        encoding = {varname: self._compression(data[varname].shape[1:]) for varname in data.data_vars}
        #print ('save_stack encoding', encoding)
    
        # Applying iterative processing to prevent Dask scheduler deadlocks.
        counter = 0
        digits = len(str(stacksize))
        # Splitting all the pairs into chunks, each containing approximately queue pairs.
        n_chunks = stacksize // queue if stacksize > queue else 1
        for chunk in np.array_split(range(stacksize), n_chunks):
            dss = [data.isel({stackvar: ind}) for ind in chunk]
            if stackvar == 'date':
                stackvals = [ds[stackvar].dt.date.values for ds in dss]
            else:
                stackvals = [ds[stackvar].item().split(' ') for ds in dss]
            # save to NetCDF file
            filenames = self.get_filenames(stackvals, name, basedir=basedir)
            #[os.remove(filename) for filename in filenames if os.path.exists(filename)]
            delayeds = xr.save_mfdataset(self.spatial_ref(dss, spatial_ref),
                                         filenames,
                                         encoding=encoding,
                                         engine=self.netcdf_engine_write,
                                         format=self.netcdf_format,
                                         compute=not is_dask)
            # process lazy chunk
            if is_dask:
                if n_chunks > 1:
                    chunk_caption = f'{caption}: {(counter+1):0{digits}}...{(counter+len(chunk)):0{digits}} from {stacksize}'
                else:
                    chunk_caption = caption
                tqdm_dask(result := dask.persist(delayeds), desc=chunk_caption)
                del delayeds, result
                # cleanup - sometimes writing NetCDF handlers are not closed immediately and block reading access
                import gc; gc.collect()
                # cleanup - release all workers memory, call garbage collector before to prevent heartbeat errors
                if timeout is not None:
                    client.restart(timeout=timeout, wait_for_workers=True)
#                 # more granular control
#                 n_workers = len(client.nthreads())
#                 client.restart(wait_for_workers=False)
#                 client.wait_for_workers(n_workers, timeout=timeout)
            # update chunks counter
            counter += len(chunk)

    def delete_stack(self, name, basedir='auto'):
        import os

        if name == '' or name is None:
            name = ''
        else:
            name = name + '_'

        filenames = self._glob_re(name + '[0-9]{8}(_[0-9]{8})*.nc', basedir=basedir)
        #print ('filenames', filenames)
        for filename in filenames:
            if os.path.exists(filename):
                os.remove(filename)
