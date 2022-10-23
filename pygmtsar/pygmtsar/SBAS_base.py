#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
from .tqdm_joblib import tqdm_joblib
from .datagrid import datagrid

class SBAS_base(tqdm_joblib, datagrid):

    def __repr__(self):
        return 'Object %s %d items\n%r' % (self.__class__.__name__, len(self.df), self.df)

    def to_dataframe(self):
        return self.df

    def dump(self, to_path=None):
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
        import pickle
        import os

        if os.path.isdir(from_path):
            sbas_pickle = os.path.join(from_path, 'SBAS.pickle')
        else:
            sbas_pickle = from_path

        print (f'NOTE: load state from file {sbas_pickle}')
        return pickle.load(open(sbas_pickle, 'rb'))


    def backup(self, backup_dir, copy=False, debug=False):
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
        Set master image date
        """
        if not master in self.df.index:
            raise Exception('Master image not found')
        self.master = master
        return self

    def get_master(self, subswath=None):
        """
        Return dataframe master record(s) for all or only selected subswath
        """
        df = self.df.loc[[self.master]]
        if not subswath is None:
            df = df[df.subswath == subswath]
        assert len(df) > 0, f'Master record for subswath {subswath} not found'
        return df

    def get_aligned(self, subswath=None, date=None):
        """
        Return dataframe aligned records (excluding master) for selected subswath
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
        import numpy as np
        # note: df.unique() returns unsorted values so it would be 21 instead of expected 12
        return np.unique(self.df.subswath)
    
    def get_subswath(self, subswath=None):
        """
        Check and return subswath or return an unique subswath to functions which work with a single subswath only.
        """
        # detect all the subswaths
        subswaths = self.get_subswaths()
        assert subswath is None or subswath in subswaths, f'ERROR: subswath {subswath} not found'
        if subswath is not None:
            return subswath
        assert len(subswaths)==1, f'ERROR: multiple subswaths {subswaths} found, merge them first using SBAS.merge_parallel()'
        # define subswath
        return subswaths[0]

    # use the function for open_grids() and save_grids()
    def get_filenames(self, subswath, pairs, name, add_subswath=True):
        import pandas as pd
        import numpy as np
        import os

        # define subswath when subswath=None and check otherwise
        subswath = self.get_subswath(subswath)
    
        if pairs is None:
            # special case to open a single grid {name}.grd or a set of subswath grids Fn_{name}.grd
            pass
        elif isinstance(pairs, pd.DataFrame):
            pairs = pairs.values
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

    def save_grids(self, grids, name, func=None, add_subswath=True, n_jobs=1, **kwargs):
        import numpy as np
        from tqdm.auto import tqdm
        import joblib
        import os

        subswaths = self.get_subswaths()
        # pipeline before merging subswaths is well defined and we do not need to use this function
        assert len(subswaths) == 1, 'ERROR: use SBAS.merge() to merge multiple subswaths before'
    
        assert len(grids.dims) in [2, 3], 'ERROR: supported 2D and 3D arrays only'
        if len(grids.dims) ==3 :
            grids = [da for da in grids]
    
        # special case for a single grid
        if not isinstance(grids, list):
            filename = self.get_filenames(None, None, name, add_subswath=add_subswath)
            grids.astype(np.float32).rename(name)\
                .to_netcdf(filename, encoding={name: self.compression}, engine=self.engine)
            return
    
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
                .to_netcdf(filename, encoding={name: self.compression}, engine=self.engine)
            return

        # test
        #for subswath in subswaths:
        #    for grid in grids:
        #        preprocess(subswath, grid)
        #return

        # process all the grids
        description = 'Saving' if func is None else 'Processing and Saving'
        with self.tqdm_joblib(tqdm(desc=description, total=len(grids)*len(subswaths))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(preprocess)(subswath, grid) \
                                                     for subswath in subswaths for grid in grids)

    # returns all grids in basedir by mask or grids by dates and name
    # Backward-compatible open_grids() returns list of grids fot the name or a single grid for a single subswath
    def open_grids(self, pairs, name, geocode=False, inverse_geocode=False,  mask=None, func=None,
                   crop_valid=False, add_subswath=True, chunks=None, n_jobs=-1):
        import pandas as pd
        import xarray as xr
        import numpy as np
        from tqdm.auto import tqdm
        import joblib

        assert not(geocode and inverse_geocode), 'ERROR: Only single geocoding option can be applied'

        if chunks is None:
            chunks = self.chunksize

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
            pairs = pairs.values
        else:
            pairs = np.asarray(pairs)

        def postprocess(da, subswath):
            if self.is_ra(da) and geocode:
                da = self.intf_ra2ll(subswath, da)
            elif self.is_geo(da) and inverse_geocode:
                da = self.intf_ll2ra(subswath, da)
            if func is not None:
                if isinstance(func, list):
                    for f in func:
                        da = f(da)
                else:
                    da = func(da)
            if mask is not None:
                assert self.is_same(mask, da), 'ERROR: mask defined in different coordinates'
                da = nanmask * da
            return da

        dass = []
        for subswath in subswaths:
            filenames = self.get_filenames(subswath, pairs=pairs, name=name, add_subswath=add_subswath)
        
            das = []
            if pairs is None:
                # special case to open a single grid {name}.grd or a set of subswath grids Fn_{name}.grd
                #print ('filename', filename)
                da = xr.open_dataarray(filenames, engine=self.engine, chunks=chunks)
                das  = postprocess(da, subswath)
            elif len(pairs.shape) == 1:
                # read all the grids from files
                for filename in filenames:
                    #print (date, filename)
                    da = xr.open_dataarray(filename, engine=self.engine, chunks=chunks)
                    das.append(da)

                # post-processing on a set of 2D rasters
                with self.tqdm_joblib(tqdm(desc='Loading', total=len(das))) as progress_bar:
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
                    da = xr.open_dataarray(filename, engine=self.engine, chunks=chunks)
                    das.append(da)

                # post-processing on a set of 2D rasters
                with self.tqdm_joblib(tqdm(desc='Loading', total=len(das))) as progress_bar:
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
