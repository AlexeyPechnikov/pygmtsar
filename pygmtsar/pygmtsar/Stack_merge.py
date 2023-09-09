# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2021, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_merge_gmtsar import Stack_merge_gmtsar
from .PRM import PRM

class Stack_merge(Stack_merge_gmtsar):

    def merge(self, grid, debug=False):
        """
        Merge the radar grids.

        Parameters
        ----------
        grid : str
            The type of grid to merge, such as 'adi'.
        debug : bool, optional
            If True, debug information will be printed. Default is False.

        Notes
        -----
        This method merges the radar coordinate grids for a given grid name. It generates the necessary
        configuration files, PRM files, and performs the merge using the merge_swath command-line tool.

        """
        import numpy as np
        import os

        fullname = lambda filename: os.path.join(self.basedir, filename)
        def format_string(row):
            dt = row.name.replace("-","")
            subswath = row["subswath"]
            return fullname(f'S1_{dt}_ALL_F{subswath}.PRM') + ':' + fullname(f'F{subswath}_{grid}.grd')

        subswaths = self.get_subswaths()    
        if len(subswaths) == 1:
            # only one subswath found, merging is not possible
            return

        filenames = [fullname(f'F{subswath}_{grid}.grd') for subswath in subswaths]
        if not np.all([os.path.exists(filename) for filename in filenames]):
            print (f'Subswath grids missed, skip merging: {grid}')
            return
        #print ('filenames', filenames)

        subswaths_str = int(''.join(map(str, subswaths)))
        grid_tofile = fullname(f'F{subswaths_str}_adi.grd')
        tmp_stem_tofile  = fullname(f'F{subswaths_str}_adi')

        config = '\n'.join(self.get_reference().apply(format_string, axis=1).tolist())
        #print ('merge_swath', config, grid_tofile, tmp_stem_tofile)
        self.merge_swath(config, grid_tofile, tmp_stem_tofile, debug=debug)
    
        # cleanup - files should exists as these are processed above
        for filename in filenames:
            if debug:
                print ('DEBUG: remove', filename)
            os.remove(filename)

    def merge_intf(self, pair, grid, debug=False):
        """
        Merge the interferograms.

        Parameters
        ----------
        pair : tuple
            A tuple containing the reference and repeat dates of the interferogram pair.
        grid : str
            The type of grid to merge, such as 'phasefilt' or 'corr'.
        debug : bool, optional
            If True, debug information will be printed. Default is False.

        Notes
        -----
        This method merges the interferograms for a given pair and grid name. It generates the necessary
        configuration files, PRM files, and performs the merge using the merge_swath command-line tool.

        """
        import os

        record2multistem = lambda record: self.multistem_stem(record.subswath, record.datetime)[0]
        fullname = lambda filename: os.path.join(self.basedir, filename)

        # extract dates from pair
        date1, date2 = pair
        #print (date1, date2)

        # records should be sorted by datetime that's equal to sorting by date and subswath
        multistems1 = self.get_repeat(None, date1).apply(record2multistem, axis=1)
        multistems2 = self.get_repeat(None, date2).apply(record2multistem, axis=1)
        if len(multistems1) == 1:
            # only one subswath found, merging is not possible
            return

        config = []
        subswaths = []
        cleanup = []
        for (multistem1, multistem2) in zip(multistems1, multistems2):
            # F2_20220702_20220714_phasefilt.grd
            prm1_filename = fullname(multistem1 + '.PRM')
            prm2_filename = fullname(multistem1 + '.PRM')
            prm_filename  = fullname(multistem1 + f'_{grid}.PRM')

            prm1 = PRM.from_file(prm1_filename)
            prm2 = PRM.from_file(prm2_filename)
            rshift = prm2.get('rshift')
            #print ('rshift', rshift)
            #assert rshift == 0, 'rshift is not equal to zero for reference PRM'
            fs1 = prm1.get('first_sample')
            fs2 = prm2.get('first_sample')
            #print ('fs1, fs2', fs1, fs2)
            #assert fs1 == fs2, 'first_sample is not equal for reference and repeat PRM'
            prm = prm1.set(rshift=rshift, first_sample=fs2 if fs2 > fs1 else fs1).to_file(prm_filename)

            subswath = int(multistem1[-1:])
            subswaths.append(subswath)
            dt1 = multistem1[3:11]
            dt2 = multistem2[3:11]
            #print (multistem1, multistem2, fullname(multistem1))
            grid_fromfile = fullname(f'F{subswath}_{dt1}_{dt2}_{grid}.grd')
            cleanup.append(grid_fromfile)
            #print (prm_filename, grid_filename)
            config.append(':'.join([prm_filename, grid_fromfile]))
        config = '\n'.join(config)
        subswaths = int(''.join(map(str,subswaths)))
        # F23_20220702_20220714_phasefilt.grd
        grid_tofile = fullname(f'F{subswaths}_{dt1}_{dt2}_{grid}.grd')
        # F23_20220702_20220714_phasefilt
        tmp_stem_tofile = fullname(f'F{subswaths}_{dt1}_{dt2}_{grid}')
        #print ('grid_tofile', grid_tofile)
        #print (config)
        # S1_20220702_ALL_F23 without extension
        stem_tofile = fullname(f'S1_{dt1}_ALL_F{subswaths}')

        # use temporary well-qualified stem file name to prevent parallel processing conflicts
        self.merge_swath(config, grid_tofile, tmp_stem_tofile, debug=debug)
        # different pairs and grids generate the same PRM file, replace it silently
        if debug:
            print ('DEBUG: replace', f'{tmp_stem_tofile}.PRM', f'{stem_tofile}.PRM')
        os.replace(f'{tmp_stem_tofile}.PRM', f'{stem_tofile}.PRM')

        # cleanup - files should exists as these are processed above
        for filename in cleanup:
            if debug:
                print ('DEBUG: remove', filename)
            os.remove(filename)

    def merge_parallel(self, pairs, intfs = ['phasefilt', 'corr'], grids=['adi'], n_jobs=1, **kwargs):
        """
        Perform parallel merging of interferograms.

        Parameters
        ----------
        pairs : list or None
            List of interferogram date pairs to merge. If None, all available pairs will be merged.
        intfs : list, optional
            List of grid names to merge, such as 'phasefilt' or 'corr'. Default is ['phasefilt', 'corr'].
        n_jobs : int, optional
            The number of parallel jobs to run. Default is -1, which uses all available CPU cores.
        **kwargs : additional keyword arguments
            Additional arguments to be passed to the merge function.

        Examples
        --------
        stack.merge_parallel(pairs)

        Notes
        -----
        This method performs parallel merging of the interferograms for the specified pairs and grids. It utilizes
        joblib.Parallel for efficient parallel processing. The merged interferograms are stored in the 'df' attribute
        of the Stack object, which is a GeoDataFrame containing information about the original or merged interferograms.
        """
        from tqdm.auto import tqdm
        import joblib
        import pandas as pd
        import geopandas as gpd

        # merging is not applicable to a single subswath
        # for this case coordinate transformation matrices already built in Stack.intf_parallel()
        subswaths = self.get_subswaths()
        if len(subswaths) == 1:
            return
        
        if n_jobs != 1:
            print ('Note: merging uses GMTSAR tool which is not robust in multi-threading mode and n_jobs should be always equal to 1')
            n_jobs = 1

        # convert pairs (list, array, dataframe) to 2D numpy array
        pairs = self.get_pairs(pairs)[['ref', 'rep']].astype(str).values

        with self.tqdm_joblib(tqdm(desc=f'Merging Interferograms', total=len(pairs)*len(intfs))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(self.merge_intf)(pair, grid, **kwargs) \
                                           for pair in pairs for grid in intfs)

        with self.tqdm_joblib(tqdm(desc=f'Merging ADI', total=len(pairs)*len(grids))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(self.merge)(grid, **kwargs) for grid in grids)

        df = self.df.groupby(self.df.index).agg({'datetime': min, 'orbit': min, 'mission': min, 'polarization':min,
                                            'subswath': lambda s: int(''.join(map(str,list(s)))),
                                            'datapath': lambda p: list(p),
                                            'metapath': lambda p: list(p),
                                            'orbitpath': min,
                                            'geometry': lambda g: g.unary_union
                                           })
        self.df = gpd.GeoDataFrame(df)
    
        # build topo_ra and geo transform matrices for the merged interferograms
        #self.topo_ra_parallel()
        #self.transforms(pairs)
