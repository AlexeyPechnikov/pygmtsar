#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
from .SBAS_merge_gmtsar import SBAS_merge_gmtsar
from .PRM import PRM

class SBAS_merge(SBAS_merge_gmtsar):

    def merge(self, pair, grid, debug=False):
        import os

        record2multistem = lambda record: self.multistem_stem(record.subswath, record.datetime)[0]
        fullname = lambda filename: os.path.join(self.basedir, filename)

        # extract dates from pair
        date1, date2 = pair
        #print (date1, date2)

        # records should be sorted by datetime that's equal to sorting by date and subswath
        multistems1 = self.get_aligned(None, date1).apply(record2multistem, axis=1)
        multistems2 = self.get_aligned(None, date2).apply(record2multistem, axis=1)
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
            #assert rshift == 0, 'rshift is not equal to zero for master PRM'
            fs1 = prm1.get('first_sample')
            fs2 = prm2.get('first_sample')
            #print ('fs1, fs2', fs1, fs2)
            #assert fs1 == fs2, 'first_sample is not equal for master and repeat PRM'
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

    def merge_parallel(self, pairs, grids = ['phasefilt', 'corr'], n_jobs=-1, **kwargs):
        from tqdm.auto import tqdm
        import joblib
        import pandas as pd
        import geopandas as gpd

        # merging is not applicable to a single subswath
        # for this case coordinate transformation matrices already built in SBAS.intf_parallel()
        subswaths = self.get_subswaths()
        if len(subswaths) == 1:
            return
        
        if isinstance(pairs, pd.DataFrame):
            pairs = pairs.values

        with self.tqdm_joblib(tqdm(desc=f'Merging Subswaths', total=len(pairs)*len(grids))) as progress_bar:
            joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(self.merge)(pair, grid, **kwargs) \
                                           for pair in pairs for grid in grids)

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
