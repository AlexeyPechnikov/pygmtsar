# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2021, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .SBAS_base import SBAS_base

class SBAS_orbits(SBAS_base):
    # for precision orbit there is only single orbit per day
    # for approximate orbit 2 and maybe more orbits per day are possible
    # so check orbit file for for each subswath
    def download_orbits(self):
        """
        Download missed orbits for all the SBAS scenes.

        Returns
        -------
        None

        Examples
        --------
        sbas.download_orbits()
        """
        from eof.download import download_eofs

	# For each date get a single product with missing orbit
        sbas_df_missingorb = sbas.df[sbas.df.orbitpath.isna()]
        sbas_df_missingorb = sbas_df_missingorb[~sbas_df_missingorb.index.duplicated(keep="first")]

        # Apply orbit to products of same date
        for record in sbas_df_missingorb.itertuples():
            orbitpath = download_eofs([record.datetime], [record.mission], save_dir=self.basedir)[0]
            self.df.at[record.Index, "orbitpath"] = orbitpath
