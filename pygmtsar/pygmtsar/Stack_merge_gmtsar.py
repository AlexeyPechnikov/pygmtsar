# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2021, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------
from .Stack_intf import Stack_intf

class Stack_merge_gmtsar(Stack_intf):

    # stem_tofile + '.PRM' generating
    def merge_swath(self, conf, grid_tofile, stem_tofile, debug=False):
        """
        Merge the swaths.

        Parameters
        ----------
        conf : str
            The configuration file content.
        grid_tofile : str
            The file path for the output NetCDF grid.
        stem_tofile : str
            The file path for the output stem.

        Returns
        -------
        stdout_data : str
            The standard output from the merge_swath process.
        stderr_data : str
            The standard error from the merge_swath process.

        Notes
        -----
        This method merges the swaths using the merge_swath command-line tool. It takes the configuration file content as input and generates the output grid and stem files.
        """
        import subprocess

        argv = ['merge_swath', '/dev/stdin', grid_tofile, stem_tofile]
        if debug:
            print ('DEBUG: argv', argv)
            print ('DEBUG: conf:', f'\n{conf}')
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             encoding='ascii')
        stdout_data, stderr_data = p.communicate(input=conf)
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: merge_swath', stderr_data)
        if len(stdout_data) > 0 and debug:
            print ('DEBUG: merge_swath', stdout_data)

        return