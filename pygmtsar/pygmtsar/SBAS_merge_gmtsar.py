#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
from .SBAS_intf import SBAS_intf

class SBAS_merge_gmtsar(SBAS_intf):

    # stem_tofile + '.PRM' generating
    def merge_swath(self, conf, grid_tofile, stem_tofile, debug=False):
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