#!/usr/bin/env python3
# Alexey Pechnikov, Sep, 2021, https://github.com/mobigroup/gmtsar
# python3 -m pip install install pandas --upgrade
# Wrapper to read, write and update PRM files and calculate Doppler orbit by calc_dop_orb command line tool
#import pytest

class PRM:
    @staticmethod
    def from_list(prm_list):
        from io import StringIO
        prm = StringIO('\n'.join(prm_list))
        return PRM._from_io(prm)

    @staticmethod
    def from_str(prm_string):
        from io import StringIO
        prm = StringIO(prm_string)
        return PRM._from_io(prm)

    @staticmethod
    def from_file(prm_filename):
        #data = json.loads(document)
        prm = PRM._from_io(prm_filename)
        prm.filename = prm_filename
        return prm

    @staticmethod
    def _from_io(prm):
        import pandas as pd
        return PRM(pd.read_csv(prm, sep='\s+=\s+', header=None, names=['name', 'value'], engine='python').set_index('name')\
                    .applymap(lambda val : pd.to_numeric(val,errors='ignore')))

    def __init__(self, prm):
        #print ('__init__')
        self.PRM = prm.reset_index().drop_duplicates(keep='last', inplace=False).set_index('name')
        self.filename = None

    def __eq__(self, other):
        return isinstance(self, PRM) and self.PRM == other.PRM

    def __str__(self):
        return self.to_str()

    def __repr__(self):
        if self.filename:
            return 'Object %s (%s) %d items\n%r' % (self.__class__.__name__, self.filename, len(self.PRM), self.PRM)
        else:
            return 'Object %s %d items\n%r' % (self.__class__.__name__, len(self.PRM), self.PRM)

    def set(self, prm=None, **kwargs):
        if isinstance(prm, PRM):
            for (key, value) in prm.PRM.itertuples():
                self.PRM.loc[key] = value
        elif prm is not None:
            raise Exception('Arguments is not a PRM object')
        for key, value in kwargs.items():
            self.PRM.loc[key] = value
        return self

    def to_dataframe(self):
        return self.PRM

    def to_file(self, prm):
        return self._to_io(prm)

    def update(self):
        if self.filename is None:
            raise Exception('PRM is not created from file, use to_file() method instead')
        return self._to_io(self.filename)

    def to_str(self):
        return self._to_io()

    def _to_io(self, prm=None):
        return self.PRM.reset_index().astype(str).apply(lambda row: (' = ').join(row), axis=1)\
            .to_csv(prm, header=None, index=None)

    def sel(self, *args):
        return PRM(self.PRM.loc[[*args]])

    def __add__(self, other):
        import pandas as pd
        if isinstance(other, PRM):
            prm = pd.concat([self.PRM, other.PRM])
            # drop duplicates
            prm = prm.groupby(prm.index).last()
        else:
            prm = self.PRM + other
        return PRM(prm)

    def get(self, *args):
        out = [self.PRM.loc[[key]].iloc[0].values[0] for key in args]
        if len(out) == 1:
            return out[0]
        return out

    def calc_dop_orb(self, earth_radius='0', doppler_centroid='0'):
        import subprocess
        import os
        cwd = os.path.dirname(self.filename) if self.filename is not None else ''
        p = subprocess.Popen(['calc_dop_orb', '/dev/stdin', '/dev/stdout', str(earth_radius), str(doppler_centroid)],
                     stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                     cwd=cwd, encoding='utf8')
        stdout_data = p.communicate(input=self.to_str())[0]
        return PRM.from_str(stdout_data)

def main():
    import sys

    if len(sys.argv) <= 2 :
        print (f"Usage: {sys.argv[0]} file.PRM parameter1 value1 parameter2 value2 ...")
        exit(0)

    prm_filename = sys.argv[1]
    #print ('prm_filename', prm_filename)
    prm = PRM.from_file(prm_filename)
    pairs = dict(zip(sys.argv[2::2], sys.argv[3::2]))
    prm = prm.set(**pairs)
    prm.to_file(prm_filename)

if __name__ == "__main__":
    # execute only if run as a script
    main()

"""
PRM.from_file('S1_20150403_ALL_F1.PRM').set(nrows=-999).update()

prm1 = PRM.from_file('S1_20150403_ALL_F1.PRM')
prm2 = PRM.from_file('S1_20150403_ALL_F1.PRM').sel('nrows') + 1000
(prm1).get('nrows'), (prm2).get('nrows'), (prm1 + prm2).get('nrows'), (prm1).get('nrows')

prm1.set(prm2).update()

(prm1 + prm2).to_file('S1_20150403_ALL_F1.PRM')

(prm1).get('nrows')

prm = PRM.from_file(...)
prm.set(prm.calc_dop_orb(0,0)).update()
"""
