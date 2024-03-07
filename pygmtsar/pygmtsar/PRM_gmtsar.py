# ----------------------------------------------------------------------------
# PyGMTSAR
# 
# This file is part of the PyGMTSAR project: https://github.com/mobigroup/gmtsar
# 
# Copyright (c) 2022, Alexey Pechnikov
# 
# Licensed under the BSD 3-Clause License (see LICENSE for details)
# ----------------------------------------------------------------------------

class PRM_gmtsar:

    @staticmethod
    def gmtsar_sharedir():
        """
        Get the location of the GMTSAR shared directory.

        Returns
        -------
        str
            The path to the GMTSAR shared directory.

        Notes
        -----
        This method calls the GMTSAR tool 'gmtsar_sharedir.csh' to obtain the location of the GMTSAR shared directory.
        The shared directory contains essential files and data used by GMTSAR.
        """
        import os
        import subprocess

        argv = ['gmtsar_sharedir.csh']
        p = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8')
        stdout_data, stderr_data = p.communicate()

        stderr_data = stderr_data.strip()
        if stderr_data is not None and len(stderr_data):
            print ('gmtsar_sharedir', stderr_data)
        data = stdout_data.strip()
        if data == '':
            return stderr_data
        return data

    def calc_dop_orb(self, earth_radius=0, doppler_centroid=0, inplace=False, debug=False):
        """
        Calculate the Doppler orbit.

        Parameters
        ----------
        earth_radius : float, optional
            The Earth radius. If set to 0, the radius will be calculated. Default is 0.
        doppler_centroid : float, optional
            The Doppler centroid. If set to 0, the Doppler will be calculated. Default is 0.
        inplace : bool, optional
            If True, the calculated Doppler orbit will be set in the current PRM object. If False, a new PRM object
            with the calculated Doppler orbit will be returned. Default is False.
        debug : bool, optional
            If True, debug information will be printed. Default is False.

        Returns
        -------
        PRM or None
            If inplace is True, returns the current PRM object with the calculated Doppler orbit.
            If inplace is False, returns a new PRM object with the calculated Doppler orbit.
        """
        import subprocess
        import os
        from pygmtsar import PRM
        
        cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        p = subprocess.Popen(['calc_dop_orb', '/dev/stdin', '/dev/stdout', str(earth_radius), str(doppler_centroid)],
                             stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             cwd=cwd, encoding='utf8')
        stdout_data, stderr_data = p.communicate(input=self.to_str())
        #print ('stdout_data', stdout_data)
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: calc_dop_orb', stderr_data)
        prm = PRM.from_str(stdout_data)
        if inplace:
            return self.set(prm)
        else:
            return prm

    def SAT_baseline(self, other, tail=None, debug=False):
        """
        Compute the satellite baseline.

        Parameters
        ----------
        other : PRM
            The PRM object for the other image.
        tail : int, optional
            The number of lines to keep from the computed baseline. Default is None, which keeps all lines.
        debug : bool, optional
            If True, debug information will be printed. Default is False.

        Returns
        -------
        PRM
            The PRM object with the computed satellite baseline.

        Notes
        -----
        This method computes the satellite baseline between the current PRM object and the specified PRM object.
        The resulting parameters are stored in the returned PRM object.
        """
        import os
        import subprocess
        from pygmtsar import PRM

        if not isinstance(other, PRM):
            raise Exception('Argument "other" should be PRM class instance')

        pipe1 = os.pipe()
        os.write(pipe1[1], bytearray(self.to_str(), 'utf8'))
        os.close(pipe1[1])
        #print ('descriptor 1', str(pipe1[0]))

        pipe2 = os.pipe()
        os.write(pipe2[1], bytearray(other.to_str(), 'utf8'))
        os.close(pipe2[1])
        #print ('descriptor 2', str(pipe2[0]))

        argv = ['SAT_baseline', f'/dev/fd/{pipe1[0]}', f'/dev/fd/{pipe2[0]}']
        if debug:
            print ('DEBUG: argv', argv)
        cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, pass_fds=[pipe1[0], pipe2[0]],
                             cwd=cwd, encoding='utf8')
        stdout_data, stderr_data = p.communicate()
        #print ('stdout_data', stdout_data)
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: SAT_baseline', stderr_data)
        prm = PRM.from_str(stdout_data)
        # replacement for SAT_baseline $1 $2 | tail -n9
        if tail is not None:
            prm.df = prm.df.tail(tail)
        return prm

    """
    coords = prm.SAT_llt2rat([-115.588333, 32.758333, -42.441303], precise=1)
    coords = prm.SAT_llt2rat([-115.588333, 32.758333, -42.441303], precise=1, binary=True)

    coords = prm.SAT_llt2rat([-115.588333, 32.758333, -42.441303], tofile='out.dat', precise=1)
    coords = prm.SAT_llt2rat([-115.588333, 32.758333, -42.441303], tofile='outd.dat', precise=1, binary=True)

    coords = prm.SAT_llt2rat(dem_data[:10], precise=1)
    [format(v, '.6f') for v in coords]
    """
    def SAT_llt2rat(self, coords=None, fromfile=None, tofile=None, precise=1, binary=False, debug=False):
        """
        Convert latitude, longitude, and elevation (LLT) coordinates to radar (RAT) coordinates.

        Parameters
        ----------
        coords : array_like or None, optional
            The LLT coordinates to convert. Should be a 2D array-like object with shape (N, 3), where N is the number of coordinates.
            Each coordinate should be in the format [longitude, latitude, elevation].
        fromfile : str or None, optional
            The file path to read the LLT coordinates from. The file should contain a space-separated list of LLT coordinates in the
            format "longitude latitude elevation". If provided, the 'coords' parameter will be ignored. Default is None.
        tofile : str or None, optional
            The file path to save the converted RAT coordinates to. If not provided, the converted RAT coordinates will be returned as a numpy array.
        precise : int, optional
            The precision level of the conversion. Set to 0 for standard back geocoding or 1 for polynomial refinement (slower). Default is 1.
        binary : bool, optional
            If True, the output coordinates will be saved in binary format. Default is False.
        debug : bool, optional
            If True, debug information will be printed. Default is False.

        Returns
        -------
        numpy.ndarray or None
            If 'tofile' is None, returns a numpy array of the converted RAT coordinates with shape (N, 5), where N is the number of coordinates.
            Each coordinate is in the format [longitude, latitude, elevation, range, azimuth]. If 'tofile' is provided, returns None.

        Notes
        -----
        This method converts LLT coordinates to RAT coordinates using the current PRM object. The converted RAT coordinates can be saved
        to a file or returned as a numpy array.
        """
        import numpy as np
        from io import StringIO, BytesIO
        import os
        import subprocess

        if coords is not None and fromfile is None:
            #lon, lat, elevation = coords
            #data=f'{lon} {lat} {elevation}'
            buffer = BytesIO()
            # to produce the same result as original command returns
            # gmt grd2xyz --FORMAT_FLOAT_OUT=%lf topo/dem.grd -s
            np.savetxt(buffer, coords, delimiter=' ', fmt='%.6f')
            stdin_data = buffer.getvalue()
        elif coords is None and fromfile is not None:
            with open(fromfile, 'rb') as f:
                stdin_data = f.read()
        else:
            raise Exception('Should be defined data source as coordinates triplet (coords) or as file (fromfile)')

        # TBD: use np.array_split() and [x for x in a if x.size > 0] for chunking processing
        pipe = os.pipe()
        # for unicode file paths
        os.write(pipe[1], bytearray(self.to_str(), 'utf8'))
        os.close(pipe[1])
        #print ('descriptor', str(pipe[0]))

        argv = ['SAT_llt2rat', f'/dev/fd/{pipe[0]}', str(precise)]
        # set binary format mode
        if binary:
            argv.append('-bod')
        if debug:
            print ('DEBUG: argv', argv)
        if debug and not binary and coords is not None:
            print ('DEBUG: coords', coords)    
        cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, pass_fds=[pipe[0]],
                             cwd=cwd, bufsize=10*1000*1000)
        stdout_data, stderr_data = p.communicate(input=stdin_data)

        stderr_data = stderr_data.decode('utf8', errors='ignore')
        if stderr_data.startswith('interpolation point outside of data constraints'):
            print ('Error: SAT_llt2rat processing stopped due to invalid coordinates for one of input points')
            return None
        if stderr_data is not None and len(stderr_data) and debug:
            print ('DEBUG: SAT_llt2rat', stderr_data)

        if tofile is not None:
            with open(tofile, 'wb') as f:
                f.write(stdout_data)
        else:
            if binary:
                out = (np.frombuffer(stdout_data, dtype=np.dtype(np.float64)))
            else:
                out = np.fromstring(stdout_data, dtype=float, sep=' ')
            return out if out.size==5 else out.reshape(-1,5)

    def resamp(self, repeatPRM, repeatSLC_tofile, interp, debug=False):
        """
        Resample the repeat image.

        Parameters
        ----------
        repeatPRM : PRM
            The PRM object for the repeat image.
        repeatSLC_tofile : str
            The file path to save the resampled repeat SLC image to.
        interp : int
            The interpolation method: 1 for nearest, 2 for bilinear, 3 for biquadratic, or 4 for bisinc.
        debug : bool, optional
            If True, debug information will be printed. Default is False.

        Returns
        -------
        PRM
            The PRM object for the resampled repeat image.

        Notes
        -----
        This method resamples the repeat image using the current PRM object and the PRM object for the repeat image.
        The resampled repeat SLC image is saved to the specified file, and the resulting PRM parameters are returned.
        """
        import os
        import subprocess
        from pygmtsar import PRM

        if not isinstance(repeatPRM, PRM):
            raise Exception('Argument should be PRM class instance')

        pipe1 = os.pipe()
        os.write(pipe1[1], bytearray(repeatPRM.to_str(), 'utf8'))
        os.close(pipe1[1])
        #print ('descriptor 1', pipe1[0])

        pipe2 = os.pipe()
        #print ('descriptor 2', pipe2[1])

        # Usage: resamp master.PRM repeat.PRM new_repeat.PRM new_repeat.SLC interp
        #
        #cmd = f'resamp /dev/stdin /dev/fd/{pipe1[0]} /dev/fd/{pipe2[1]} /dev/stdout {intrp} | sponge ___'
        #cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        #p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
        #                     stderr=subprocess.PIPE, pass_fds=[pipe1[0], pipe2[1]],
        #                     cwd=cwd, encoding='utf8', shell=True)
        argv = ['resamp', f'/dev/stdin', f'/dev/fd/{pipe1[0]}',
                f'/dev/fd/{pipe2[1]}', '/dev/stdout', str(interp)]
        if debug:
            print ('DEBUG: argv', argv)
        cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, pass_fds=[pipe1[0], pipe2[1]],
                             cwd=cwd)
        stdout_data, stderr_data = p.communicate(input=bytearray(self.to_str(), 'utf8'))

        # print errors and notifications
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: resamp', stderr_data.decode('utf8'))

        # save big SLC binary file
        with open(repeatSLC_tofile, 'wb') as f:
            f.write(stdout_data)

        # PRM file should be small text
        data = os.read(pipe2[0],int(10e6))
        return PRM.from_str(data.decode('utf8'))

    def SAT_look(self, coords=None, fromfile=None, tofile=None, binary=False, debug=False):
        """
        Compute the satellite look vector.

        Parameters
        ----------
        coords : array_like or None, optional
            The LLT coordinates to compute the look vector for. Should be a 2D array-like object with shape (N, 3),
            where N is the number of coordinates. Each coordinate should be in the format [longitude, latitude, elevation].
        fromfile : str or None, optional
            The file path to read the LLT coordinates from. The file should contain a space-separated list of LLT coordinates
            in the format "longitude latitude elevation". If provided, the 'coords' parameter will be ignored. Default is None.
        tofile : str or None, optional
            The file path to save the computed look vectors to. If not provided, the computed look vectors will be returned as a numpy array.
        binary : bool, optional
            If True, the output look vectors will be saved in binary format. Default is False.
        debug : bool, optional
            If True, debug information will be printed. Default is False.

        Returns
        -------
        numpy.ndarray or None
            If 'tofile' is None, returns a numpy array of the computed look vectors with shape (N, 6), where N is the number of coordinates.
            Each look vector is in the format [longitude, latitude, elevation, look_E, look_N, look_U].
            If 'tofile' is provided, returns None.

        Notes
        -----
        This method computes the look vectors for LLT coordinates using the current PRM object. The computed look vectors can be saved
        to a file or returned as a numpy array.
        """
        import numpy as np
        from io import StringIO, BytesIO
        import os
        import subprocess

        if coords is not None and fromfile is None:
            #lon, lat, elevation = coords
            #data=f'{lon} {lat} {elevation}'
            buffer = BytesIO()
            np.savetxt(buffer, coords, delimiter=' ', fmt='%.6f')
            stdin_data = buffer.getvalue()
        elif coords is None and fromfile is not None:
            with open(fromfile, 'rb') as f:
                stdin_data = f.read()
        else:
            raise Exception('Should be defined data source as coordinates triplet (coords) or as file (fromfile)')

        pipe = os.pipe()
        os.write(pipe[1], bytearray(self.to_str(), 'utf8'))
        os.close(pipe[1])
        #print ('descriptor', str(pipe[0]))

        argv = ['SAT_look', f'/dev/fd/{pipe[0]}']
        # set binary format mode
        if binary:
            argv.append('-bod')
        if debug:
            print ('DEBUG: argv', argv)
        cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, pass_fds=[pipe[0]],
                             cwd=cwd, bufsize=10*1000*1000)
        stdout_data, stderr_data = p.communicate(input=stdin_data)

        stderr_data = stderr_data.decode('utf8')
        if stderr_data is not None and len(stderr_data) and debug:
            print ('DEBUG: SAT_look', stderr_data)
            return None

        if tofile is not None:
            with open(tofile, 'wb') as f:
                f.write(stdout_data)
        else:
            if binary:
                out = (np.frombuffer(stdout_data, dtype=np.dtype(np.float64)))
            else:
                out = np.fromstring(stdout_data, dtype=float, sep=' ')
            return out if out.size==6 else out.reshape(-1,6)
