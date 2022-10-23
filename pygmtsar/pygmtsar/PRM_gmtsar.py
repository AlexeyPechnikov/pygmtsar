#!/usr/bin/env python3
# Alexey Pechnikov, Oct, 2022, https://github.com/mobigroup/gmtsar

class PRM_gmtsar:

    @staticmethod
    def gmtsar_sharedir():
        """
        Call GMTSAR tool to obtain GMTSAR shared directory location.
        Note: Docker builds on Apple Silicon for linux/amd64 platform produce weird issue
        when gmtsar_sharedir.csh output switches to stderr in joblib calls with n_jobs > 1
        """
        import os
        import subprocess

        argv = ['gmtsar_sharedir.csh']
        p = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout_data, stderr_data = p.communicate()

        stderr_data = stderr_data.decode('utf8').strip()
        if stderr_data is not None and len(stderr_data):
            print ('gmtsar_sharedir', stderr_data)
        data = stdout_data.decode('utf8').strip()
        if data == '':
            return stderr_data
        return data

    def calc_dop_orb(self, earth_radius=0, doppler_centroid=0, inplace=False, debug=False):
        """
        Usage: calc_dop_orb  file.PRM  added.PRM  earth_radius  [doppler_centroid]
            file.PRM     - input name of PRM file 
            new.PRM      - output additional parameters to add to the PRM file 
            earth_radius - input set earth radius, 0 calculates radius 
            [doppler_centroid] - no parameter calculates doppler 
            [doppler_centroid] - use value (e.g. 0.0) to force doppler 
        """
        import subprocess
        import os
        from pygmtsar import PRM
        
        cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        p = subprocess.Popen(['calc_dop_orb', '/dev/stdin', '/dev/stdout', str(earth_radius), str(doppler_centroid)],
                             stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             cwd=cwd, encoding='ascii')
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
        SAT_baseline

        Usage: (two modes)

        mode 1:

            SAT_baseline PRM_master PRM_master

            This is used to compute height information
            (writes out height information for appending to PRM file)

        mode 2:

            SAT_baseline PRM_master PRM_aligned 

            PRM_master     PRM file for reference image
            PRM_aligned        PRM file of secondary image
            Please make sure the orbit file data is in PRM 
            Program runs through repeat orbit to find nearest point 
            to the start, center and end on the reference orbit
            (writes out parameters for appending to PRM file)
        """
        import os
        import subprocess
        from pygmtsar import PRM

        if not isinstance(other, PRM):
            raise Exception('Argument "other" should be PRM class instance')

        pipe1 = os.pipe()
        os.write(pipe1[1], bytearray(self.to_str(), 'ascii'))
        os.close(pipe1[1])
        #print ('descriptor 1', str(pipe1[0]))

        pipe2 = os.pipe()
        os.write(pipe2[1], bytearray(other.to_str(), 'ascii'))
        os.close(pipe2[1])
        #print ('descriptor 2', str(pipe2[0]))

        argv = ['SAT_baseline', f'/dev/fd/{pipe1[0]}', f'/dev/fd/{pipe2[0]}']
        cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, pass_fds=[pipe1[0], pipe2[0]],
                             cwd=cwd, encoding='ascii')
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
         Usage: SAT_llt2rat master.PRM prec [-bo[s|d]] < inputfile > outputfile

             master.PRM   -  parameter file for master image and points to LED orbit file.
             precise      -  (0) standard back geocoding, (1) - polynomial refinenent (slower).
             inputfile    -  lon, lat, elevation [ASCII].
             outputfile   -  range, azimuth, elevation(ref to radius in PRM), lon, lat [ASCII default].
             -bos or -bod -  binary single or double precision output.

             Note: -bos mode support deleted as obsolete
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
        os.write(pipe[1], bytearray(self.to_str(), 'ascii'))
        os.close(pipe[1])
        #print ('descriptor', str(pipe[0]))

        argv = ['SAT_llt2rat', f'/dev/fd/{pipe[0]}', str(precise)]
        # set binary format mode
        if binary:
            argv.append('-bod')
        #print (argv)
        cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, pass_fds=[pipe[0]],
                             cwd=cwd, bufsize=10*1000*1000)
        stdout_data, stderr_data = p.communicate(input=stdin_data)

        stderr_data = stderr_data.decode('ascii')
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

    def resamp(self, alignedPRM, alignedSLC_tofile, interp, debug=False):
        """
        Usage: resamp master.PRM aligned.PRM new_aligned.PRM new_aligned.SLC interp
        master.PRM         - PRM for master imagea
        aligned.PRM        - PRM for aligned image
        new_aligned.PRM    - PRM for aligned aligned image
        new_aligned.SLC    - SLC for aligned aligned image
        interp             - interpolation method: 1-nearest; 2-bilinear; 3-biquadratic; 4-bisinc

        master  = PRM.from_file('master.PRM')
        aligned = PRM.from_file('aligned.PRM')
        master.resamp(aligned,'new_aligned.SLC',interp).to_file('new_aligned.PRM')
        """
        import os
        import subprocess
        from pygmtsar import PRM

        if not isinstance(alignedPRM, PRM):
            raise Exception('Argument should be PRM class instance')

        pipe1 = os.pipe()
        os.write(pipe1[1], bytearray(alignedPRM.to_str(), 'ascii'))
        os.close(pipe1[1])
        #print ('descriptor 1', pipe1[0])

        pipe2 = os.pipe()
        #print ('descriptor 2', pipe2[1])

        # Usage: resamp master.PRM aligned.PRM new_aligned.PRM new_aligned.SLC interp
        #
        #cmd = f'resamp /dev/stdin /dev/fd/{pipe1[0]} /dev/fd/{pipe2[1]} /dev/stdout {intrp} | sponge ___'
        #cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        #p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
        #                     stderr=subprocess.PIPE, pass_fds=[pipe1[0], pipe2[1]],
        #                     cwd=cwd, encoding='ascii', shell=True)
        argv = ['resamp', f'/dev/stdin', f'/dev/fd/{pipe1[0]}',
                f'/dev/fd/{pipe2[1]}', '/dev/stdout', str(interp)]
        cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, pass_fds=[pipe1[0], pipe2[1]],
                             cwd=cwd)
        stdout_data, stderr_data = p.communicate(input=bytearray(self.to_str(), 'ascii'))

        # print errors and notifications
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: resamp', stderr_data.decode('ascii'))

        # save big SLC binary file
        with open(alignedSLC_tofile, 'wb') as f:
            f.write(stdout_data)

        # PRM file should be small text
        data = os.read(pipe2[0],int(10e6)).decode('ascii')
        return PRM.from_str(data)

    # TODO: add topo_ra argument processing
    # two binary files real.grd and imag.grd will be created
    # TBD: update phasediff tool to allow output files basename argument
    def phasediff(self, other, topo_ra_fromfile, imag_tofile, real_tofile, debug=False):
        """
        phasediff [GMTSAR] - Compute phase difference of two images

        Usage: phasediff ref.PRM rep.PRM [-topo topo_ra.grd] [-model modelphase.grd] [-imag imag.grd] [-real real.grd]
            (topo_ra and model in GMT grd format)
        """
        import os
        import subprocess
        from pygmtsar import PRM

        if not isinstance(other, PRM):
            raise Exception('Argument "other" should be PRM class instance')

        pipe1 = os.pipe()
        os.write(pipe1[1], bytearray(self.to_str(), 'ascii'))
        os.close(pipe1[1])
        #print ('descriptor 1', str(pipe1[0]))

        pipe2 = os.pipe()
        os.write(pipe2[1], bytearray(other.to_str(), 'ascii'))
        os.close(pipe2[1])
        #print ('descriptor 2', str(pipe2[0]))

        cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        argv = ['phasediff', f'/dev/fd/{pipe1[0]}', f'/dev/fd/{pipe2[0]}',
                '-imag', os.path.relpath(imag_tofile, cwd),
                '-real', os.path.relpath(real_tofile, cwd)]
        if topo_ra_fromfile is not None:
            argv.append('-topo')
            argv.append(os.path.relpath(topo_ra_fromfile,cwd))
        #print ('argv', argv)
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, pass_fds=[pipe1[0], pipe2[0]],
                             cwd=cwd, encoding='ascii')
        stdout_data, stderr_data = p.communicate()
        #print ('stdout_data', stdout_data)
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: phasediff', stderr_data)
        return

    # Usage: make_gaussian_filter name_of_PRM_file RNG_DEC AZI_DEC WAVELENGTH(m)
    # make_gaussian_filter S1_20150121_ALL_F1.PRM 2 1 400
    def make_gaussian_filter(self, range_dec, azi_dec, wavelength=200, debug=False):
        """
        make_gaussian_filter

        Usage: make_gaussian_filter name_of_PRM_file RNG_DEC AZI_DEC WAVELENGTH(m) [output_filename]
        Example: make_gaussian_filter IMG-HH-ALPSRP211830620-H1.0__A.PRM 2 4 200
        Output: gauss_WAVELENGTH or output_filename
        """
        import os
        import subprocess
        import numpy as np

        pipe1 = os.pipe()
        #print ('descriptor 1', pipe1[1])

        cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        argv = ['make_gaussian_filter', '/dev/stdin',
                str(range_dec), str(azi_dec), str(wavelength), f'/dev/fd/{pipe1[1]}']
        #print ('argv', argv)
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, pass_fds=[pipe1[1]],
                             cwd=cwd, encoding='ascii')
        stdout_data, stderr_data = p.communicate(self.to_str())
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: make_gaussian_filter', stderr_data)
            print ('DEBUG: make_gaussian_filter', stdout_data)

        data = os.read(pipe1[0],int(10e6)).decode('ascii')
        return [np.fromstring(stdout_data, dtype=int, sep=' '), data]

    # the command line tool supports only GMT binary format input (=bf) while NetCDF output is possible too
    def conv(self, idec, jdec, output_file, filter_file=None, filter_string=None, input_file=None, debug=False):
        """
        conv [GMTSAR] - 2-D image convolution

        Usage: conv idec jdec filter_file input output
           idec           - row decimation factor
           jdec           - column decimation factor
           filter_file    - eg. filters/gauss17x5
           input          - name of file to be filtered (I*2 or R*4)
           output         - name of filtered output file (R*4 only)

           examples:
           conv 4 2 filters/gauss9x5 IMG-HH-ALPSRP109430660-H1.0__A.PRM test.grd
           (makes and filters amplitude file from an SLC-file)

           conv 4 2 filters/gauss5x5 infile.grd outfile.grd
           (filters a float file)
        """
        import os
        import subprocess

        if filter_file is None and filter_string is None:
            raise Exception('Should be defined argument filter_file or filter_string')

        if input_file is None and self.filename is None:
            raise Exception('PRM object should be created from file or defined argument input_file')

        cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        if input_file is None:
            input_filename = os.path.basename(self.filename)
        else:
            input_filename = os.path.relpath(input_file, cwd)

        # open descriptors
        fds = []
        if filter_string is not None:
            pipe2 = os.pipe()
            os.write(pipe2[1], bytearray(filter_string, 'ascii'))
            os.close(pipe2[1])
            input_filter_filename = f'/dev/fd/{pipe2[0]}'
            fds.append(pipe2[0])
            #print ('descriptor 2', str(pipe2[0]))
        else:
            input_filter_filename = os.path.relpath(filter_file, cwd)

        if output_file is None:
            output_filename = '/dev/stdout'
        else:
            output_filename = os.path.relpath(output_file, cwd)

        argv = ['conv', str(idec), str(jdec), input_filter_filename, input_filename, output_filename]
        #print ('argv', argv)
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, pass_fds=fds,
                             cwd=cwd)
        stdout_data, stderr_data = p.communicate()
        #print ('stdout_data', len(stdout_data))
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: conv', stderr_data.decode('ascii'))
        return

    # actually, -alpha command line argument is useless
    # if amp files are defined set alpha = -1 to denote coherence-based alpha
    def phasefilt(self, imag_fromfile, real_fromfile, amp1_fromfile, amp2_fromfile,
                  phasefilt_tofile, corrfilt_tofile, psize=32, debug=False):
        """
        phasefilt [GMTSAR] - Apply adaptive non-linear phase filter

        USAGE:
            phasefilt -imag imag.grd -real real.grd [-alpha alpha][-psize size][-amp1 amp1.grd -amp2 amp2.grd][-diff][-v]
             applies Goldstein adaptive filter to phase [output: filtphase.grd]
             or applies modified Goldstein adaptive filter to phase [output: filtphase.grd, corrfilt.grd]
            -imag [required] GMT format file of imaginary component
            -real [required] GMT format file of real component
            -alpha  exponent for filter - usually between 0.0 and 1.5 (0.0 should not filter).
                    default: 0.5 [Goldstein filter] (anything above 1.0 may be excessive)
                    alpha < 0 will set alpha = (1 - coherence) [modified Goldstein]
            -psize patch size for filtering. Must be power of two.
                    default: 32
            -amp1 GMT format file of amplitude image of image 1. Needed (and applies) modified filter.
            -amp2 GMT format file of amplitude image of image 2. Needed (and applies) modified filter.
        """
        import os
        import subprocess

        if corrfilt_tofile is None and alpha<0:
            raise Exception('For argument alpha<0 shoukd be defined argument corrfilt_tofile')

        cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        # -v - only args after verbose (if set) will be echoed
        argv = ['phasefilt',
                '-imag', os.path.relpath(imag_fromfile, cwd),
                '-real', os.path.relpath(real_fromfile, cwd),
                '-amp1', os.path.relpath(amp1_fromfile, cwd),
                '-amp2', os.path.relpath(amp2_fromfile, cwd),
                '-phasefilt', os.path.relpath(phasefilt_tofile, cwd),
                '-corrfilt', os.path.relpath(corrfilt_tofile, cwd),
                '-psize', str(psize)]
        #print ('argv', argv)

        p = subprocess.Popen(argv, stderr=subprocess.PIPE, cwd=cwd)
        stderr_data = p.communicate()[1]
        if len(stderr_data) > 0 and debug:
            print ('DEBUG: phasefilt', stderr_data.decode('ascii'))
        return

    def SAT_look(self, coords=None, fromfile=None, tofile=None, binary=False, debug=False):
        """
        Usage: SAT_look master.PRM [-bo[s|d]] < inputfile > outputfile

            master.PRM   -  parameter file for master image and points to LED orbit file
            inputfile    -  lon, lat, elevation [ASCII]
            outputfile   -  lon, lat, elevation look_E look_N look_U [ASCII default]
            -bos or -bod -  binary single or double precision output

                Note: -bos mode does not work

        example: SAT_look master.PRM < topo.llt > topo.lltn

        Note that the output elevation is the one above reference radius specified in the PRM file
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
        os.write(pipe[1], bytearray(self.to_str(), 'ascii'))
        os.close(pipe[1])
        #print ('descriptor', str(pipe[0]))

        argv = ['SAT_look', f'/dev/fd/{pipe[0]}']
        # set binary format mode
        if binary:
            argv.append('-bod')
        #print (argv)
        cwd = os.path.dirname(self.filename) if self.filename is not None else '.'
        p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, pass_fds=[pipe[0]],
                             cwd=cwd, bufsize=10*1000*1000)
        stdout_data, stderr_data = p.communicate(input=stdin_data)

        stderr_data = stderr_data.decode('ascii')
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
