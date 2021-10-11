//
//  split_spectrum.c
//
//  Created by Xiaohua Xu on 7/18/18.
//
//  Used to estimate ionospheric delay in interferograms.
//

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * DATE                                                                    *
 ***************************************************************************/

#include "gmt.h"
#include "gmtsar.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "tiffio.h"
//#include <complex.h>

// void right_shift(fcomplex *, int);
// void left_shift(fcomplex *, int);
void circ_shift(fcomplex *, fcomplex *, int, int);
// void bandpass_filter(double, double, int, int, int, double, fcomplex *);
// void kaiser(double, int, double *);
// double bessi0(double);
int fliplr(double *, int);
int cos_window(double, double, double, int, double *);
int split1(int, char **);
int split2(int, char **);
fcomplex fmean(short *, int);

char *USAGE1 = "\nUSAGE: split_spectrum prm1 [split_half] \n\n"
               "   program used to split range spectrum for SLC using a modified cosine "
               "filter\n\n"
               "   SLCs are bandpassed and then shifted to the center of the spectrum\n\n"
               "   split_half is build for ALOS FBD FBS cases, put 1 for using half the spectrum\n\n"
               "   outputs are SLCH SLCL, for TOPS data, outputs are high.tiff low.tiff\n\n";


int main(int argc, char **argv) {
	if (argc != 3 && argc != 2)
		die("", USAGE1);
	split1(argc, argv);
	return (1);
}

int split1(int argc1, char **argv1) {

	FILE *SLC_file1, *SLCH, *SLCL;
    TIFF *tif, *tifh, *tifl;
    uint32 width, height, nii;
    uint16 s = 0;
	struct PRM p1;
	int ii, jj, nffti, nc;
	double bc, bw, cf, fh, fl;
	double rng_samp_rate, chirp_slope, pulse_dur, rng_bandwidth, wavelength;
	double SPEED_OF_LIGHT = 299792458.0;
	double *filterh, *filterl;
	short *buf1;
	fcomplex *c1, *c1h, *c1l, *tmp, fm;
	double f_h=0, f_l=0, w_h=0, w_l=0;
	
	void *API = NULL; /* GMT API control structure */
	// double complex zl,zh,zp,z1h,z2h,z1l,z2l;
	// double *ph;
    int half_band = 0;

	if (argc1 != 2 && argc1 != 3)
		die("", USAGE1);
	if ((API = GMT_Create_Session(argv1[0], 0U, 0U, NULL)) == NULL)
		return EXIT_FAILURE;
    if (argc1 == 3) half_band = atoi(argv1[2]);

	/* read in prm files */
	get_prm(&p1, argv1[1]);

	rng_samp_rate = p1.fs;
	chirp_slope = p1.chirp_slope;
	pulse_dur = p1.pulsedur;
	rng_bandwidth = fabs(pulse_dur * chirp_slope);
	wavelength = p1.lambda;
	cf = SPEED_OF_LIGHT / wavelength;

    if (p1.SC_identity == 10) {
        bc = 42000000.0 / 3.0;
    }
    else if (half_band == 1) {
        bc = rng_bandwidth / 3.0/2;
    }
    else {
	    bc = rng_bandwidth / 3.0;
    }
	bw = bc;

	fh = cf + bc;
	fl = cf - bc;

	nffti = find_fft_length(p1.num_rng_bins);
	nc = (int)fabs(round(bc / rng_samp_rate * nffti));

	filterh = (double *)malloc(nffti * sizeof(double));
	filterl = (double *)malloc(nffti * sizeof(double));
	cos_window(bc, bw, rng_samp_rate, nffti, filterh);
	cos_window(-bc, bw, rng_samp_rate, nffti, filterl);

	//fprintf(stderr, "%.12f %.12f %.12f %.12f %.12f\n", bc / 1e6, bw / 1e6, fh / 1e6, fl / 1e6, cf / 1e6);
	// test on how the window look like
	// for (ii=0;ii<p1.num_rng_bins;ii++){
	//    fprintf(stdout,"%.12f\n",filterl[ii]);
	//}
	// exit(0);

	// read in SLC and run split spectrum
    if (p1.SC_identity != 10) {
	    if ((SLC_file1 = fopen(p1.SLC_file, "rb")) == NULL)
		    die("Can't open ", p1.SLC_file);
	    // if ((SLC_file2 = fopen(p2.SLC_file,"rb")) == NULL) die("Can't open
	    // ",p2.input_file); if ((p_out = fopen("phase_output","wb")) == NULL)
	    // die("Can't open ","phase_output");

	    if ((SLCH = fopen("SLCH", "wb")) == NULL)
		    die("Can't open ", "SLCH");
	    // if ((SLCH2 = fopen("SLCH2","wb")) == NULL) die("Can't open ","SLCH2");
	    if ((SLCL = fopen("SLCL", "wb")) == NULL)
		    die("Can't open ", "SLCL");
	    // if ((SLCL2 = fopen("SLCL2","wb")) == NULL) die("Can't open ","SLCL2");
    }
    else {
        TIFFSetWarningHandler(NULL);
        if ((tif = TIFFOpen(p1.input_file,"rb")) == NULL)
            die("Couldn't open tiff file: \n", p1.input_file);
        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
        if ((tifh = TIFFOpen("high.tiff","wb")) == NULL)
            die("Couldn't open tiff file: \n", "high.tiff");
        TIFFSetField(tifh, TIFFTAG_IMAGEWIDTH, width);
        TIFFSetField(tifh, TIFFTAG_IMAGELENGTH, height);
        TIFFSetField(tifh, TIFFTAG_BITSPERSAMPLE, sizeof(short) * 8 * 2);
        TIFFSetField(tifh, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_COMPLEXINT);
        TIFFSetField(tifh, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
        if ((tifl = TIFFOpen("low.tiff","wb")) == NULL)
            die("Couldn't open tiff file: \n", "low.tiff");
        TIFFSetField(tifl, TIFFTAG_IMAGEWIDTH, width);
        TIFFSetField(tifl, TIFFTAG_IMAGELENGTH, height);
        TIFFSetField(tifl, TIFFTAG_BITSPERSAMPLE, sizeof(short) * 8 * 2);
        TIFFSetField(tifl, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_COMPLEXINT);
        TIFFSetField(tifl, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    }
        
    if (p1.SC_identity != 10) {
	    buf1 = (short *)malloc(nffti * 2 * sizeof(short));
    }
    else {
        buf1 = (short *)_TIFFmalloc(TIFFScanlineSize(tif) * 2);
    }
	// buf2 = (short *) malloc (nffti*2*sizeof(short));
	c1 = (fcomplex *)malloc(nffti * sizeof(fcomplex));
	// c2 = (fcomplex *) malloc (nffti*sizeof(fcomplex));
	c1h = (fcomplex *)malloc(nffti * sizeof(fcomplex));
	c1l = (fcomplex *)malloc(nffti * sizeof(fcomplex));
	// c2h = (fcomplex *) malloc (nffti*sizeof(fcomplex));
	// c2l = (fcomplex *) malloc (nffti*sizeof(fcomplex));
	tmp = (fcomplex *)malloc(nffti * sizeof(fcomplex));
	// ph = (double *) malloc (p1.num_rng_bins*sizeof(double));

	for (ii = 0; ii < nffti; ii++) {
		buf1[ii * 2] = 0;
		buf1[ii * 2 + 1] = 0;
		// buf2[ii*2] = 0;
		// buf2[ii*2+1] = 0;
		c1[ii].r = 0.0;
		c1[ii].i = 0.0;
		// c2[ii].r = 0.0; c2[ii].i = 0.0;
		c1h[ii].r = 0.0;
		c1h[ii].i = 0.0;
		c1l[ii].r = 0.0;
		c1l[ii].i = 0.0;
		// c2h[ii].r = 0.0; c2h[ii].i = 0.0;
		// c2l[ii].r = 0.0; c2l[ii].i = 0.0;
	}

	// fprintf(stderr,"%.6f
	// %.6f\n",cimag(cpow((1.0+2.0*I),(2.0-1.0*I))),creal(cpow((1.0+2.0*I),(2.0-1.0*I))));

	fprintf(stderr, "Number of NFFT is %d\n", nffti);

	fprintf(stderr, "Writing lines ");
    if (p1.SC_identity != 10) height = p1.num_valid_az * p1.num_patches;
	for (ii = 0; ii < height; ii++) {
        if (p1.SC_identity != 10) {
		    fread(buf1, sizeof(short), p1.num_rng_bins * 2, SLC_file1);
        }
        else {
            nii = ii;
            TIFFReadScanline(tif, buf1, nii, s);
        }
        fm = fmean(buf1,p1.num_rng_bins);
		// fread(buf2,sizeof(short),p1.num_rng_bins*2,SLC_file2);
		for (jj = 0; jj < nffti; jj++) {
			if (jj < p1.num_rng_bins) {
				c1[jj].r = (float)(buf1[2 * jj]) - fm.r;
				c1[jj].i = (float)(buf1[2 * jj + 1]) - fm.i;
				// c2[jj].r = (float)(buf2[2*jj]);
				// c2[jj].i = (float)(buf2[2*jj+1]);
			}
			else {
				c1[jj].r = 0.0;
				c1[jj].i = 0.0;
				// c2[jj].r = 0.0;
				// c2[jj].i = 0.0;
			}
		}

		// 1-D fourier transform
		GMT_FFT_1D(API, (float *)c1, nffti, GMT_FFT_FWD, GMT_FFT_COMPLEX);
		// GMT_FFT_1D (API, (float *)c2, nffti, GMT_FFT_FWD, GMT_FFT_COMPLEX);

		// assign them to arrays with band passing
		for (jj = 0; jj < nffti; jj++) {
			c1h[jj].r = c1[jj].r * filterh[jj];
			c1h[jj].i = c1[jj].i * filterh[jj];
			c1l[jj].r = c1[jj].r * filterl[jj];
			c1l[jj].i = c1[jj].i * filterl[jj];
			
			w_h += c1h[jj].r*c1h[jj].r+c1h[jj].i*c1h[jj].i;
			w_l += c1l[jj].r*c1l[jj].r+c1l[jj].i*c1l[jj].i;
			f_h += (c1h[jj].r*c1h[jj].r+c1h[jj].i*c1h[jj].i)*p1.fs/nffti*jj;
			f_l += (c1l[jj].r*c1l[jj].r+c1l[jj].i*c1l[jj].i)*p1.fs/nffti*(jj-nffti);
			
			// c2h[jj].r = c2[jj].r*filterh[jj];
			// c2h[jj].i = c2[jj].i*filterh[jj];
			// c2l[jj].r = c2[jj].r*filterl[jj];
			// c2l[jj].i = c2[jj].i*filterl[jj];

			/*
			 c1h[jj].r = c1[jj].r;
			 c1h[jj].i = c1[jj].i;
			 c1l[jj].r = c1[jj].r;
			 c1l[jj].i = c1[jj].i;
			 c2h[jj].r = c2[jj].r;
			 c2h[jj].i = c2[jj].i;
			 c2l[jj].r = c2[jj].r;
			 c2l[jj].i = c2[jj].i;
			 */
		}

		// shift back to the center
		circ_shift(c1h, tmp, nffti, -nc);
		circ_shift(c1l, tmp, nffti, nc);
		// circ_shift(c2h,tmp,nffti,-nc);
		// circ_shift(c2l,tmp,nffti,nc);

		// 1-D inverse fourier transform
		GMT_FFT_1D(API, (float *)c1h, nffti, GMT_FFT_INV, GMT_FFT_COMPLEX);
		GMT_FFT_1D(API, (float *)c1l, nffti, GMT_FFT_INV, GMT_FFT_COMPLEX);
		// GMT_FFT_1D (API, (float *)c2h, nffti, GMT_FFT_INV, GMT_FFT_COMPLEX);
		// GMT_FFT_1D (API, (float *)c2l, nffti, GMT_FFT_INV, GMT_FFT_COMPLEX);

		// compute ionospheric phase ?
		for (jj = 0; jj < p1.num_rng_bins; jj++) {
			// the commented part won't work in cases ionospheric signal is big
			// fl*fh/(f_H^2-f_L^2)*angle(C_L^(f_H/f_0)*conj(C_H^(f_L/f_0)))

			// zh = (c1h[jj].r*c2h[jj].r + c1h[jj].i*c2h[jj].i) + (c1h[jj].i*c2h[jj].r
			// - c1h[jj].r*c2h[jj].i) * I; zl = (c1l[jj].r*c2l[jj].r +
			// c1l[jj].i*c2l[jj].i) + (c1l[jj].i*c2l[jj].r - c1l[jj].r*c2l[jj].i) * I;
			// z1h = c1h[jj].r + c1h[jj].i * I;
			// z2h = c2h[jj].r + c2h[jj].i * I;
			// z1l = c1l[jj].r + c1l[jj].i * I;
			// z2l = c2l[jj].r + c2l[jj].i * I;
			// zh = z1h*conj(z2h);
			// zl = z1l*conj(z2l);
			// zp = cpow(zl,fh/cf+0.0*I)*conj(pow(zh,fl/cf+0.0*I));
			// ph[jj] = fl*fh/(fh*fh-fl*fl)*atan2(creal(zp),cimag(zp));

			// ph[jj] = fmod(ph[jj]+PI,PI)-PI;

			buf1[2 * jj] = (short)round(c1h[jj].r);
			buf1[2 * jj + 1] = (short)round(c1h[jj].i);
			// buf2[2*jj] = (short)round(c2h[jj].r);
			// buf2[2*jj+1] = (short)round(c2h[jj].i);
		}
        
        if (p1.SC_identity != 10) {
		    fwrite((void *)buf1, sizeof(short), p1.num_rng_bins * 2, SLCH);
        }
        else {
            TIFFWriteScanline(tifh, buf1, nii, s);
        }
        
		// fwrite((void *)buf2,sizeof(short),p1.num_rng_bins*2,SLCH2);
		for (jj = 0; jj < p1.num_rng_bins; jj++) {
			buf1[2 * jj] = (int)round(c1l[jj].r);
			buf1[2 * jj + 1] = (int)round(c1l[jj].i);
			// buf2[2*jj] = (int)round(c2l[jj].r);
			// buf2[2*jj+1] = (int)round(c2l[jj].i);
		}
        if (p1.SC_identity != 10) {
		    fwrite((void *)buf1, sizeof(short), p1.num_rng_bins * 2, SLCL);
        }
        else {
            TIFFWriteScanline(tifl, buf1, nii, s);
        }
    
		// fwrite((void *)buf2,sizeof(short),p1.num_rng_bins*2,SLCL2);

		// fwrite((void *)ph,sizeof(double),p1.num_rng_bins,p_out);

		if (ii % 1000 == 0)
			fprintf(stderr, "%d ", ii);
	}
	fprintf(stderr, "...\n");
	if (p1.SC_identity == 10) {
		printf("low_wavelength = %.12f\n",SOL/(cf+f_l/w_l));
    	printf("center_wavelength = %.12f\n",SOL/cf);
    	printf("high_wavelength = %.12f\n",SOL/(cf+f_h/w_h));
		printf("low_freq = %.12f\n",cf+f_l/w_l);
    	printf("center_freq = %.12f\n",cf);
    	printf("high_freq = %.12f\n",cf+f_h/w_h);
	}
	else {
		printf("low_wavelength = %.12f\n",SOL/fl);
    	printf("center_wavelength = %.12f\n",SOL/cf);
    	printf("high_wavelength = %.12f\n",SOL/fh);
    	printf("low_freq = %.12f\n",fl);
    	printf("center_freq = %.12f\n",cf);
    	printf("high_freq = %.12f\n",fh);
    }
    
    printf("low_bandwidth = %.12f\n",bw);
    printf("center_bandwidth = %.12f\n",rng_bandwidth);
    printf("high_bandwidth = %.12f\n",bw);

	free(filterh);
	free(filterl);
	
	// free(buf2);
	free(c1);
	// free(c2);
	free(c1h);
	free(c1l);
	// free(c2h);
	// free(c2l);
	// free(ph);
    if (p1.SC_identity != 10) {
	    fclose(SLC_file1);
	    // fclose(SLC_file2);
	    // fclose(p_out);
	    fclose(SLCH);
        fclose(SLCL);
        free(buf1);
    }
    else {
        TIFFClose(tif);
        TIFFClose(tifh);
        TIFFClose(tifl);
        free(buf1);
    }
	// fclose(SLCH2);
	// fclose(SLCL2);
	free(tmp);
	return (1);
}

char *USAGE2 = "\nUSAGE: split_spectrum prm1 prm2\n\n"
               "   program used to split range spectrum for coregistered SLC using a "
               "modified cosine filter\n\n"
               "   SLCs are bandpassed and then shifted to the center of the spectrum\n\n"
               "   output is SLCH1 SLCH2 SLCL1 SLCL2\n";

int split2(int argc2, char **argv2) {

	FILE *SLC_file1, *SLC_file2, *SLCH1, *SLCH2, *SLCL1, *SLCL2;
	// FILE *p_out;
	struct PRM p1, p2;
	int ii, jj, nffti, nc;
	double bc, bw, cf, fh, fl;
	double rng_samp_rate, chirp_slope, pulse_dur, rng_bandwidth, wavelength;
	double SPEED_OF_LIGHT = 299792458.0;
	double *filterh, *filterl;
	short *buf1, *buf2;
	fcomplex *c1, *c2, *c1h, *c1l, *c2h, *c2l, *tmp;
	void *API = NULL; /* GMT API control structure */
	// double complex zl,zh,zp,z1h,z2h,z1l,z2l;
	// double *ph;

	if (argc2 != 3)
		die("", USAGE2);
	if ((API = GMT_Create_Session(argv2[0], 0U, 0U, NULL)) == NULL)
		return EXIT_FAILURE;

	/* read in prm files */
	get_prm(&p1, argv2[1]);
	get_prm(&p2, argv2[2]);

	rng_samp_rate = p1.fs;
	chirp_slope = p1.chirp_slope;
	pulse_dur = p1.pulsedur;
	rng_bandwidth = fabs(pulse_dur * chirp_slope);
	wavelength = p1.lambda;
	cf = SPEED_OF_LIGHT / wavelength;

	bc = rng_bandwidth / 3.0;
	bw = bc;

	fh = cf + bc;
	fl = cf - bc;

	nffti = find_fft_length(p1.num_rng_bins);
	nc = (int)fabs(round(bc / rng_samp_rate * nffti));

	filterh = (double *)malloc(nffti * sizeof(double));
	filterl = (double *)malloc(nffti * sizeof(double));
	cos_window(bc, bw, rng_samp_rate, nffti, filterh);
	cos_window(-bc, bw, rng_samp_rate, nffti, filterl);

	fprintf(stderr, "%.12f %.12f %.12f %.12f %.12f\n", bc / 1e6, bw / 1e6, fh / 1e6, fl / 1e6, cf / 1e6);
	// test on how the window look like
	// for (ii=0;ii<p1.num_rng_bins;ii++){
	//    fprintf(stdout,"%.12f\n",filterl[ii]);
	//}
	// exit(0);

	// read in SLC and run split spectrum
	if ((SLC_file1 = fopen(p1.SLC_file, "rb")) == NULL)
		die("Can't open ", p1.input_file);
	if ((SLC_file2 = fopen(p2.SLC_file, "rb")) == NULL)
		die("Can't open ", p2.input_file);
	// if ((p_out = fopen("phase_output","wb")) == NULL) die("Can't open
	// ","phase_output");

	if ((SLCH1 = fopen("SLCH1", "wb")) == NULL)
		die("Can't open ", "SLCH1");
	if ((SLCH2 = fopen("SLCH2", "wb")) == NULL)
		die("Can't open ", "SLCH2");
	if ((SLCL1 = fopen("SLCL1", "wb")) == NULL)
		die("Can't open ", "SLCL1");
	if ((SLCL2 = fopen("SLCL2", "wb")) == NULL)
		die("Can't open ", "SLCL2");

	buf1 = (short *)malloc(nffti * 2 * sizeof(short));
	buf2 = (short *)malloc(nffti * 2 * sizeof(short));
	c1 = (fcomplex *)malloc(nffti * sizeof(fcomplex));
	c2 = (fcomplex *)malloc(nffti * sizeof(fcomplex));
	c1h = (fcomplex *)malloc(nffti * sizeof(fcomplex));
	c1l = (fcomplex *)malloc(nffti * sizeof(fcomplex));
	c2h = (fcomplex *)malloc(nffti * sizeof(fcomplex));
	c2l = (fcomplex *)malloc(nffti * sizeof(fcomplex));
	tmp = (fcomplex *)malloc(nffti * sizeof(fcomplex));
	// ph = (double *) malloc (p1.num_rng_bins*sizeof(double));

	for (ii = 0; ii < nffti; ii++) {
		buf1[ii * 2] = 0;
		buf1[ii * 2 + 1] = 0;
		buf2[ii * 2] = 0;
		buf2[ii * 2 + 1] = 0;
		c1[ii].r = 0.0;
		c1[ii].i = 0.0;
		c2[ii].r = 0.0;
		c2[ii].i = 0.0;
		c1h[ii].r = 0.0;
		c1h[ii].i = 0.0;
		c1l[ii].r = 0.0;
		c1l[ii].i = 0.0;
		c2h[ii].r = 0.0;
		c2h[ii].i = 0.0;
		c2l[ii].r = 0.0;
		c2l[ii].i = 0.0;
	}

	// fprintf(stderr,"%.6f
	// %.6f\n",cimag(cpow((1.0+2.0*I),(2.0-1.0*I))),creal(cpow((1.0+2.0*I),(2.0-1.0*I))));

	fprintf(stderr, "Number of NFFT is %d\n", nffti);

	fprintf(stderr, "Writing lines ");
	for (ii = 0; ii < p1.num_valid_az * p1.num_patches; ii++) {
		fread(buf1, sizeof(short), p1.num_rng_bins * 2, SLC_file1);
		fread(buf2, sizeof(short), p1.num_rng_bins * 2, SLC_file2);
		for (jj = 0; jj < nffti; jj++) {
			if (jj < p1.num_rng_bins) {
				c1[jj].r = (float)(buf1[2 * jj]);
				c1[jj].i = (float)(buf1[2 * jj + 1]);
				c2[jj].r = (float)(buf2[2 * jj]);
				c2[jj].i = (float)(buf2[2 * jj + 1]);
			}
			else {
				c1[jj].r = 0.0;
				c1[jj].i = 0.0;
				c2[jj].r = 0.0;
				c2[jj].i = 0.0;
			}
		}

		// 1-D fourier transform
		GMT_FFT_1D(API, (float *)c1, nffti, GMT_FFT_FWD, GMT_FFT_COMPLEX);
		GMT_FFT_1D(API, (float *)c2, nffti, GMT_FFT_FWD, GMT_FFT_COMPLEX);

		// assign them to arrays with band passing
		for (jj = 0; jj < nffti; jj++) {
			c1h[jj].r = c1[jj].r * filterh[jj];
			c1h[jj].i = c1[jj].i * filterh[jj];
			c1l[jj].r = c1[jj].r * filterl[jj];
			c1l[jj].i = c1[jj].i * filterl[jj];
			c2h[jj].r = c2[jj].r * filterh[jj];
			c2h[jj].i = c2[jj].i * filterh[jj];
			c2l[jj].r = c2[jj].r * filterl[jj];
			c2l[jj].i = c2[jj].i * filterl[jj];
			/*
			c1h[jj].r = c1[jj].r;
			c1h[jj].i = c1[jj].i;
			c1l[jj].r = c1[jj].r;
			c1l[jj].i = c1[jj].i;
			c2h[jj].r = c2[jj].r;
			c2h[jj].i = c2[jj].i;
			c2l[jj].r = c2[jj].r;
			c2l[jj].i = c2[jj].i;
			*/
		}

		// shift back to the center
		circ_shift(c1h, tmp, nffti, -nc);
		circ_shift(c1l, tmp, nffti, nc);
		circ_shift(c2h, tmp, nffti, -nc);
		circ_shift(c2l, tmp, nffti, nc);

		// 1-D inverse fourier transform
		GMT_FFT_1D(API, (float *)c1h, nffti, GMT_FFT_INV, GMT_FFT_COMPLEX);
		GMT_FFT_1D(API, (float *)c1l, nffti, GMT_FFT_INV, GMT_FFT_COMPLEX);
		GMT_FFT_1D(API, (float *)c2h, nffti, GMT_FFT_INV, GMT_FFT_COMPLEX);
		GMT_FFT_1D(API, (float *)c2l, nffti, GMT_FFT_INV, GMT_FFT_COMPLEX);

		// compute ionospheric phase
		for (jj = 0; jj < p1.num_rng_bins; jj++) {
			// the commented part won't work in cases ionospheric signal is big
			// fl*fh/(f_H^2-f_L^2)*angle(C_L^(f_H/f_0)*conj(C_H^(f_L/f_0)))

			// zh = (c1h[jj].r*c2h[jj].r + c1h[jj].i*c2h[jj].i) + (c1h[jj].i*c2h[jj].r
			// - c1h[jj].r*c2h[jj].i) * I; zl = (c1l[jj].r*c2l[jj].r +
			// c1l[jj].i*c2l[jj].i) + (c1l[jj].i*c2l[jj].r - c1l[jj].r*c2l[jj].i) * I;
			// z1h = c1h[jj].r + c1h[jj].i * I;
			// z2h = c2h[jj].r + c2h[jj].i * I;
			// z1l = c1l[jj].r + c1l[jj].i * I;
			// z2l = c2l[jj].r + c2l[jj].i * I;
			// zh = z1h*conj(z2h);
			// zl = z1l*conj(z2l);
			// zp = cpow(zl,fh/cf+0.0*I)*conj(pow(zh,fl/cf+0.0*I));
			// ph[jj] = fl*fh/(fh*fh-fl*fl)*atan2(creal(zp),cimag(zp));

			// ph[jj] = fmod(ph[jj]+PI,PI)-PI;

			buf1[2 * jj] = (short)round(c1h[jj].r);
			buf1[2 * jj + 1] = (short)round(c1h[jj].i);
			buf2[2 * jj] = (short)round(c2h[jj].r);
			buf2[2 * jj + 1] = (short)round(c2h[jj].i);
		}
		fwrite((void *)buf1, sizeof(short), p1.num_rng_bins * 2, SLCH1);
		fwrite((void *)buf2, sizeof(short), p1.num_rng_bins * 2, SLCH2);
		for (jj = 0; jj < p1.num_rng_bins; jj++) {
			buf1[2 * jj] = (int)round(c1l[jj].r);
			buf1[2 * jj + 1] = (int)round(c1l[jj].i);
			buf2[2 * jj] = (int)round(c2l[jj].r);
			buf2[2 * jj + 1] = (int)round(c2l[jj].i);
		}
		fwrite((void *)buf1, sizeof(short), p1.num_rng_bins * 2, SLCL1);
		fwrite((void *)buf2, sizeof(short), p1.num_rng_bins * 2, SLCL2);

		// fwrite((void *)ph,sizeof(double),p1.num_rng_bins,p_out);

		if (ii % 1000 == 0)
			fprintf(stderr, "%d ", ii);
	}
	fprintf(stderr, "...\n");

	free(filterh);
	free(filterl);
	free(buf1);
	free(buf2);
	free(c1);
	free(c2);
	free(c1h);
	free(c1l);
	free(c2h);
	free(c2l);
	// free(ph);
	fclose(SLC_file1);
	fclose(SLC_file2);
	// fclose(p_out);
	fclose(SLCH1);
	fclose(SLCL1);
	fclose(SLCH2);
	fclose(SLCL2);
	free(tmp);
	return (1);
}

fcomplex fmean(short *c, int N) {

    fcomplex m;
    int i;
    
    m.r = 0.0;
    m.i = 0.0;
    for (i=0;i<N;i++) {
        m.r += (float)c[2*i];
        m.i += (float)c[2*i+1];
    }
    m.r = m.r/(float)N;
    m.i = m.i/(float)N;
  
   return(m);
}

int cos_window(double fc, double fb, double fs, int N, double *filter) {
	// create a cosine window filter in fourier domain

	int i, nc, nb, flat_nb, cos_nb;

	if (fabs(fc) > fs / 3) {
		die("center frequency too big!", "");
	}

	for (i = 0; i < N; i++)
		filter[i] = 0.0;

	nc = (int)fabs(round(fc / fs * N));
	nb = (int)round(fb / fs * N);
	// fprintf(stderr,"%d  %d %d\n",nc,nb,N);
	// modify the section below to tune the band filter parameters ,
	// 2*flat_nb+2*cos_nb = 1
	flat_nb = (int)round(nb / 4);
	cos_nb = (int)round(nb / 4);

	for (i = nc - flat_nb - 1; i < nc + flat_nb; i++)
		filter[i] = 1;
	for (i = nc - flat_nb - cos_nb - 1; i < nc - flat_nb - 1; i++)
		filter[i] = 0.5 - cos(PI * (i - (nc - flat_nb - cos_nb - 1) + 1) / cos_nb) / 2.0;
	for (i = nc + flat_nb; i < nc + flat_nb + cos_nb; i++)
		filter[i] = cos(PI * (i - (nc + flat_nb) + 1) / cos_nb) / 2.0 + 0.5;

	if (fc < 0) {
		fliplr(filter, N);
	}
	return (1);
}

int fliplr(double *a, int N) {
	int i;
	double d;

	for (i = 0; i < (int)(N / 2); i++) {
		d = a[i];
		a[i] = a[N - 1 - i];
		a[N - 1 - i] = d;
	}

	return (1);
}

void circ_shift(fcomplex *in, fcomplex *tmp, int na, int nc) {

	int i, ii;
	// int ncm;
	// ncm = nc%na;

	if (nc % na != 0) {

		for (i = 0; i < na; i++) {
			tmp[i].r = in[i].r;
			tmp[i].i = in[i].i;
		}

		for (i = 0; i < na; i++) {
			ii = (i + nc) % na;
			if (ii < 0)
				ii = ii + na;
			if (ii < 0 || ii > na) {
				fprintf(stderr, "%d", ii);
				die("index outside of range", "");
			}
			in[ii].r = tmp[i].r;
			in[ii].i = tmp[i].i;
		}
		// for(i = 0; i < abs(ncm); i++)
		//    left_shift(in, na);
	}
}

/*

// modified bessel function
double bessi0(double x) {

    double ax,ans;
    double y;

    if ((ax=fabs(x)) < 3.75) {
        y=x/3.75;
        y*=y;
        ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
    } else {
        y=3.75/ax;
        ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1+y*0.392377e-2))))))));
    }

    return ans;
}

// kaiser window, maximizes energy concentration
void kaiser(double alpha, int n, double *coef){

    int i;
    int hn;
    double a;

    hn = (n - 1) / 2;

    for(i = -hn; i<=hn; i++){
        a = 1.0 - 4.0 * i * i / (n - 1.0) / (n - 1.0);
        coef[i+hn] = bessi0(alpha * sqrt(a)) / bessi0(alpha);
    }
}

// create a band pass filter
void bandpass_filter(double bw, double bc, int n, int nfft, int ncshift, double
alpha, fcomplex *filter) {

    int i;
    double *kw;
    int hn;
    fcomplex bwx,bcx;

    hn = (n-1)/2;
    if(n > nfft){
        die("Error: fft length too small!\n","");
    }
    if(abs(ncshift) > nfft){
        die("","Error: fft length too small or shift too big!\n\n");
    }

    //set filter to zero
    for(i = 0; i < nfft; i++){
        filter[i].r = 0.0;
        filter[i].i = 0.0;
    }


    kw = (double *) malloc ((hn*2+1) * sizeof(double));
    kaiser(alpha,n,kw);

    //compute filter
    for(i=0;i<2*hn+1;i++) {
        bcx.r = cos(bc*2.0*PI*(i-hn));
        bcx.i = sin(bc*2.0*PI*(i-hn));

        if(i==hn+1){
            bwx.r = 1.0;
            bwx.i = 0.0;
        }
        else{
            bwx.r = sin(bw*PI*(i-hn))/(bw*PI*(i-hn));
            bwx.i = 0.0;
        }
        filter[i] = Cmul(bcx,bwx);
        filter[i].r = bw*kw[i]*filter[i].r;
        filter[i].i = bw*kw[i]*filter[i].i;
    }

    // shift the filter
    circ_shift(filter, nfft, -abs(ncshift));

    free(kw);
}

void circ_shift(fcomplex *in, int na, int nc){

    int i;
    int ncm;

    ncm = nc%na;

    if(ncm < 0){
        for(i = 0; i < abs(ncm); i++)
            left_shift(in, na);
    }
    else if(ncm > 0){
        for(i = 0; i < ncm; i++)
            right_shift(in, na);
    }
    else{ //ncm == 0, no need to shift
        i = 0;
    }
}

void left_shift(fcomplex *in, int na){

    int i;
    fcomplex x;

    if(na < 1){
        fprintf(stderr, "Error: array size < 1\n\n");
        exit(1);
    }
    else if(na > 1){
        x.r = in[0].r;
        x.i = in[0].i;
        for(i = 0; i <= na - 2; i++){
            in[i].r = in[i+1].r;
            in[i].i = in[i+1].i;
        }
        in[na-1].r = x.r;
        in[na-1].i = x.i;
    }
    else{ //na==1, no need to shift
        i = 0;
    }
}

void right_shift(fcomplex *in, int na){

    int i;
    fcomplex x;

    if(na < 1){
        fprintf(stderr, "Error: array size < 1\n\n");
        exit(1);
    }
    else if(na > 1){
        x.r = in[na-1].r;
        x.i = in[na-1].i;
        for(i = na - 1; i >= 1; i--){
            in[i].r = in[i-1].r;
            in[i].i = in[i-1].i;
        }
        in[0].r = x.r;
        in[0].i = x.i;
    }
    else{ //na==1, no need to shift
        i = 0;
    }
}

*/
