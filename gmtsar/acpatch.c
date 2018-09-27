/************************************************************************
 * acpatch is a subroutine to perform azimuth compression of radar echos.*
 * It is part of the esar SAR processor.  Echos should be range          *
 * compressed and migrated before applying azimuth compression.          *
 ************************************************************************/
/************************************************************************
 * Creator: Evelyn J. Price	(Scripps Institution of Oceanography)   *
 * Date   : 11/18/96                                                     *
 ************************************************************************/
/************************************************************************
 * Modification history:                                                 *
 * 34MAR2000 - Modified to account for azimuth stretch with azimuth      *
 *                                                                       *
 *                                                                       *
 ************************************************************************/
#include "gmt.h"
#include "siocomplex.h"
#include "soi.h"

void acpatch(void *API, fcomplex **data, int nrows, double delr, double fd, double fdd, double fddd) {
	int i, j, k;
	int n, nfc, nf0;
	int *np;
	static int firsttime = 1;
	double dxsamp1, a2, a4, phase, rd0;
	double *r, dx, v1, az, y2;
	float t, sinsq;
	float *y, *f0, *f_rate;
	fcomplex cpha, *ref, *fft_vec;

	/* allocate memory */

	if ((np = (int *)malloc(num_rng_bins * sizeof(int))) == NULL) {
		printf("sorry, couldn't allocate mem for np.\n");
		exit(-1);
	}
	if ((r = (double *)malloc(num_rng_bins * sizeof(double))) == NULL) {
		printf("sorry, couldn't allocate mem for r.\n");
		exit(-1);
	}
	if ((y = (float *)malloc(num_rng_bins * sizeof(float))) == NULL) {
		printf("sorry, couldn't allocate mem for y.\n");
		exit(-1);
	}
	if ((f0 = (float *)malloc(num_rng_bins * sizeof(float))) == NULL) {
		printf("sorry, couldn't allocate mem for f0.\n");
		exit(-1);
	}
	if ((f_rate = (float *)malloc(num_rng_bins * sizeof(float))) == NULL) {
		printf("sorry, couldn't allocate mem for f_rate.\n");
		exit(-1);
	}
	if ((ref = (fcomplex *)malloc(nrows * sizeof(fcomplex))) == NULL) {
		printf("sorry, couldn't allocate mem for ref.\n");
		exit(-1);
	}
	if ((fft_vec = (fcomplex *)malloc(nrows * sizeof(fcomplex))) == NULL) {
		printf("sorry, couldn't allocate mem for fft_vec.\n");
		exit(-1);
	}

	/* this is the actual x_pix_size that we want to output */
	dxsamp1 = vel1 / prf1;

	a4 = lambda / (2.0 * dxsamp1 * az_res);

	if ((strcmp(deskew, "y") == 0) || (strcmp(deskew, "Y") == 0))
		a2 = -2.0 * PI / (dxsamp1 * (float)nrows);
	else
		a2 = 0.0;

	dx = vel1 / prf1;
	v1 = pow(lambda, 2.0) / (8.0 * pow(dx, 2.0));

	/* convert a resampling coefficients to a function of range instead of pixel
	 */

	if (firsttime == 1) {
		sub_int_a = sub_int_a - stretch_a * near_range / delr;
		stretch_a = stretch_a / delr;
		firsttime = 0;
	}

	for (i = 0; i < num_rng_bins; i++) {
		r[i] = near_range + ((float)i) * delr;
		f0[i] = fd + fdd * r[i] + fddd * r[i] * r[i];
		rd0 = r[i] / (1 + v1 * pow((f0[i] / prf1), 2.0));
		f_rate[i] = -2.0 * pow(vel1, 2.0) * pow((rd0 / r[i]), 2.0) / (lambda * r[i]);
		np[i] = (int)(r[i] * a4 / 2);
		az = stretch_a * r[i] + sub_int_a;
		y2 = PI2 * az / ((float)nrows);

		sinsq = lambda * f0[i] / (2.0 * vel1);
		y[i] = r[i] * a2 * sinsq + y2;

		/* create reference function */
		for (j = 0; j < nrows; j++) {
			ref[j].r = 0.0f;
			ref[j].i = 0.0f;
		}

		phase = PI * pow(f0[i], 2.0) / f_rate[i];
		ref[0] = Cexp(phase);
		ref[0] = RCmul((1.0 / nrows), ref[0]);
		for (j = 0; j < np[i]; j++) {
			t = ((float)(j + 1)) / prf1;
			phase = PI * f_rate[i] * t * t + PI2 * f0[i] * t;
			ref[j + 1] = Cexp(phase);
			ref[j + 1] = RCmul((1.0 / nrows), ref[j + 1]);
			phase = PI * f_rate[i] * t * t - PI2 * f0[i] * t;
			ref[-j + nrows - 1] = Cexp(phase);
			ref[-j + nrows - 1] = RCmul((1.0 / nrows), ref[-j + nrows - 1]);
		}

		/*  transform the reference */
		// dir = -1;
		// cfft1d_(&nrows,ref,&dir);
		GMT_FFT_1D(API, (float *)ref, nrows, GMT_FFT_FWD, GMT_FFT_COMPLEX);
		/*  multiply the reference by the data */
		n = (int)((f0[i] / prf1) + 0.5);
		nf0 = nrows * (f0[i] - n * prf1) / prf1;
		nfc = nf0 + nrows / 2;
		if (nfc > nrows)
			nfc = nfc - nrows;
		phase = -y[i] * nf0;
		for (k = 0; k < nfc; k++) {
			ref[k] = Conjg(ref[k]);
			data[k][i] = Cmul(data[k][i], ref[k]);
			cpha = Cexp(phase);
			data[k][i] = Cmul(data[k][i], cpha);
			phase = phase + y[i];
		}
		phase = -y[i] * nf0;
		for (k = nrows - 1; k >= nfc; k--) {
			ref[k] = Conjg(ref[k]);
			data[k][i] = Cmul(data[k][i], ref[k]);
			cpha = Cexp(phase);
			data[k][i] = Cmul(data[k][i], cpha);
			phase = phase - y[i];
		}
		/*  inverse transform the product */
		for (j = 0; j < nrows; j++) {
			fft_vec[j] = data[j][i];
		}
		// dir = 1;
		// cfft1d_(&nrows,fft_vec,&dir);
		GMT_FFT_1D(API, (float *)fft_vec, nrows, GMT_FFT_INV, GMT_FFT_COMPLEX);
		for (j = 0; j < nrows; j++) {
			data[j][i] = fft_vec[j];
		}
	}
	free((char *)np);
	free((char *)r);
	free((char *)y);
	free((char *)f0);
	free((char *)f_rate);
	free((char *)ref);
	free((char *)fft_vec);
}
