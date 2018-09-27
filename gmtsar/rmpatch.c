/************************************************************************
 *  rmpatch performs range migration on range compressed radar echos. 	*
 ************************************************************************/
/************************************************************************
 * Creator: Evelyn J. Price	(Scripps Institution of Oceanography)   *
 * Date   : 11/18/96							*
 ************************************************************************/
/************************************************************************
 * Modification History							*
 * 24MAR2000 - modified to handle the range stretch with abs. azimuth    *
 * 									*
 * 									*
 ************************************************************************/

#include "lib_functions.h"
#include "siocomplex.h"
#include "soi.h"

void rmpatch(fcomplex **data, int nrows, double delr, double fd, double fdd, double fddd) {

	int *nvtmp;
	int na, i, j, ifrac, n, k, nfilter = 8192, ideskew;
	static int firsttime = 1;
	float *f0, *f_rate, *bdel, *vtmp;
	static float xintp[73728];
	float frac, ratio, freq;
	float c_xintp[8];
	double *r, *rd0, dx, v1, tmpd;
	fcomplex c_ctmpb[8], tmp[8], *c_ctmpa;

	/* initializations */
	if ((strcmp(deskew, "y") == 0) || (strcmp(deskew, "Y") == 0)) {
		ideskew = 1;
	}
	else {
		ideskew = 0;
	}

	if ((nvtmp = (int *)malloc(num_rng_bins * sizeof(int))) == NULL) {
		printf("Sorry, can't allocate memory for nvtmp.\n");
		exit(-1);
	}
	if ((rd0 = (double *)malloc(num_rng_bins * sizeof(double))) == NULL) {
		printf("Sorry, can't allocate memory for rd0.\n");
		exit(-1);
	}
	if ((f0 = (float *)malloc(num_rng_bins * sizeof(float))) == NULL) {
		printf("Sorry, can't allocate memory for f0.\n");
		exit(-1);
	}
	if ((f_rate = (float *)malloc(num_rng_bins * sizeof(float))) == NULL) {
		printf("Sorry, can't allocate memory for f_rate.\n");
		exit(-1);
	}
	if ((bdel = (float *)malloc(num_rng_bins * sizeof(float))) == NULL) {
		printf("Sorry, can't allocate memory for bdel.\n");
		exit(-1);
	}
	if ((vtmp = (float *)malloc(num_rng_bins * sizeof(float))) == NULL) {
		printf("Sorry, can't allocate memory for vtmp.\n");
		exit(-1);
	}
	if ((r = (double *)malloc(num_rng_bins * sizeof(double))) == NULL) {
		printf("Sorry, can't allocate memory for r.\n");
		exit(-1);
	}
	if ((c_ctmpa = (fcomplex *)malloc(num_rng_bins * sizeof(fcomplex))) == NULL) {
		printf("Sorry, can't allocate memory for rd0.\n");
		exit(-1);
	}

	/* load the interpolation array */
	/* convert r resampling coefficients to a function of range instead of pixel
	 */

	if (firsttime == 1) {
		intp_coef(nfilter, xintp);
		sub_int_r = sub_int_r - stretch_r * near_range / delr;
		stretch_r = stretch_r / delr;
		firsttime = 0;
	}

	dx = vel1 / prf1;
	v1 = pow(lambda, 2.0) / (8.0 * pow(dx, 2.0));
	for (i = 0; i < num_rng_bins; i++) {
		r[i] = near_range + i * delr;
		f0[i] = fd + fdd * r[i] + fddd * r[i] * r[i];
		rd0[i] = r[i] / (1 + v1 * pow((f0[i] / prf1), 2.0));
		f_rate[i] = -2.0 * pow(vel1, 2.0) * pow((rd0[i] / r[i]), 2.0) / (lambda * r[i]);
		bdel[i] = stretch_r * r[i] + sub_int_r;
	}

	for (na = 0; na < nrows; na++) {

		/* get the interpolation amounts for a given azimuth pixel na as f(line) */
		freq = ((float)na) / ((float)nrows) * prf1;
		for (i = 0; i < num_rng_bins; i++) {

			/*     frequencies must be within 0.5*prf of centroid */
			ratio = (freq - f0[i]) / prf1;
			n = (int)(ratio + 0.5);
			freq = freq - n * prf1;

			/*     range of a pixel at freq f, bdel is range correction for
			 * interferogram */

			if (ideskew == 1) {
				tmpd = bdel[i] + ((r[i] - (lambda / 4.0) * pow(f0[i], 2.0) / f_rate[i]) - r[0]) / delr +
				       rd0[i] * (v1 / delr) * (pow(freq, 2.0) - pow(f0[i], 2.0)) / pow(prf1, 2.0);
			}
			else {
				tmpd = i + rd0[i] * (v1 / delr) * (pow(freq, 2.0) - pow(f0[i], 2.0)) / pow(prf1, 2.0) + bdel[i];
			}
			nvtmp[i] = (int)tmpd;
			vtmp[i] = tmpd - nvtmp[i];
		}

		/*  interpolate that line according to coeffs determined above */
		for (i = 0; i < num_rng_bins; i++) {
			c_ctmpa[i].r = 0.0;
			c_ctmpa[i].i = 0.0;
			if ((nvtmp[i] >= 3) && (nvtmp[i] < (num_rng_bins - 5))) {
				frac = vtmp[i];
				ifrac = 8 * ((int)(frac * ((float)nfilter) + 0.5));
				for (k = 0; k < 8; k++) {
					c_xintp[k] = xintp[ifrac + k];
					c_ctmpb[k] = data[na][nvtmp[i] - 2 + k];
				}

				for (j = 0; j < 8; j++) {
					tmp[j] = RCmul(c_xintp[j], c_ctmpb[j]);
				}
				for (j = 0; j < 8; j++) {
					c_ctmpa[i] = Cadd(c_ctmpa[i], tmp[j]);
				}
			}
		}
		for (i = 0; i < num_rng_bins; i++) {
			data[na][i] = c_ctmpa[i];
		}
	}
	free((char *)nvtmp);
	free((char *)f0);
	free((char *)rd0);
	free((char *)f_rate);
	free((char *)bdel);
	free((char *)vtmp);
	free((char *)r);
	free((char *)c_ctmpa);
}
