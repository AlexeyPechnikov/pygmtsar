/*	$Id: do_freq_xcorr.c 109 2015-01-19 23:01:24Z sandwell $	*/
/*-------------------------------------------------------*/
#include "gmt.h"
#include "gmtsar.h"
#include "sarleader_fdr.h"
#include "siocomplex.h"
#include "xcorr.h"
#include <math.h>

/*-------------------------------------------------------------------------------*/
void fft_multiply(void *API, int N, int M, struct FCOMPLEX *c1, struct FCOMPLEX *c2, struct FCOMPLEX *c3) {
	int i, j, isign;

	/* do forward fft 					*/
	GMT_FFT_2D(API, (float *)c1, N, M, GMT_FFT_FWD, GMT_FFT_COMPLEX);
	GMT_FFT_2D(API, (float *)c2, N, M, GMT_FFT_FWD, GMT_FFT_COMPLEX);

	/* multiply a with conj(b)				*/
	/* the isign should shift the results appropriately	*/
	/* normalize by absolute values				*/
	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			isign = (int)lrint(pow(-1.0, (double)(i + j)));
			c3[i * N + j] = RCmul(isign, Cmul((c1[i * N + j]), (Conjg(c2[i * N + j]))));
			/*
			                        if (norm_flag) {
			                                c3[i*N + j] = RCmul(1.0/Cabs(c3[i*N +
			   j]),c3[i*N + j]); if isnan(c3[i*N + j].r) c3[i*N + j].r = 0.0; if
			   isnan(c3[i*N + j].i) c3[i*N + j].i = 0.0;
			                                }
			*/
		}
	}

	/* inverse fft  for cross-correlation (c matrix) 	*/
	GMT_FFT_2D(API, (float *)c3, N, M, GMT_FFT_INV, GMT_FFT_COMPLEX);
}
/*-------------------------------------------------------------------------------*/
void do_freq_corr(void *API, struct xcorr *xc, int iloc) {
	int i, j, ii;
	float ipeak, jpeak;
	double cmax, cave, max_corr;

	max_corr = -999.0;
	jpeak = ipeak = -999.0f;
	cmax = -1;
	cave = 0.0;

	/* d1, c1 is the master					*/
	/* d2, c2 is the aligned					*/
	if (debug)
		print_complex(xc->c1, xc->npy, xc->npx, 1);
	if (debug)
		print_complex(xc->c2, xc->npy, xc->npx, 1);

	/* multiply c1 and c2 uisng fft */
	fft_multiply(API, xc->npx, xc->npy, xc->c1, xc->c2, xc->c3);

	/* transfer results into correlation matrix		*/
	for (i = 0; i < xc->nyc; i++) {
		for (j = 0; j < xc->nxc; j++) {

			ii = (i + xc->ny_corr / 2) * xc->npx + j + (xc->nx_corr / 2);
			xc->corr[i * xc->nxc + j] = Cabs(xc->c3[ii]);
			cave += xc->corr[i * xc->nxc + j];

			if (xc->corr[i * xc->nxc + j] > cmax) {
				cmax = xc->corr[i * xc->nxc + j];
				jpeak = j - xc->nxc / 2.0f;
				ipeak = i - xc->nyc / 2.0f;
			}
		}
	}

	if ((ipeak == -999.0) || (jpeak == -999.0)) {
		fprintf(stderr, "error! jpeak %f ipeak %f cmax %lf xc %f \n", jpeak, ipeak, cmax, xc->corr[100]);
		exit(1);
	}

	/* calculate normalized correlation at best point */
	cave /= (xc->nxc * xc->nyc);

	/* estimate maximum correlation using frequency - poor */
	max_corr = (cmax / cave);

	/* use normalized time correlation rather than frequency estimate for value */
	if (ipeak != -999.0)
		max_corr = calc_time_corr(xc, (int)ipeak, (int)jpeak);

	/* put values into correlation matrix and scale to max correlation */
	for (i = 0; i < xc->nxc * xc->nyc; i++)
		xc->corr[i] = max_corr * xc->corr[i] / cmax;

	if (debug)
		fprintf(stderr, " (freq) jpeak %f xoffset %d corr %4.2lf \n", jpeak, xc->x_offset, max_corr);
	if (debug)
		fprintf(stderr, " (freq) ipeak %f yoffset %d corr %4.2lf \n", ipeak, xc->y_offset, max_corr);

	xc->loc[iloc].xoff = -1 * jpeak;
	xc->loc[iloc].yoff = -1 * ipeak;
	xc->loc[iloc].corr = (float)max_corr;
}
