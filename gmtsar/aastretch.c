/************************************************************************
 * aastretch performs an azimuth dependent azimuth interpolation         *
 ************************************************************************/
/************************************************************************
 * Creator: Meng Wei, Scripps Institution of Oceanography)     		*
 * Date   : 09/29/06                                                     *
 ************************************************************************/
/************************************************************************
 * Modification History                                                  *
 * 02/23/11  DTS  modified interpolation to use cubic spline             *
 * Date                                                                  *
 ************************************************************************/

#include "gmt.h"
#include "gmtsar.h"
#include "siocomplex.h"
#include "soi.h"
#include <math.h>

void aastretch(fcomplex **fdata, int ipatch, int nrows, int num_valid_az, int num_rng_bins, float coef) {
	int i, j, low_ind, istart, nsmx;
	double *xs, *real, *imag, *ss, *as;
	double yss, test;

	low_ind = (nrows - num_valid_az) / 2;
	nsmx = low_ind / 2;
	/* check to see if the total shift will be larger than 1/4 the length of the
	 * aperture */
	if (coef * num_valid_az >= nsmx || coef * num_valid_az <= -nsmx) {
		fprintf(stderr, "The a_stretch_a is too large.\n");
		exit(-1);
	}

	/* allocate memory for vectors */
	xs = (double *)malloc((nrows) * sizeof(double));
	ss = (double *)malloc((nrows) * sizeof(double));
	as = (double *)malloc((nrows) * sizeof(double));
	real = (double *)malloc((nrows) * sizeof(double));
	imag = (double *)malloc((nrows) * sizeof(double));

	/* calculate the original and shifted positions */
	for (i = 0; i < nrows; i++) {
		xs[i] = (ipatch - 1) * num_valid_az + i - low_ind;
	}

	/* process each column */
	for (j = 0; j < num_rng_bins; j++) {

		/* get a column of data from the 2-D array */
		for (i = 0; i < nrows; i++) {
			real[i] = fdata[i][j].r;
			imag[i] = fdata[i][j].i;
		}

		/* interpolate each column - real then imaginary */
		spline_(&istart, &nrows, xs, real, ss, as);
		for (i = low_ind - nsmx; i < nrows - nsmx; i++) {
			yss = xs[i] * (1 + coef);
			evals_(&istart, &yss, &nrows, xs, real, ss, &test);
			fdata[i][j].r = (float)test;
		}
		spline_(&istart, &nrows, xs, imag, ss, as);
		for (i = low_ind - nsmx; i < nrows - nsmx; i++) {
			yss = xs[i] * (1 + coef);
			evals_(&istart, &yss, &nrows, xs, imag, ss, &test);
			fdata[i][j].i = (float)test;
		}
	}
	free((double *)xs);
	free((double *)real);
	free((double *)imag);
	free((double *)ss);
	free((double *)as);
}
