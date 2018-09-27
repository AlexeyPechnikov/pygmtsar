/************************************************************************
 * rng_filter applies a low-pass filter in the fourier domain by         *
 *	zeroing in wavenumber space                                     *
 ************************************************************************/
/************************************************************************
 * Creator: David T. Sandwell (Scripps Institution of Oceanography)	*
 * Date   : 10/09/14							*
 ************************************************************************/
/************************************************************************
 * Modification History							*
 * 01/25/18 DTS modified to use the GMT fft routines			*
 * Date									*
 ************************************************************************/

#include "gmt.h"
#include "siocomplex.h"
#include "soi.h"
#include <math.h>

void rng_filter(void *API, fcomplex *cin, int nffti, fcomplex *cout) {
	int i, nf, nt;
	nf = .70 * nffti / 2;
	nt = nffti - 1;

	/* do the forward fft */

	// dir = -1;
	// cfft1d_(&nffti,cin,&dir);
	GMT_FFT_1D(API, (float *)cin, nffti, GMT_FFT_FWD, GMT_FFT_COMPLEX);

	/* first zero the output array */

	for (i = 0; i < nffti; i++) {
		cout[i].r = 0.;
		cout[i].r = 0.;
	}

	/* only keep the lower frequencies */

	for (i = 0; i < nf; i++) {
		cout[i].r = cin[i].r;
		cout[i].i = cin[i].i;
		cout[nt - i].r = cin[nt - i].r;
		cout[nt - i].i = cin[nt - i].i;
	}

	/* now inverse fft */

	// dir = 1;
	// cfft1d_(&nffti,cout,&dir);
	GMT_FFT_1D(API, (float *)cout, nffti, GMT_FFT_INV, GMT_FFT_COMPLEX);
}
