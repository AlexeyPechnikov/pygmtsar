/************************************************************************
 * rng_cmp performs the range compression operation on raw radar echos.  *
 *	using a precomputed range reference function.			*
 ************************************************************************/
/************************************************************************
 * Creator: Evelyn J. Price	(Scripps Institution of Oceanography)	*
 * Date   : 11/18/96							*
 ************************************************************************/
/************************************************************************
 * Modification History							*
 * 									*
 * Date									*
 ************************************************************************/

#include "gmt.h"
#include "siocomplex.h"
#include "soi.h"
#include <math.h>

void rng_cmp(void *API, int ranfft, fcomplex *data, fcomplex *ref) {
	int i;
	GMT_FFT_1D(API, (float *)data, ranfft, GMT_FFT_FWD, GMT_FFT_COMPLEX);
	for (i = 0; i < ranfft; i++) {
		data[i] = Cmul(ref[i], data[i]);
	}
	GMT_FFT_1D(API, (float *)data, ranfft, GMT_FFT_INV, GMT_FFT_COMPLEX);
}
