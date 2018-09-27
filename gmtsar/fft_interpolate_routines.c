/*	$Id: fft_interpolate_routines.c 109 2015-01-19 23:01:24Z sandwell $
 */
#include "gmtsar.h"
/*--------------------------------------------------------------------------------------*/
/* some fft-based interpolation routines - rjm July 2009
 */
/*											*/
/* fft_interpolate_1d(vin, N, vout, factor)	uses fft to interpolate a
 * signal	*/
/*				vin - complex vector length N
 */
/*				vout - complex vector length M = factor*N
 */
/*				factor - interpolation factor
 */
/*				N must be even
 */
/*				adjust amplitude
 */
/*											*/
/* fft_interpolate_2d(min, N, mout, factor)	2D but uses 1D fft
 */
/*				min - complex matrix length NxN
 */
/*				mout - complex matrix length MxM = factor*N
 */
/*											*/
/* fft_arrange_interpolate(in, N, out, factor)	rearrange terms of transformed
 * signal	*/
/*											*/
/*	all vectors and matrices must be pre-allocated
 */
/*	all use the SIO fcomplex and call GMT_FFT_1D or GMT_FFT_2D
 */
/* 	fcomplex - struct {float r; float i} where r is real and i imag
 */
/*--------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
void fft_arrange_interpolate(struct FCOMPLEX *in, int N, struct FCOMPLEX *out, int M, int ifactor) {
	int i, N2;

	N2 = N / 2;

	for (i = 0; i < M; i++)
		out[i].r = out[i].i = 0.0f;

	for (i = 0; i < N2; i++) {
		out[i].r = in[i].r;
		out[i].i = in[i].i;
		out[i + (2 * ifactor - 1) * N2].r = in[i + N2].r;
		out[i + (2 * ifactor - 1) * N2].i = in[i + N2].i;
	}
}
/*------------------------------------------------------------------------*/
void fft_interpolate_1d(void *API, struct FCOMPLEX *in, int N, struct FCOMPLEX *out, int ifactor) {
	int i, M;

	M = ifactor * N;

	/* forward fft */
	// dir = -1;
	// cfft1d_(&N, in, &dir);
	GMT_FFT_1D(API, (float *)in, N, GMT_FFT_FWD, GMT_FFT_COMPLEX);

	/* re-arrange values in fft */
	fft_arrange_interpolate(in, N, out, M, ifactor);

	/* backward fft */
	// dir = 1;
	// cfft1d_(&M, out, &dir);
	GMT_FFT_1D(API, (float *)out, M, GMT_FFT_INV, GMT_FFT_COMPLEX);

	/* scale amplitude */
	for (i = 0; i < M; i++) {
		out[i].r = ((float)ifactor) * out[i].r;
		out[i].i = ((float)ifactor) * out[i].i;
	}
}
/*--------------------------------------------------------------------------------------*/
void fft_interpolate_2d(void *API, struct FCOMPLEX *in, int N1, int M1, struct FCOMPLEX *out, int N, int M, int ifactor) {
	int i, j;
	struct FCOMPLEX *tmp1, *tmp2, *tmp3;
#if 0
	/* sanity checks */
	if (N != (N1 * ifactor)) error_flag = 1;
	if (M != (M1 * ifactor)) error_flag = 1;
#endif
	tmp1 = (struct FCOMPLEX *)malloc(N1 * M * sizeof(struct FCOMPLEX));
	tmp2 = (struct FCOMPLEX *)malloc(N1 * sizeof(struct FCOMPLEX));
	tmp3 = (struct FCOMPLEX *)malloc(N * sizeof(struct FCOMPLEX));

	for (i = 0; i < N1 * M; i++)
		tmp1[i].i = tmp1[i].r = 0.0f;

	if (debug)
		print_complex(in, N1, M1, 0);

	for (i = 0; i < N1; i++)
		fft_interpolate_1d(API, &in[i * M1], M1, &tmp1[i * M], ifactor);

	if (debug)
		print_complex(tmp1, N1, M, 0);

	/* now do columns - need to transpose */
	for (i = 0; i < M; i++) {
		for (j = 0; j < N1; j++)
			tmp2[j] = tmp1[j * M + i];
		fft_interpolate_1d(API, tmp2, N1, tmp3, ifactor);
		for (j = 0; j < N; j++)
			out[j * M + i] = tmp3[j];
	}

	if (debug)
		print_complex(out, N, M, 0);

	free((char *)tmp1);
	free((char *)tmp2);
	free((char *)tmp3);
}
/*--------------------------------------------------------------------------------------*/
