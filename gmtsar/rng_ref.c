/************************************************************************
 * rng_ref computes the range reference function used for range * compression of
 *raw radar echoes					*
 ************************************************************************/
/************************************************************************
 * Creator: Evelyn J. Price 	(Scripps Institution of Oceanography)	*
 * Date   : 11/18/96							*
 ************************************************************************/
/************************************************************************
 * Modification History							*
 *									*
 * Date									*
 * 08/06/08 DTS - secondary range migration code added  from HZ          *
 *									*
 ************************************************************************/

#include "gmt.h"
#include "siocomplex.h"
#include "soi.h"
#define RW 0666

void rng_ref(void *API, int ranfft, float delr, fcomplex *ref1) {
	float ts, t, phase, rhww1, win1;
	float wgt;
	fcomplex cref;
	fcomplex junk;
	int npts, k1start, k1end, icaltone;
	int i, k, index = 0;
	double delw, rc, kx, omega, omega0;

	/* Compute range reference function */

	/* Compute number of points in reference function. */
	if ((strcmp(off_vid, "y") == 0) || (strcmp(off_vid, "Y") == 0)) {
		npts = 2.0 * fs * pulsedur;
		ts = 1.0 / (2.0 * fs);
	}
	else {
		npts = fs * pulsedur;
		ts = 1. / fs;
	}

	if (fmod(npts, 2.0) == 0.0)
		npts = npts + 1;

	k1start = (int)fabs(pctbw) * npts;
	k1end = npts - (int)fabs(pctbw) * npts;

	/* compute reference function */

	k = 0;
	for (i = -npts / 2; i <= npts / 2; i++) {
		t = i * ts;
		if ((strcmp(off_vid, "y") == 0) || (strcmp(off_vid, "Y") == 0)) {
			phase = PI * slope * t * t + PI * fs * t;
			if ((k >= k1start) && (k < k1end)) {
				ref1[i + npts / 2].r = cosf(phase);
				ref1[i + npts / 2].i = 0.0f;
			}
			else {
				ref1[i + npts / 2].r = 0.0f;
				ref1[i + npts / 2].i = 0.0f;
			}
		}
		else {
			phase = PI * slope * t * t;
			if ((k >= k1start) && (k < k1end)) {
				index = i + npts / 2;
				ref1[index] = Cexp(phase);
			}
			else {
				ref1[index].r = 0.0f;
				ref1[index].i = 0.0f;
			}
		}
		k = k + 1;
	}

	/* window reference function */

	rhww1 = 1.0 - rhww;
	for (i = 0; i < ranfft; i++) {
		if (i < npts) {
			win1 = rhww - rhww1 * cos((2.0 * PI * ((float)i)) / ((float)(npts - 1)));
			junk = ref1[i];
			ref1[i] = RCmul(win1, junk);
		}
		else {
			ref1[i].r = 0.0f;
			ref1[i].i = 0.0f;
		}
	}

	/*  extend the chirp */

	if (nextend > 0) {
		for (i = 0; i < npts; i++) {
			if ((k = i - nextend) < 0)
				k = k + ranfft;
			ref1[k] = ref1[i];
		}
		for (i = 0; i < nextend; i++) {
			if ((k = npts - nextend + i) < 0)
				k = k + ranfft;
			ref1[k].r = 0.0f;
			ref1[k].i = 0.0f;
		}
	}

	/* Calculate fft of range reference function */

	// dir = -1;
	// cfft1d_(&ranfft,ref1,&dir);
	GMT_FFT_1D(API, (float *)ref1, ranfft, GMT_FFT_FWD, GMT_FFT_COMPLEX);

	/* zero out dc and caltone location */

	icaltone = (int)(caltone * ranfft + 0.5);

	for (i = 0; i < 6; i++) {
		wgt = 0.5 - 0.5 * cos(((float)i) / 5.0 * PI);
		ref1[i] = RCmul((float)wgt, ref1[i]);
		if ((strcmp(off_vid, "y") == 0) || (strcmp(off_vid, "Y") == 0)) {
			ref1[i + ranfft / 2] = RCmul((float)wgt, ref1[i + ranfft / 2]);
			ref1[ranfft / 2 - i - 1] = RCmul((float)wgt, ref1[ranfft / 2 - 1 - i]);
		}

		ref1[i + icaltone] = RCmul((float)wgt, ref1[i + icaltone]);
		k = icaltone - i;
		if (k < 0)
			k = k + ranfft;
		ref1[k] = RCmul((float)wgt, ref1[k]);

		if ((strcmp(off_vid, "y") == 0) || (strcmp(off_vid, "Y") == 0)) {
			k = i + ranfft - icaltone;
			if (k >= ranfft)
				k = k - ranfft;
			ref1[k] = RCmul((float)wgt, ref1[k]);
			ref1[ranfft - icaltone - 1 - i] = RCmul((float)wgt, ref1[ranfft - icaltone - 1 - i]);
		}

		if (i > 0)
			ref1[ranfft - i] = RCmul((float)wgt, ref1[ranfft - i]);
	}

	/* apply secondary range migration */

	if ((strcmp(srm, "y") == 0) || (strcmp(srm, "Y") == 0)) {
		delw = PI2 * fs / ranfft;
		rc = near_range + delr * num_rng_bins / 2;
		kx = PI2 * (fd1 + fdd1 * rc + pow(fddd1 * rc, 2.0)) / vel1;
		omega0 = PI2 * SOL / lambda;
		/*	   printf("delw= %f rc=%f kx=%f omega0 = %f \n",delw,rc,kx,omega0);*/
		phase = -0.25 * rc * (lambda / PI2) * pow(kx, 2.0) * pow(((omega0 - ranfft / 2. * delw) / omega0), 2.0);
		/*	   printf("phase = %f \n",phase);*/

		/* loop over frequencies and update the reference function */

		for (i = 0; i < ranfft; i++) {
			k = i;
			if (i > ranfft / 2)
				k = i - ranfft;
			omega = omega0 + k * delw;
			phase = -0.25 * rc * (lambda / PI2) * pow(kx, 2.0) * pow(((omega) / omega0), 2.0);
			/*	      printf("i= %d k= %d phase = %f \n",i,k,phase);*/
			junk = Cexp(phase);
			ref1[i] = Cmul(ref1[i], junk);
		}
	}

	for (i = 0; i < ranfft; i++) {
		cref = Conjg(ref1[i]);
		ref1[i] = RCmul((float)(1.0 / ranfft), cref);
	}
}
