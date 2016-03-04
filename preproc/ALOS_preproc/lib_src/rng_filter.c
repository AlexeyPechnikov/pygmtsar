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
* 									*
* Date									*
************************************************************************/ 

#include"image_sio.h"
#include"siocomplex.h"

void cfft1d_(int *, fcomplex *, int *);

void rng_filter(cin,nffti,cout)
int nffti;
fcomplex *cin, *cout;
{
	int i, dir, nf, nt;
        nf = .70*nffti/2;
        nt = nffti-1;

/* do the forward fft */

        dir = -1;
	cfft1d_(&nffti,cin,&dir);

/* first zero the output array */

	for(i=0;i<nffti;i++){
		cout[i].r=0.;
		cout[i].r=0.;
	}

/* only keep the lower frequencies */
		
	for(i=0;i<nf;i++){
		cout[i].r=cin[i].r;
		cout[i].i=cin[i].i;
		cout[nt-i].r=cin[nt-i].r;
		cout[nt-i].i=cin[nt-i].i;
	}

/* now inverse fft */
	
        dir = 1;
	cfft1d_(&nffti,cout,&dir);
}
