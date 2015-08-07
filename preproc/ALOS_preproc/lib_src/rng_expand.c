/************************************************************************
* rng_expand interpolates a complex array to twice the length by        *
*	zero padding in the wavenumber space                            *
************************************************************************/
/************************************************************************
* Creator: David T. Sandwell (Scripps Institution of Oceanography)	*
* Date   : 06/21/07							*
************************************************************************/
/************************************************************************
* Modification History							*
* 									*
* Date									*
************************************************************************/ 

#include"image_sio.h"
#include"siocomplex.h"

void cfft1d_(int *, fcomplex *, int *);

void rng_expand(cin,nffti,cout,nffto)
int nffti, nffto;
fcomplex *cin, *cout;
{
	int i, dir, n2;
        n2 = nffti/2;

/* do the forward fft */

        dir = -1;
	cfft1d_(&nffti,cin,&dir);

/* first zero the output array */

	for(i=0;i<nffto;i++){
		cout[i].r=0.;
		cout[i].r=0.;
	}

/* then move the input to the output 1 to 1 and 2 to 4 
   multiply by 2 because all the zeros were added */
		
	for(i=0;i<n2;i++){
		cout[i].r=2.*cin[i].r;
		cout[i].i=2.*cin[i].i;
		cout[i+3*n2].r=2.*cin[i+n2].r;
		cout[i+3*n2].i=2.*cin[i+n2].i;
	}

/* now inverse fft */
	
        dir = 1;
	cfft1d_(&nffto,cout,&dir);
}
