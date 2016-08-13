/************************************************************************
* cfft1d is a subroutine to do a 1-D fft using FFTW routines            *
************************************************************************/
/************************************************************************
* Creator: David T. Sandwell    (Scripps Institution of Oceanography    *
* Date   : 12/27/96                                                     *
************************************************************************/
/************************************************************************
* Modification history:                                                 *
*  04/01/98  - changed to have arguments be pointers (Fotran callable)  *
*  10/16/03  - changed to call fftw instead of the sun perflib          *
************************************************************************/

#include "/sw/include/fftw3.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

cfft1d_(np,c,dir)
int *np,*dir;
fftwf_complex *c;
{

        static int nold = 0;
        static fftwf_plan pf, pi;
        int i,n;

/* Make the plans for FFTW and destroy the old ones if they exist.
   This is done when the length of the FFT has changed or when *dir == 0. */

        n = *np;
        if((n != nold) || (*dir == 0)){
          if(nold != 0) {
          fftwf_destroy_plan(pf);
          fftwf_destroy_plan(pi);
          }
          pf = fftwf_plan_dft_1d(n,c,c,-1,FFTW_MEASURE);
          pi = fftwf_plan_dft_1d(n,c,c, 1,FFTW_MEASURE);
	  printf(" reset plan \n");
          nold = n;
        }

/* Do forward transform with NO normalization. */

        if(*dir == -1){
          fftwf_execute(pf);
        }

/* Do inverse transform with normalization. */

        if(*dir == 1){
          fftwf_execute(pi);
          for(i=0;i<n;i++){
            c[i][0] = c[i][0]/((float) n);
            c[i][1] = c[i][1]/((float) n);
          }
        }
}
