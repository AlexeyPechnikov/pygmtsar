/************************************************************************
* cfft1d is a subroutine used to call and initialize perflib Fortran FFT *
* routines in a Sun computer.                                           *
************************************************************************/
/************************************************************************
* Creator: David T. Sandwell	(Scripps Institution of Oceanography    *
* Date   : 12/27/96                                                     *
************************************************************************/
/************************************************************************
* Modification history:                                                 *
*  04/01/98  - changed to have arguments be pointers (Fotran callable)  *
*                                                                       *
* DATE                                                                  *
************************************************************************/
 
#include "../include/soi.h"
#include <malloc.h>

cfft1d_(np,c,dir)
int *np,*dir;
fcomplex *c;
{

	static float *work;
	static int nold = 0;
	int i,n;

/* Initialize work array with sines and cosines to save CPU time later 
   This is done when the length of the FFT has changed or when *dir == 0. */

        n = *np;
	if((n != nold) || (*dir == 0)){
	  if(nold != 0) free((char *) work);
	  if((work = (float *) malloc((4*n+30)*sizeof(float))) == NULL){
	    fprintf(stderr,"Sorry, can't allocate mem.\n");
	    return(-1);
	  }
	  cffti_(np,work);
	  nold = n;
	}

/* Do forward transform with NO normalization.  Forward is exp(+i*k*x) */

	if(*dir == -1){
	  cfftf_(np,c,work);
	}

/* Do inverse transform with normalization.  Inverse is exp(-i*k*x) */

	if(*dir == 1){
	  cfftb_(np,c,work);
          for(i=0;i<n;i++){
	    c[i].r = c[i].r/((float) n);
	    c[i].i = c[i].i/((float) n);
	  }
	}

}
