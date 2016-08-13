/************************************************************************
* shift performs a range shift using the shift property of the FFT      *
************************************************************************/
/************************************************************************
* Creator: David T. SandwellScripps Institution of Oceanography)	*
* Date   : 08/09/06							*
************************************************************************/
/************************************************************************
* Modification History							*
* 									*
* Date									*
************************************************************************/ 

#include "image_sio.h"
#include "siocomplex.h"

void shift(ranfft,data,shift)
int ranfft;
fcomplex *data;
double shift;
{
        float arg;
	int i, dir, n2;
        fcomplex cshift;
        dir = -1;
        n2=ranfft/2;
	cfft1d_(&ranfft,data,&dir);
	for(i=0;i<ranfft;i++){
          arg = -2.*PI*shift*i/ranfft;
          if(i > n2) arg = -2.*PI*shift*(i-ranfft)/ranfft;
          cshift=Cexp(arg);
	  data[i] = Cmul(cshift,data[i]);
	}
        dir = 1;
	cfft1d_(&ranfft,data,&dir);
}
