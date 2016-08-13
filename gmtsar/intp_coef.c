/************************************************************************
* intp_coef calculates an 8 pt. sinc function for range migration       *
* 	interpolation.							*
************************************************************************/
/************************************************************************
* Creator: Evelyn J. Price	(Scripps Institution of Oceanography)   *
* Date   : 11/18/96							*
************************************************************************/
/************************************************************************
* Modification History							*
* 									*
* Date									*
************************************************************************/

#include "soi.h"
#include <math.h>

void intp_coef(int nfilter, float *xintp)
{
      
	int i,j;
	float x,y;
/* compute the interpolation factors */
	for(i=0;i<=nfilter;i++){
	  j=i*8;
	  x = ((float) i)/((float) nfilter);
	  y = sin(PI*x)/PI;
	  if((x!=0.0f) && (x!=1.0f)){
            xintp[j  ] = -y/(3.0f+x);
            xintp[j+1] =  y/(2.0f+x);
            xintp[j+2] = -y/(1.0f+x);
            xintp[j+3] =  y/x;
            xintp[j+4] =  y/(1.0f-x);
            xintp[j+5] = -y/(2.0f-x);
            xintp[j+6] =  y/(3.0f-x);
            xintp[j+7] = -y/(4.0f-x);
	  }
	  else if(x == 0.0f){
            xintp[j  ] = 0.0f;
            xintp[j+1] = 0.0f;
            xintp[j+2] = 0.0f;
            xintp[j+3] = 1.0f;
            xintp[j+4] = 0.0f;
            xintp[j+5] = 0.0f;
            xintp[j+6] = 0.0f;
            xintp[j+7] = 0.0f;
	  }
	  else if(x == 1.0f){
            xintp[j  ] = 0.0f;
            xintp[j+1] = 0.0f;
            xintp[j+2] = 0.0f;
            xintp[j+3] = 0.0f;
            xintp[j+4] = 1.0f;
            xintp[j+5] = 0.0f;
            xintp[j+6] = 0.0f;
            xintp[j+7] = 0.0f;
	  }
	}
}
