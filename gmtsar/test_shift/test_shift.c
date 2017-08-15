/*  program to test the fft shift routine */
#include "../include/soi.h"
#include "../include/siocomplex.h"
#include "../include/gmtsar.h"
#include <math.h>
void main(int argc, char *argv[]) {

int k,nd=1024;
fcomplex *datai, *datao;
double arg;

	void    *API = NULL; /* GMT API control structure */

       /* Begin: Initializing new GMT session */
        API = GMT_Create_Session (argv[0], 0U, 0U, NULL);

datai  = (fcomplex *) malloc(nd*sizeof(fcomplex));
datao  = (fcomplex *) malloc(nd*sizeof(fcomplex));

/*  fill the data array with a gaussian function */

  for (k=0; k<nd; k++){
    datai[k].r=0.0;
    datai[k].i=0.0;
    arg=(k-nd/2)*(k-nd/2)/4.;
    datai[k].r=(k-nd/2)*exp(-arg);
    datai[k].i=-(k-nd/2)*exp(-arg);
    datao[k].r=datai[k].r;
    datao[k].i=datai[k].i;
  }
  shift(API,nd,datao,50.5);
  //shift(nd,datao,50.5);
  for (k=0; k<nd; k++){
    printf(" %d %lf %lf %lf %lf \n",k,datai[k].r,datai[k].i,datao[k].r,datao[k].i);
  }
}
