/*	$Id: do_time_int_xcorr.c 39 2013-04-07 00:49:34Z pwessel $	*/
/*-------------------------------------------------------*/
#include <math.h>
#include "gmtsar.h"
#include "xcorr.h"
/*-------------------------------------------------------------------------------*/
double calc_time_corr(struct xcorr xc, int ioff, int joff)
{
	int	ip, jp;
	long long	a, b;
	long long	gamma_num, gamma_denom1, gamma_denom2;
	double	gamma, gamma_denom;

	gamma_denom = gamma_num = 0.0;
	gamma_num = gamma_denom1 = gamma_denom2 = 0;

	/* calculate normalized correlation 	*/
	/* template (b) stays the same		*/
	/* a is master				*/

	for (ip=0; ip<xc.ny_corr; ip++){
		for (jp=0; jp<xc.nx_corr; jp++){

			/* pixel values */
			a = xc.i1[(xc.ysearch+ip+ioff)*xc.npx + jp + joff + xc.xsearch];
			b = xc.i2[(xc.ysearch+ip)*xc.npx + jp + xc.xsearch]; 	

			/* standard correlation */
			gamma_num += (a * b);
			gamma_denom1 += a * a;
			gamma_denom2 += b * b;
			}
		}

	gamma_denom = sqrt(1.0 * gamma_denom1 * gamma_denom2);

	if (gamma_denom == 0.0) {
		if (verbose) fprintf(stderr,"calc_corr: denominator = zero: setting corr to 0 \n");
		gamma = 0.0;
		} else {  gamma = 100.0 * fabs(gamma_num / gamma_denom);}

	if (debug) fprintf(stdout," corr %6.2lf \n", gamma);

	return(gamma);
}
/*-------------------------------------------------------------------------------*/
double calc_time_corr_hat(struct xcorr xc, int ioff, int joff)
{
	int	ip, jp;
	long long	a, b;
	long long	gamma_num, gamma_denom1, gamma_denom2;
	double		gamma_denom;
	double		gamma;

	gamma_num = gamma_denom1 = gamma_denom2 = 0;
	gamma_denom = gamma = 0.0;
	
	for (ip=0; ip<xc.ny_corr; ip++){
		for (jp=0; jp<xc.nx_corr; jp++){

			/* pixel values */
			a = xc.i1[(xc.ysearch+ip+ioff)*xc.npx + jp + joff + xc.xsearch];
			b = xc.i2[(xc.ysearch+ip)*xc.npx + jp + xc.xsearch]; 	

			/* frequency independent */
			gamma_num += ((a * a * b * b));
			gamma_denom1 += ((a * a * a * a));
			gamma_denom2 += ((b * b * b * b));
			}
		}

	gamma_denom = sqrtl(1.0 * gamma_denom1 * gamma_denom2);

	if (gamma_denom == 0.0) {
		fprintf(stderr,"calc_corr: division by zero \n");
		gamma = 0.0;
		} else { gamma = fabs(gamma_num / gamma_denom);} 

	if ( gamma <= 0.5) {
	 	gamma = 0.0;
	} else {
		gamma = 100.0 * sqrt( (gamma * 2.0) - 1.0);
		}

	if (debug) fprintf(stdout," corr %lf \n", gamma);

	return(gamma);
}
/*-------------------------------------------------------------------------------*/
void do_time_corr(struct xcorr xc, int iloc)
{
	int	ioff, joff; 
	int	ic, jc;
	float	ipeak, jpeak;
	float	max_corr;

	/* set parameters */
	max_corr = -1;
	ipeak = jpeak = -9999;

	/* loops to calculate correlation at various offsets 	*/
	/* correlation window may not be the same as offset	*/
	/* ioff, joff specifies the offset			*/
	/* ic, jc is offset as non-negative for matrix		*/

	for (ioff=-xc.ysearch; ioff<xc.ysearch; ioff++){
		for (joff=-xc.xsearch; joff<xc.xsearch; joff++){

			ic = ioff + xc.ysearch;
			jc = joff + xc.xsearch;

			/* calculate the correlation for each patch 	*/
			/* 0 means standard correlation			*/
			if (xc.corr_flag == 0) xc.corr[ic*xc.nxc + jc] = calc_time_corr(xc, ioff, joff); 
			/* 1 means Gatelli correlation			*/
			if (xc.corr_flag == 1) xc.corr[ic*xc.nxc + jc] = calc_time_corr_hat(xc, ioff, joff); 

			if (fabs(xc.corr[ic*xc.nxc + jc]) > max_corr){
				max_corr = (float)fabs(xc.corr[ic*xc.nxc + jc]);
				jpeak = joff;	
				ipeak = ioff;	
				}
			}
		}

	/* do not do further interpolation in time domain */
	if (debug) fprintf(stdout," (time) jpeak %f xoffset %d \n", jpeak, xc.x_offset);
	if (debug) fprintf(stdout," (time) ipeak %f yoffset %d \n", ipeak, xc.y_offset);

	xc.loc[iloc].xoff = -1 * jpeak;
	xc.loc[iloc].yoff = -1 * ipeak;
	xc.loc[iloc].corr = max_corr;
}
