/*	$Id: print_results.c 33 2013-04-06 05:37:15Z pwessel $	*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gmtsar.h"
#include "siocomplex.h"
#include "xcorr.h"
/*-------------------------------------------------------------------------------*/
void print_results(struct xcorr *xc, int iloc)
{ 
int	ishft;
int	interp;
float	xoff, xfrac;
float	yoff, yfrac;
float	corr;
	
	xoff = xc->loc[iloc].xoff;
	xfrac = xc->loc[iloc].xfrac;
	yoff = xc->loc[iloc].yoff;
	yfrac = xc->loc[iloc].yfrac;
	corr = xc->loc[iloc].corr;

	if (xc->interp_flag) {
		interp = xc->interp_factor;
		} else {
		interp = 1;
		}

	ishft = (int) xc->loc[iloc].y * xc->astretcha;
	
	/* account for range interpolation (xc->ri) */
	if (debug) fprintf(stdout, " xoff %f xfrac %f xoff/ri %f rshift %d yoff %f yfrac %f ashift %d\n", xoff, xfrac, xoff / (float) xc->ri, xc->x_offset, yoff, yfrac, xc->y_offset);

	xoff = (xoff / (float) xc->ri) - (xfrac / (float) xc->ri) + xc->x_offset;
	yoff = yoff - yfrac + xc->y_offset + ishft;

	if (verbose) {
		fprintf(stdout, " location %d (%3d,%3d) interpolation (range %d corr %d) correlation %6.2f offset (%6.3f,%6.3f) \n"
		,iloc, xc->loc[iloc].x, xc->loc[iloc].y, xc->ri, interp, corr, xoff, yoff);
		}

	fprintf(xc->file," %d %6.3f %d %6.3f %6.2f \n",xc->loc[iloc].x,xoff,xc->loc[iloc].y,yoff,corr);
}
/*-------------------------------------------------------------------------------*/
void print_complex(struct FCOMPLEX *a, int ny, int nx, int real_flag)
{
int	i, j;

	if (real_flag == 0) fprintf(stdout,"\ncomplex: \n");
	if (real_flag == 1) fprintf(stdout,"\ncomplex (real only): \n");

	for (i=0; i<ny; i++) {
		for (j=0; j<nx; j++) {
			if (real_flag == 0) fprintf(stdout,"(%6.2f,%6.2f) ", a[i*nx+j].r, a[i*nx+j].i);
			if (real_flag == 1) fprintf(stdout,"%3.1f ", a[i*nx+j].r);
			}
		fprintf(stdout,"\n");
		}
	fprintf(stdout,"\n");
}
/*-------------------------------------------------------------------------------*/
void print_float(float *a, int ny, int nx)
{
int	i, j;

	for (i=0; i<ny; i++) {
		for (j=0; j<nx; j++) fprintf(stdout, " %4.2f ", a[i*nx+j]); 
		fprintf(stdout,"\n");
		}
}
/*-------------------------------------------------------------------------------*/
void print_double(double *a, int ny, int nx)
{
int	i, j;

	for (i=0; i<ny; i++) {
		for (j=0; j<nx; j++) fprintf(stdout, " %4.2lf ", a[i*nx+j]); 
		fprintf(stdout,"\n");
		}
}
/*-------------------------------------------------------------------------------*/
void print_int(int *a, int ny, int nx)
{
int	i, j;

	for (i=0; i<ny; i++) {
		for (j=0; j<nx; j++) fprintf(stdout, " %d ", a[i*nx+j]); 
		fprintf(stdout,"\n");
		}
}
/*-------------------------------------------------------------------------------*/
