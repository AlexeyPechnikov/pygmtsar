/************************************************************************
 * various utilities needed for gmtsar                                  *
 * modified from readgrd                       		                *
 ************************************************************************/
/************************************************************************
 * Creator: Rob J. Mellors, San Deigo State University                  *
 * Date   : December 18, 2007                                           *
 *                                                                      *
 * Modification history:                                                *
 *   some slight format changes - Dec 18, 2007 - RJM 		        *
 *   comments added - Jan 14, 2010 - DTS                                *
 ************************************************************************/

#include "gmtsar.h"
#include "lib_functions.h"

/*---------------------------------------------------------------*/
/* check endian of machine 	*/
/* 1 if big; -1 if little	*/
int is_big_endian_()
{
	union
	{
	long l;
	char c[sizeof (long) ];
	} u;
	u.l = 1;
	return( u.c[sizeof(long) - 1] ==  1 ? 1 : -1);
}
int is_big_endian__()
{
	return is_big_endian_();
}
/*---------------------------------------------------------------*/
/* write out error message and exit 				*/
/* use two strings to allow more complicated error messages 	*/
void die (char *s1, char *s2)
{
	fprintf(stderr," %s %s \n",s1,s2);
	exit(1);
}
/*---------------------------------------------------------------*/
/* cross3 is a routine to take the cross product of 3-D vectors  */
void cross3(double *a, double *b, double *c)

/* input and output vectors  having 3 elements */ 

{

       c[0] =  (a[1]*b[2]) - (a[2]*b[1]);
       c[1] = (-a[0]*b[2]) + (a[2]*b[0]);
       c[2] = (a[0]*b[1]) - (a[1]*b[0]);

}
/*----------------------------------------------------------------------*/
/* find_length returns length of a 3 component vector 	*/
double find_length(double *x)
{
double	length;

	length = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

	return(length);
}
/*----------------------------------------------------------------------*/
/* find_distance returns distance between two 3 component vector 	*/
double find_distance(double *xs, double *x)
{
int	i;
double	dx[3];

	for (i=0; i<3; i++) dx[i] = xs[i] - x[i];

	return(find_length(dx));
}
/*----------------------------------------------------------------------*/
/* find_distance returns distance between two 3 component vector 	*/
/* sometimes cannot specify as vectors easily				*/
double find_distance3(double x11, double x12, double x13, double x21, double x22, double x23)
{
double	dx[3], distance;

	dx[0] = x11 - x21;
	dx[1] = x12 - x22;
	dx[2] = x13 - x23;

	distance = find_length(dx);

	return(distance);
}
/*---------------------------------------------------------------------------*/
void find_unit_vector(double *x, double *xn)
{
int	i;
double	r;

	r = find_length(x);
	for (i=0; i<3; i++) xn[i] = x[i]/r;
} 
/*---------------------------------------------------------------------------*/
/* find seconds		*/
void	get_seconds(struct PRM p, double *start, double *end)
{
int	m;
double	doy;
double	n_secs_day;

n_secs_day = 24.0*60.0*60.0;
doy = p.clock_start;
m = p.nrows - p.num_valid_az;

/* also correct the time for the azimuth shift */
if(p.sub_int_a < 0.) p.sub_int_a = 0.;
*start = n_secs_day*doy + (p.ashift+p.sub_int_a)/(p.prf) + (1.0*m)/(2.0*p.prf);
*end = *start + p.num_patches * p.num_valid_az/p.prf;

}
/*---------------------------------------------------------------------------*/
int	find_fft_length(int n)
{
int	nfft;

	nfft = 2;

	while (nfft < n) nfft = 2*nfft;

	if (debug) fprintf(stderr,"find_fft_length:\n...data length n %d nfft %d \n\n",n,nfft);

	return(nfft);
}
/*-----------------------------------------------------------------------*/
int	geo2latlon(double *x, double *llt, struct PRM r)
{
	double fll;
	fll = (r.ra - r.rc) / r.ra;

	xyz2plh(x, llt, r.ra, fll);	

	if (verbose) fprintf(stderr, " geo2latlon: %8.3f %8.3f %8.0f \n", x[0], x[1], x[2]);

	return(EXIT_SUCCESS);
}
/*-----------------------------------------------------------------------*/
