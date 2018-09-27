/************************************************************************
 * radopp modifies doppler frequencies to make them a function of	*
 * 	range.								*
 ************************************************************************/
/************************************************************************
 * Creator: Evelyn J. Price (Scripps Institution of Oceanography)   	*
 * Date   : 11/18/96							*
 ************************************************************************/
/************************************************************************
 * Modification History: *
 *									*
 * Date: 11/18/96							*
 *									*
 ************************************************************************/
#include <math.h>

void radopp(double *fd, double *fdd, double *fddd, double r, double del) {

	double temp1, temp2, temp3, fac;

	fac = r / del;

	temp1 = *fd - *fdd * (fac) + *fddd * pow(fac, 2.0);
	temp2 = *fdd / del - 2.0 * (*fddd) * (fac) / del;
	temp3 = *fddd / (pow(del, 2.0));

	*fd = temp1;
	*fdd = temp2;
	*fddd = temp3;
}
