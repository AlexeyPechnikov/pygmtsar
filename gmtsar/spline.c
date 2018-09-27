#include <stdio.h>
#include <stdlib.h>

/* February 11, 2010 converted from fortran routines spline.f providrd by Bob
   Parker the conversion was first performed using f2c.  Then Rob Mellors made
   some changes to eliminate the f2c.h file so the code has no dependencies.
   David Sandwell tested the cod against the fortran version and the differences
   were al the level of the machine precision except when the interpolated point
   is outside the bounds of the unput data.  In this case the fortran code and
   C-code diverge.  THis is not yet understood. */

/*  comments from original fortran code follow
   Finds array  s  for spline interpolator  eval.
   Based on an alogorithm by Richard F. Thompson in 'Spline
   interpolation on a digital computer' (x-692-70-261) 1970 Goddard
   Space Flight Center.  Note these splines are not 'natural' since the
   2nd derivative is not zero x(1) and x(n).
   nn  number of data points supplied (may be negative, see below)
   x   array of x-coordinates where function is sampled.  The sequence
       xx(1),xx(2)...  must be a strictly increasing sequence.
       THIS IS NOT CHECKED by spline.
   u   array of sample values that are to be interpolated.
   s   output array of 2nd derivative at sample points.
   a   working space array of dimension at least  nn.

   If the user wishes derivatives at the ends of series to assume
   specified values, du(1)/dx, du(nn)/dx  must be placed in  s1, s2
   and in the call nn=-number of terms in series.  Normally a parabola
   is fitted through the 1st and last 3 points to find the slopes.

   If 2 points are supplied, a straight lines is fitted.
   If 3 points are supplied, a parabola is fitted.

   At evaluation time, straight-line extrapolation is provided
   outside the interval (x(1), x(n)).    */

/*------------------------------------------------------------*/
int spline_(int *istart, int *nn, double *x, double *u, double *s, double *a) {
	int n, n1, i1, i, j;
	double r1, r2, r3, r4, r5, r6;
	double q1, qn, c;

	/* decrease arrays so that the fortran indexing works */
	--a;
	--s;
	--u;
	--x;

	*istart = 0;
	n = abs(*nn);

	if (n <= 2) {
		s[1] = 0.0;
		s[2] = 0.0;
		return (0);
	}

	r1 = u[2] - u[1];
	r2 = x[2] - x[1];
	r3 = u[3] - u[1];
	r4 = x[3] - x[1];
	r5 = r2;
	r6 = r4;
	q1 = (r1 / (r5 * r5) - r3 / (r6 * r6)) / (1.0 / r2 - 1.0 / r4);

	r1 = u[n - 1] - u[n];
	r2 = x[n - 1] - x[n];
	r3 = u[n - 2] - u[n];
	r4 = x[n - 2] - x[n];
	r5 = r2;
	r6 = r4;
	qn = (r1 / (r5 * r5) - r3 / (r6 * r6)) / (1.0 / r2 - 1.0 / r4);

	if (n == 3) {
		s[1] = (qn - q1) / (x[3] - x[1]);
		s[2] = s[1];
		s[3] = s[1];
		return (0);
	}

	if (*nn <= 0) {
		q1 = s[1];
		qn = s[2];
	}

	s[1] = ((u[2] - u[1]) / (x[2] - x[1]) - q1) * 6.0;
	n1 = n - 1;
	i1 = n1;

	for (i = 2; i <= i1; ++i) {
		s[i] = (u[i - 1] / (x[i] - x[i - 1]) - u[i] * (1.0 / (x[i] - x[i - 1]) + 1.0 / (x[i + 1] - x[i])) +
		        u[i + 1] / (x[i + 1] - x[i])) *
		       6.0;
	}

	s[n] = (qn + (u[n1] - u[n]) / (x[n] - x[n1])) * 6.0;
	a[1] = (x[2] - x[1]) * 2.f;
	a[2] = (x[2] - x[1]) * 1.5 + (x[3] - x[2]) * 2.0;
	s[2] -= s[1] * 0.5;
	i1 = n1;

	for (i = 3; i <= i1; ++i) {
		c = (x[i] - x[i - 1]) / a[i - 1];
		a[i] = (x[i + 1] - x[i - 1]) * 2.f - c * (x[i] - x[i - 1]);
		s[i] -= c * s[i - 1];
	}

	c = (x[n] - x[n1]) / a[n1];
	a[n] = (2.f - c) * (x[n] - x[n1]);
	s[n] -= c * s[n1];

	s[n] /= a[n];
	i1 = n1;

	for (j = 1; j <= i1; ++j) {
		i = n - j;
		s[i] = (s[i] - (x[i + 1] - x[i]) * s[i + 1]) / a[i];
	}

	return (EXIT_SUCCESS);
}
/*------------------------------------------------------------*/
/*  comments from original fortran code follow
   Performs cubic spline interpolation of a function sampled unequally
   in  x.  The routine  spline  must be called once before  eval
   is used, to set up the array s.  See comments in  spline.
   y   the coordinate at which function value is desired.
   nn  number of samples of original function.
   x   array of sample coordinates.  The sequence x(1),x(2).....x(nn)
       must be strictly increasing.  THIS IS NOT CHECKED!
   u   array of samples of function at the coordinates x(1),x(2)...
   s   array of  2nd derivatives at sample points.  Found by  spline
       which must be called once before beginning interpolation.
   If  y  falls outside  (x(1), x(nn))  extrapolation is based upon
   the tangent straight line at the endpoint.a   */

int evals_(int *istart, double *y, int *nn, double *x, double *u, double *s, double *eval) {
	int n, l0 = 0, l1, i1;
	double h, r1, h6, xi;

	--s;
	--u;
	--x;

	n = abs(*nn);
	if (*istart > n || *istart < 1) {
		*istart = (n - 1) * (*y - x[1]) / (x[n] - x[1]) + 1;
	}

	if (*y <= x[1]) {
		goto L3000;
	}

	if (*y >= x[n]) {
		goto L3100;
	}

	if (*y < x[*istart]) {
		goto L1200;
	}

	i1 = n;
	for (l1 = *istart; l1 <= i1; ++l1) {
		if (x[l1] > *y) {
			goto L1150;
		}
	}

L1150:
	l0 = l1 - 1;
	goto L1500;

L1200:
	i1 = *istart;
	for (l1 = 1; l1 <= i1; ++l1) {
		l0 = *istart - l1;
		if (x[l0] <= *y) {
			goto L1350;
		}
	}

L1350:
	l1 = l0 + 1;
L1500:
	*istart = l0;
	xi = (*y - x[l0]) / (x[l1] - x[l0]);

	r1 = x[l1] - x[l0];
	h6 = r1 * r1 / 6.f;
	*eval = u[l0] + xi * (u[l1] - u[l0] - h6 * (1.0 - xi) * (s[l1] * (xi + 1.f) + s[l0] * (2.0 - xi)));
	return 0;

L3000:
	h = x[2] - x[1];
	*eval = u[1] + (*y - x[1]) * ((u[2] - u[1]) / h - h * (s[2] - s[1] * 2.0) / 6.0);
	return 0;
L3100:
	h = x[n] - x[n - 1];
	*eval = u[n] + (*y - x[n]) * ((u[n] - u[n - 1]) / h + h * (s[n] * 2.0 - s[n - 1]) / 6.0);

	return (EXIT_SUCCESS);
}
/*------------------------------------------------------------*/
