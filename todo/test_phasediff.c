#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// from file https://github.com/gmtsar/gmtsar/blob/master/gmtsar/phasediff.c
void die(char *s1, char *s2) {
  fprintf(stderr, "%s: %s\n", s1, s2);
  exit(-1);
}
void calc_drho(int xdim, double *range, double *topo, double avet, double re, double height, double B, double alpha, double Bx,
               double *drho) {
	int k;
	/* EX: changing to long double for better precision */
	long double rho, sint, cost, cosa, sina, b;
	// long double term1,term2,c,c2,ret,ret2;
	long double term1, c, c2, ret, ret2;

	sina = sin(alpha);
	cosa = cos(alpha);
	c = re + height;
	c2 = c * c;
	b = B;
	for (k = 0; k < xdim; k++) {

		/* compute the look angle using equation (C26) in Appendix C */
		rho = range[k];
		ret = re + avet + topo[k];
		ret2 = ret * ret;
		cost = ((rho * rho + c2 - ret2) / (2. * rho * c));
		// thet = acos(cost);
		if (cost >= 1.)
			die("calc_drho", "cost >= 0");
		sint = sqrtl(1. - cost * cost);

		/* compute the range change using equation (c23) in Appendic C */
		// term1 = -B*(sint*cosa-cost*sina);
		// term2 = B*B*(cost*cosa+sint*sina)*(cost*cosa+sint*sina)/(2.*rho);
		// drho[k] = term1 + term2;

		/* New (Eric Lindsey, April 2015): compute the range change using the full
		 * nonlinear equation */
		// term1 = rho*rho + b*b - 2*rho*b*(sint*cosa-cost*sina);
		// term1 = rho*rho + b*b - 2*rho*b*sin(thet-alpha);
		// drho[k] = -rho + sqrtl(term1);

		/* Compute the offset effect from non-parallel orbit */
		term1 = rho * rho + b * b - 2 * rho * b * (sint * cosa - cost * sina) - Bx * Bx;
		// term1 = rho*rho + b*b - 2*rho*b*(sint*cosa-cost*sina);
		drho[k] = -rho + sqrtl(term1);
	}
}

int main() {
    int xdim = 10;
    double range[] = {750000, 760000, 770000, 780000, 790000, 800000, 810000, 820000, 830000, 840000};
    double topo[] = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
    double avet = 1.0;
    double re = 6371000.0;
    double height = 800.0;
    double B = 150.0;
    double alpha = M_PI / 6;
    double Bx = 75.0;

    double drho[10];
    calc_drho(xdim, range, topo, avet, re, height, B, alpha, Bx, drho);

    printf("drho results:\n");
    for (int i = 0; i < xdim; i++) {
        printf("drho[%d] = %f\n", i, drho[i]);
    }

    return 0;
}
