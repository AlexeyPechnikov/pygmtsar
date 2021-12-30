/************************************************************************
 * extend the state vector for LED files from Sentinel1 and Radarsat2   *
 ***********************************************************************/
/************************************************************************
 * Creator: David T. Sandwell & Xiaohua Xu				*
 * Date   : 03/05/2015                                                   *
 ************************************************************************/
/************************************************************************
 * Modification History:                                                 *
 *                                                                       *
 * Date   : 			                                        *
 *                                                                       *
 ************************************************************************/

#include "gmtsar.h"
#include "lib_functions.h"
#include "orbit.h"

void read_orb(FILE *, struct SAT_ORB *);
void write_orb(FILE *, struct SAT_ORB *);
void polyfit(double *, double *, double *, int *, int *);

int main(int argc, char **argv) {
	struct SAT_ORB *orb;
	struct SAT_ORB *orb_out;
	FILE *infile, *outfile;
	int nd_in, nd_add, nd_out;
	int n, k;
	double *pt, *px, *py, *pz, *ppt, *ppx, *ppvx, *ppy, *ppvy, *ppz, *ppvz;
	double frac, idsec;

	void extrapolate_poly(double *, double *, double *, double *, double *, int, int, int);
	int Nodr;

	if ((argc != 4)) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: extend_orbit file_in.LED file_out.LED frac\n");
		fprintf(stderr, "    orig.LED   - original LED file \n");
		fprintf(stderr, "    extend.LED - extended LED file \n");
		fprintf(stderr, "    frac       - fraction of frame length to add to each "
		                "end (e.g. 0.25) \n");
		fprintf(stderr, "\n");
		return (0);
	}

	frac = atof(argv[3]);

	/*  get the orbit data */
	if ((infile = fopen(argv[1], "r")) == NULL)
		die("Can't open ", argv[1]);
	orb = (struct SAT_ORB *)malloc(sizeof(struct SAT_ORB));
	orb_out = (struct SAT_ORB *)malloc(sizeof(struct SAT_ORB));
	read_orb(infile, orb);

	/* make a longer structure for the output LED file */
	nd_in = orb->nd;
	nd_add = nd_in * frac;
	nd_out = nd_in + 2 * nd_add;
	orb_out->nd = nd_out;
	orb_out->iy = orb->iy;
	orb_out->id = orb->id;
	idsec = orb->dsec;
	orb_out->dsec = orb->dsec;
	orb_out->sec = orb->sec - nd_add * idsec;
	orb_out->points = (struct ORB_XYZ *)malloc(orb_out->nd * sizeof(struct ORB_XYZ));

	/* allocate the memory for the vectors to be fit and fill the vectors */
	pt = malloc(nd_in * sizeof(double));
	px = malloc(nd_in * sizeof(double));
	py = malloc(nd_in * sizeof(double));
	pz = malloc(nd_in * sizeof(double));
	ppt = malloc(nd_add * 2 * sizeof(double));
	ppx = malloc(nd_add * 2 * sizeof(double));
	ppvx = malloc(nd_add * 2 * sizeof(double));
	ppy = malloc(nd_add * 2 * sizeof(double));
	ppvy = malloc(nd_add * 2 * sizeof(double));
	ppz = malloc(nd_add * 2 * sizeof(double));
	ppvz = malloc(nd_add * 2 * sizeof(double));
	for (n = 0; n < nd_in; n++) {
		pt[n] = orb->points[n].pt;
		px[n] = orb->points[n].px;
		py[n] = orb->points[n].py;
		pz[n] = orb->points[n].pz;
	}

	/* extrapolate the orbit into the longer structure */
	Nodr = 4;
	extrapolate_poly(pt, px, ppt, ppx, ppvx, nd_in, nd_add, Nodr);
	extrapolate_poly(pt, py, ppt, ppy, ppvy, nd_in, nd_add, Nodr);
	extrapolate_poly(pt, pz, ppt, ppz, ppvz, nd_in, nd_add, Nodr);
	for (n = 0; n < nd_add; n++) {
		orb_out->points[n].pt = ppt[n];
		orb_out->points[nd_out - n - 1].pt = ppt[2 * nd_add - n - 1];
		orb_out->points[n].px = ppx[n];
		orb_out->points[nd_out - n - 1].px = ppx[2 * nd_add - n - 1];
		orb_out->points[n].py = ppy[n];
		orb_out->points[nd_out - n - 1].py = ppy[2 * nd_add - n - 1];
		orb_out->points[n].pz = ppz[n];
		orb_out->points[nd_out - n - 1].pz = ppz[2 * nd_add - n - 1];
		orb_out->points[n].vx = ppvx[n];
		orb_out->points[nd_out - n - 1].vx = ppvx[2 * nd_add - n - 1];
		orb_out->points[n].vy = ppvy[n];
		orb_out->points[nd_out - n - 1].vy = ppvy[2 * nd_add - n - 1];
		orb_out->points[n].vz = ppvz[n];
		orb_out->points[nd_out - n - 1].vz = ppvz[2 * nd_add - n - 1];
	}

	/* replace the original state vector in the interior of the output vector */
	k = 0;
	for (n = nd_add; n < nd_in + nd_add; n++) {
		orb_out->points[n].pt = orb->points[k].pt;
		orb_out->points[n].px = orb->points[k].px;
		orb_out->points[n].py = orb->points[k].py;
		orb_out->points[n].pz = orb->points[k].pz;
		orb_out->points[n].vx = orb->points[k].vx;
		orb_out->points[n].vy = orb->points[k].vy;
		orb_out->points[n].vz = orb->points[k].vz;
		k++;
	}

	/* write the extrapolated LED file */
	if ((outfile = fopen(argv[2], "w")) == NULL)
		die("Can't open ", argv[2]);
	write_orb(outfile, orb_out);

	fclose(infile);
	fclose(outfile);
	exit(0);
}

void extrapolate_poly(double *T, double *Y, double *TT, double *YY, double *DYY, int nd_in, int nd_add, int Nodr) {

	int i;
	double dt, mt, my, coef[6], t_tmp, dt_tmp, y_tmp, y_tmp1, y_tmp2;
	double polyval(double, double *, int);
	double sum_array(double *, int);

	dt_tmp = 0.2;

	if (nd_in < 5) {
		fprintf(stderr, "Not enough points to estimate polynomial coefficitnets\n");
		exit(0);
	}

	mt = sum_array(T, nd_in) / nd_in;
	my = sum_array(Y, nd_in) / nd_in;
	for (i = 0; i < nd_in; i++) {
		T[i] = T[i] - mt;
		Y[i] = Y[i] - my;
	}

	dt = T[nd_in / 2 + 1] - T[nd_in / 2];
	polyfit(T, Y, coef, &nd_in, &Nodr);

	/* interpolate nd_add points in the front and at the back */
	for (i = 0; i < nd_add; i++) {
		t_tmp = T[0] - dt * (nd_add - i);
		y_tmp = polyval(t_tmp, coef, Nodr);
		y_tmp1 = polyval(t_tmp - dt_tmp, coef, Nodr);
		y_tmp2 = polyval(t_tmp + dt_tmp, coef, Nodr);

		TT[i] = t_tmp + mt;
		YY[i] = y_tmp + my;
		DYY[i] = (y_tmp2 - y_tmp1) / 2 / dt_tmp;

		t_tmp = T[nd_in - 1] + dt * (i + 1);
		y_tmp = polyval(t_tmp, coef, Nodr);
		y_tmp1 = polyval(t_tmp - dt_tmp, coef, Nodr);
		y_tmp2 = polyval(t_tmp + dt_tmp, coef, Nodr);

		TT[nd_add + i] = t_tmp + mt;
		YY[nd_add + i] = y_tmp + my;
		DYY[nd_add + i] = (y_tmp2 - y_tmp1) / 2 / dt_tmp;
	}

	for (i = 0; i < nd_in; i++) {
		T[i] = T[i] + mt;
		Y[i] = Y[i] + my;
	}
}

double polyval(double T, double *C, int Nodr) {

	int j;
	double y_tmp;

	y_tmp = 0;
	for (j = 0; j < Nodr; j++) {
		y_tmp = y_tmp + C[j] * pow(T, j);
	}

	return (y_tmp);
}

double sum_array(double *a, int num_elements) {

	int i;
	double sum = 0;
	for (i = 0; i < num_elements; i++) {
		sum = sum + a[i];
	}
	return (sum);
}
