/****************************************************************************
 *  Subroutine to project a longitude, latitude, and topography
 *  into a file of range, azimuth, and topography.
 *  The basic approach is to :
 *   1 Read the header of the master radar image.  This supplies
 *     the radar co-ordinate system as well as the start and stop
 *     times for the orbit calculation.
 *   2 Read the topography data and convert to xyz positions.
 *   3 Fly the satellite along its orbit and determine the time
 *     of closest approach to each of the xyz points.   The program
 *     must have local access to the LED-file for the master.
 *   4 For each of the xyz points, calculate range, azimuth and topography
 ****************************************************************************/

/***************************************************************************
 * Creator:  Xiaopeng Tong and David T. Sandwell                           *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  08/10/2006                                                    *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * DATE                                                                    *
 * 12/15/07 - modified to work with complete grids                         *
 * 01/29/08 - modified to increase speed using the Golden Section Search   *
 * 04/13/08 - modified to give muiltiple choice of output                  *
 *            (single,double of binary or acsii)                           *
 * The algorithm is called Golden Section Search in One                    *
 * dimensional. Refer to Numerical Recipes for further details.            *
 * There is a deadly error in the NR edit 2nd:                             *
 * SHFT3(x0,x1,x2,R*x1+C*x3)                                               *
 * should be: SHFT3(x0,x1,x2,R*x3+C*x1)                                    *
 * Same error occurs at the nearby line.(2 lines below)                    *
 * The inputs for the program are 3 coordinates of a point in space        *
 * The function of the distance is in the same file.                       *
 * The function of the orbit is achieved by getorb_alos_                   *
 * the outputs are minimum range from the orbit to the point and the time. *
 * 09/17/08 - modified to read the orbit position all in once into an      *
 * array to speed up. That is to say, modify both the goldop subrountine   *
 * to use the orb_position array and the getorb_alos in the  main program  *
 * to read the array.                                                      *
 * 06/15/09 - modifed to correct the contradiction of range sampling rate  *
 * between FBS and FBD mode.
 ****************************************************************************/

#include "gmtsar.h"
#include "llt2xyz.h"
#include "orbit.h"

#define R 0.61803399
#define C 0.382
#define SHFT2(a, b, c)                                                                                                           \
	(a) = (b);                                                                                                                   \
	(b) = (c);
#define SHFT3(a, b, c, d)                                                                                                        \
	(a) = (b);                                                                                                                   \
	(b) = (c);                                                                                                                   \
	(c) = (d);
#define TOL 3
void set_prm_defaults(struct PRM *);
void read_orb(FILE *, struct SAT_ORB *);
void hermite_c(double *, double *, double *, int, int, double, double *, int *);

void llt2rat_sub(struct PRM *prm, double *target_llt, double *target_rat) {

	double rln, rlt, rht, dr, t1, t2, tm;
	// double ts,rng, thet, relp, telp;
	double ts, rng;
	double xp[3];
	double xt[3];
	double rp[3];
	/* double r0,rf,a0,af; */
	// double rad=PI/180.;
	double fll, rdd, daa, drr;
	int j, nrec, npad = 8000;
	int goldop();
	int stai, endi, midi;
	double **orb_pos;
	//struct PRM prm;
	struct SAT_ORB *orb;
	FILE *ldrfile;
	//FILE *fprm1;
	int calorb_alos(struct SAT_ORB *, double **orb_pos, double ts, double t1, int nrec);
	double rsr;
	char value[128], name[128];

	/*  get the orbit data */
	ldrfile = fopen(prm->led_file, "r");
	if (ldrfile == NULL)
		die("can't open LED file", prm->led_file);
	orb = (struct SAT_ORB *)malloc(sizeof(struct SAT_ORB));
	read_orb(ldrfile, orb);

	dr = 0.5 * SOL / prm->fs;
#if 0
        r0 = -10.;
        rf = prm->num_rng_bins + 10.;
        a0 = -20.;
        af = prm->num_patches * prm->num_valid_az + 20.;
#endif
	/* compute the flattening */
	fll = (prm->ra - prm->rc) / prm->ra;

	/* compute the start time, stop time and increment */
	t1 = 86400. * prm->clock_start + (prm->nrows - prm->num_valid_az) / (2. * prm->prf);
	t2 = t1 + prm->num_patches * prm->num_valid_az / prm->prf;

	/* sample the orbit only every 2th point or about 8 m along track */
	ts = 2. / prm->prf;
	nrec = (int)((t2 - t1) / ts);

	/* allocate memory for the orbit postion into a 2-dimensional array. It's
	   about nrec+-npad * 4 * 4bytes = 320 KB the first "4" means space and time
	   dimension : t, x, y, z         */

	/* allocate storage for an array of pointers  */
	orb_pos = malloc(4 * sizeof(double *));

	/* for each pointer, allocate storage for an array of floats  */
	for (j = 0; j < 4; j++) {
		orb_pos[j] = malloc((nrec + 2 * npad) * sizeof(double));
	}

	/* read in the postion of the orbit */
	(void)calorb_alos(orb, orb_pos, ts, t1, nrec);

	/* read the llt points and convert to xyz.  */
	rlt = target_llt[0];
	rln = target_llt[1];
	rht = target_llt[2];
	rp[0] = rlt;
	rp[1] = rln;
	rp[2] = rht;
	plh2xyz(rp, xp, prm->ra, fll);
	xt[0] = -1.0;

	/* compute the topography due to the difference between the local radius and
	 * center radius */
	// thet = rlt * rad;
	// if(prm.rc > 6350000. && prm.ra > 6350000. && prm.RE > 6350000.) {
	// relp=1./sqrt((cos(thet)/prm.ra)*(cos(thet)/prm.ra)+(sin(thet)/prm.rc)*(sin(thet)/prm.rc));
	// telp=relp-prm.RE;
	//}
	// else {
	// telp=0.;
	//}
	// rp[2]=rht + telp;
	rp[2] = sqrt(xp[0] * xp[0] + xp[1] * xp[1] + xp[2] * xp[2]) - prm->RE;

	/* minimum for each point */
	stai = 0;
	endi = nrec + npad * 2 - 1;
	midi = (stai + (endi - stai) * C);
	(void)goldop(ts, t1, orb_pos, stai, endi, midi, xp[0], xp[1], xp[2], &rng, &tm);
	xt[0] = rng;
	xt[1] = tm;

	/* compute the range and azimuth in pixel space and correct for an azimuth
	 * bias*/
	xt[0] = (xt[0] - prm->near_range) / dr - (prm->rshift + prm->sub_int_r) + prm->chirp_ext;
	xt[1] = prm->prf * (xt[1] - t1) - (prm->ashift + prm->sub_int_a);

	/* compute the azimuth and range correction if the Doppler is not zero */

	if (prm->fd1 != 0.) {
		rdd = (prm->vel * prm->vel) / rng;
		daa = -0.5 * (prm->lambda * prm->fd1) / rdd;
		drr = 0.5 * rdd * daa * daa / dr;
		daa = prm->prf * daa;
		xt[0] = xt[0] + drr;
		xt[1] = xt[1] + daa;
	}

	target_rat[0] = xt[0];
	target_rat[1] = xt[1];
	target_rat[2] = rp[2];

	free(orb_pos);
	free(orb);
}

/*    subfunctions    */

int goldop(double ts, double t1, double **orb_pos, int ax, int bx, int cx, double xpx, double xpy, double xpz, double *rng,
           double *tm) {

	/* use golden section search to find the minimum range between the target and
	 * the orbit */
	/* xpx, xpy, xpz is the position of the target in cartesian coordinate */
	/* ax is stai; bx is endi; cx is midi it's easy to tangle */

	double f1, f2;
	int x0, x1, x2, x3;
	int xmin;
	double dist();

	x0 = ax;
	x3 = bx;
	// if (fabs(bx-cx) > fabs(cx-ax)) {
	if (abs(bx - cx) > abs(cx - ax)) {
		x1 = cx;
		x2 = cx + (int)fabs((C * (bx - cx)));
	}
	else {
		x2 = cx;
		x1 = cx - (int)fabs((C * (cx - ax))); /* make x0 to x1 the smaller segment */
	}

	f1 = dist(xpx, xpy, xpz, x1, orb_pos);
	f2 = dist(xpx, xpy, xpz, x2, orb_pos);

	while ((x3 - x0) > TOL) {
		if (f2 < f1) {
			SHFT3(x0, x1, x2, (int)(R * x3 + C * x1));
			SHFT2(f1, f2, dist(xpx, xpy, xpz, x2, orb_pos));
		}
		else {
			SHFT3(x3, x2, x1, (int)(R * x0 + C * x2));
			SHFT2(f2, f1, dist(xpx, xpy, xpz, x1, orb_pos));
		}
	}
	if (f1 < f2) {
		xmin = x1;
		*tm = orb_pos[0][x1];
		*rng = f1;
	}
	else {
		xmin = x2;
		*tm = orb_pos[0][x2];
		*rng = f2;
	}

	return (xmin);
}

double dist(double x, double y, double z, int n, double **orb_pos) {

	double d, dx, dy, dz;

	dx = x - orb_pos[1][n];
	dy = y - orb_pos[2][n];
	dz = z - orb_pos[3][n];
	d = sqrt(dx * dx + dy * dy + dz * dz);

	return (d);
}

int calorb_alos(struct SAT_ORB *orb, double **orb_pos, double ts, double t1, int nrec)
/* function to calculate every position in the orbit   */

{
	int i, k, nval;
	int npad = 8000;   /* number of buffer points to add before and after the
	                      acquisition */
	int ir;            /* return code: 0 = ok; 1 = interp not in center; 2 = time out of
	                      range */
	double xs, ys, zs; /* position at time */
	double *pt, *px, *py, *pz, *pvx, *pvy, *pvz;
	double pt0;
	double time;

	px = (double *)malloc(orb->nd * sizeof(double));
	py = (double *)malloc(orb->nd * sizeof(double));
	pz = (double *)malloc(orb->nd * sizeof(double));
	pvx = (double *)malloc(orb->nd * sizeof(double));
	pvy = (double *)malloc(orb->nd * sizeof(double));
	pvz = (double *)malloc(orb->nd * sizeof(double));
	pt = (double *)malloc(orb->nd * sizeof(double));

	pt0 = 86400. * orb->id + orb->sec;
	for (k = 0; k < orb->nd; k++) {
		pt[k] = pt0 + k * orb->dsec;
		px[k] = orb->points[k].px;
		py[k] = orb->points[k].py;
		pz[k] = orb->points[k].pz;
		pvx[k] = orb->points[k].vx;
		pvy[k] = orb->points[k].vy;
		pvz[k] = orb->points[k].vz;
	}

	nval = 6;

	/* loop to get orbit position of every point and store them into orb_pos */
	for (i = 0; i < nrec + npad * 2; i++) {
		time = t1 - npad * ts + i * ts;
		orb_pos[0][i] = time;

		hermite_c(pt, px, pvx, orb->nd, nval, time, &xs, &ir);
		hermite_c(pt, py, pvy, orb->nd, nval, time, &ys, &ir);
		hermite_c(pt, pz, pvz, orb->nd, nval, time, &zs, &ir);

		orb_pos[1][i] = xs;
		orb_pos[2][i] = ys;
		orb_pos[3][i] = zs;
	}

	free((double *)px);
	free((double *)py);
	free((double *)pz);
	free((double *)pt);
	free((double *)pvx);
	free((double *)pvy);
	free((double *)pvz);

	return orb->nd;
}
