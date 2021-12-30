/*	$Id$	*/
/***************************************************************************
 *  Program to project a longitude, latitude, and topography
 *  into a file of longitude latitude, and topography and unit look vector.
 *  The basic approach is to :
 *   1 Read the header of the master radar image.  This supplies
 *     the radar co-ordinate system as well as the start and stop
 *     times for the orbit calculation.
 *   2 Read the topography data and convert to xyz positions.
 *   3 Fly the satellite along its orbit and determine the time
 *     of closest approach to each of the xyz points.   The program
 *     must have local access to the LED-file for the master.
 *   4 For each of the xyz points, calculate range, azimuth and topography and
 *     unit look vector then transform the unit look vector from global to
 *     local coordinate
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
 * 11/28/08 - modified from ALOS_llt2rat.c to compute accurate look vector.*
 * need to specify the latitude and longitude of the pole to               *
 * transform from global elipsoid cartesian to local East-North-           *
 * up coordinate                                                           *
 *            The start and end clock of the orbit is changed to allow for *
 *            long swath, the start time should be in  the first frame so  *
 *            PRM file.                                                    *
 * 12/03/10 - modified to work with ENVISAT, Matt Wei                      *
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

char *USAGE = " \n Usage: "
              "SAT_look master.PRM [-bo[s|d]] < inputfile > outputfile  \n\n"
              "             master.PRM   -  parameter file for master image and points "
              "to LED orbit file \n"
              "             inputfile    -  lon, lat, elevation [ASCII] \n"
              "             outputfile   -  lon, lat, elevation look_E look_N look_U "
              "[ASCII default] \n"
              "             -bos or -bod -  binary single or double precision output \n\n"
              "  example: SAT_look master.PRM < topo.llt > topo.lltn    \n\n"
              "  Note that the output elevation is the one above reference radius "
              "specified in the PRM file\n";

void read_orb(FILE *, struct SAT_ORB *);
void hermite_c(double *, double *, double *, int, int, double, double *, int *);
void set_prm_defaults(struct PRM *);

int main(int argc, char **argv) {

	FILE *fprm1;
	int otype;
	double rln, rlt, rht, t1, t2, tm;
	double ts, rng, thet, relp, telp;
	double xp[3];
	/* double xt[3],dr; */
	double rp[3];
	double dd[6]; /* dummy for output  double precision */
	float ds[6];  /* dummy for output  single precision */
	/* double r0,rf,a0,af; */
	double rad = PI / 180.;
	double fll;
	int i, j, k, nrec, npad = 8000;
	int goldop();
	int stai, endi, midi, xmin;
	double **orb_pos;
	struct PRM prm;
	struct SAT_ORB *orb;
	FILE *ldrfile;
	int calorb_alos(struct SAT_ORB *, double **orb_pos, double ts, double t1, int nrec);
	double len, unit_x, unit_y, unit_z;
	double dist();

	double Rx[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
	double Rz[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
	double A[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
	double local_look[3] = {0, 0, 0};
	double pole_lon, pole_lat, b, g;

	/* Make sure usage is correct and files can be opened  */
	if (argc < 2 || argc > 3) {
		fprintf(stderr, "%s\n", USAGE);
		exit(-1);
	}

	/* otype:    1 -- ascii; 2 -- single precision binary; 3 -- double precision
	 * binary    */
	otype = 1;
	if (argc == 3) {
		if (!strcmp(argv[2], "-bos"))
			otype = 2;
		else if (!strcmp(argv[2], "-bod"))
			otype = 3;
		else {
			fprintf(stderr, " %s *** option not recognized ***\n\n", argv[2]);
			fprintf(stderr, "%s", USAGE);
			exit(1);
		}
	}

	/*  open and read the parameter file */
	if ((fprm1 = fopen(argv[1], "r")) == NULL) {
		fprintf(stderr, "couldn't open master.PRM \n");
		fprintf(stderr, "%s\n", USAGE);
		exit(-1);
	}
	/* initialize the prm file   */
    null_sio_struct(&prm);
	set_prm_defaults(&prm);

	get_sio_struct(fprm1, &prm);

	/* dubugging
	        fprintf(stderr,"debug\n");
	        fprintf(stderr,"%s\n",prm.led_file);
	        fprintf(stderr,"%lf\n",prm.fs);
	        fprintf(stderr,"%d\n",prm.num_rng_bins);
	        fprintf(stderr,"%d\n",prm.num_valid_az);
	        fprintf(stderr,"%lf\n",prm.rc);
	        fprintf(stderr,"%lf\n",prm.ra);
	        fprintf(stderr,"%lf\n",prm.prf);
	        fprintf(stderr,"%lf\n",prm.RE);
	        fprintf(stderr,"%lf\n",prm.sub_int_a);
	        fprintf(stderr,"%lf\n",prm.SC_clock_start);
	        fprintf(stderr,"%lf\n",prm.clock_start);
	        fprintf(stderr,"%d\n",prm.nrows);
	        fprintf(stderr,"%d\n",prm.num_patches);
	        fprintf(stderr,"%lf\n",prm.near_range);
	        fprintf(stderr,"rshift = %d\n",prm.rshift);
	        fprintf(stderr,"sub_int_r = %lf\n",prm.sub_int_r);
	        fprintf(stderr,"chirp_ext = %d\n",prm.chirp_ext);
	*/

	fclose(fprm1);

	/*  get the orbit data */
	ldrfile = fopen(prm.led_file, "r");
	if (ldrfile == NULL)
		die("can't open ", prm.led_file);
	orb = (struct SAT_ORB *)malloc(sizeof(struct SAT_ORB));
	read_orb(ldrfile, orb);

#if 0
	dr = 0.5*SOL/prm.fs;
        r0 = -10.;
        rf = prm.num_rng_bins + 10.;
        a0 = -20.;
        af = prm.num_patches * prm.num_valid_az + 20.;
#endif
	/* compute the flattening */

	fll = (prm.ra - prm.rc) / prm.ra;

	/* compute the start time, stop time and increment */

	t1 = 86400. * prm.clock_start + (prm.nrows - prm.num_valid_az) / (2. * prm.prf);
	t2 = t1 + prm.num_patches * prm.num_valid_az / prm.prf;

	/* sample the orbit only every 2th point or about 8 m along track */

	ts = 2. / prm.prf;
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
	/*        fprintf(stderr,"ts=%lf,t1=%lf,nrec=%d\n",ts,t1,nrec);   */

	(void)calorb_alos(orb, orb_pos, ts, t1, nrec);

	/* look at the orb_pos
	        for(i=0;i<4;i++) {
	          for (j=0;j<nrec+2*npad;j=j+100) {
	            fprintf(stderr,"%d,  %lf\n",j,orb_pos[i][j]); }
	        }
	*/

	/* read the llt points and convert to xyz.  */

	k = 0;
	while (scanf(" %lf %lf %lf ", &rln, &rlt, &rht) == 3) {
		rp[0] = rlt;
		rp[1] = rln;
		rp[2] = rht;
		plh2xyz(rp, xp, prm.ra, fll);
		if (rp[1] > 180.)
			rp[1] = rp[1] - 360.;
		/* xt[0]=-1.0;  */

		/* compute the topography due to the difference between the local radius and
		 * center radius */

		thet = rlt * rad;
		if (prm.rc > 6350000. && prm.ra > 6350000. && prm.RE > 6350000.) {
			relp = 1. / sqrt((cos(thet) / prm.ra) * (cos(thet) / prm.ra) + (sin(thet) / prm.rc) * (sin(thet) / prm.rc));
			telp = relp - prm.RE;
		}
		else {
			telp = 0.;
		}
		rp[2] = rht + telp;

		/* minimum for each point */

		stai = 0;
		endi = nrec + npad * 2 - 1;
		midi = (stai + (endi - stai) * C);

		xmin = goldop(ts, t1, orb_pos, stai, endi, midi, xp[0], xp[1], xp[2], &rng, &tm);

		/* xt[0]=rng;
		xt[1]=tm; */

		/*    fprintf(stderr,"rng = %f, tm = %f",rng,tm);   */

		/*    calculate the unit look vector originated from ground point to
		   satellite in the global cartesian coordinate              */
		/*
		  fprintf(stderr,"x y z of the satellite minimum: %f, %f,
		  %f\n",orb_pos[1][xmin],orb_pos[2][xmin],orb_pos[3][xmin]);   */

		len = dist(xp[0], xp[1], xp[2], xmin, orb_pos);
		unit_x = (orb_pos[1][xmin] - xp[0]) / len;
		unit_y = (orb_pos[2][xmin] - xp[1]) / len;
		unit_z = (orb_pos[3][xmin] - xp[2]) / len;

		/*      fprintf(stderr,"unit look vector: %f, %f,
		 * %f\n",unit_x,unit_y,unit_z);   */

		/*  there is a way to test the sign is right  */
		/* compute the look vector in local coordinate  */

		/* rotate the global coordinate to a local coordinate system  */

		/*   location of the rotation pole   */
		/*   refer to a matlab version rot3d.m to check the formula  */

		pole_lon = fmod((rln + 360.0), 360.0);

		/*     fprintf(stderr,"pole_lon %f\n",pole_lon);

		*/
		pole_lat = rlt;

		/*   construct a rotation matrix   */
		b = (pole_lat - 90) / 180. * 3.14159;
		Rx[0][0] = 1;
		Rx[0][1] = 0;
		Rx[0][2] = 0;
		Rx[1][0] = 0;
		Rx[1][1] = cos(b);
		Rx[1][2] = -sin(b);
		Rx[2][0] = 0;
		Rx[2][1] = sin(b);
		Rx[2][2] = cos(b);

		g = (pole_lon + 90) / 180. * 3.14159;
		Rz[0][0] = cos(g);
		Rz[0][1] = sin(g);
		Rz[0][2] = 0;
		Rz[1][0] = -sin(g);
		Rz[1][1] = cos(g);
		Rz[1][2] = 0;
		Rz[2][0] = 0;
		Rz[2][1] = 0;
		Rz[2][2] = 1;

		/* zero A */
		for (i = 0; i < 3; i++)
			for (j = 0; j < 3; j++)
				A[i][j] = 0;

		for (i = 0; i < 3; i++)
			for (j = 0; j < 3; j++)
				for (k = 0; k < 3; k++)
					A[i][j] += Rx[i][k] * Rz[k][j];

		/* print out the matrix                         */

		/*
		         for (i=0;i<3;i++)
		         {
		           for (j=0;j<3;j++)
		             fprintf(stderr,"%f  ",Rx[i][j]);
		           putchar('\n');
		         }

		        for (i=0;i<3;i++)
		         {
		           for (j=0;j<3;j++)
		             fprintf(stderr,"%f  ",Rz[i][j]);
		           putchar('\n');
		         }

		        for (i=0;i<3;i++)
		         {
		           for (j=0;j<3;j++)
		             fprintf(stderr,"%f  ",A[i][j]);
		           putchar('\n');
		         }
		*/

		for (i = 0; i < 3; i++)
			local_look[i] = A[i][0] * unit_x + A[i][1] * unit_y + A[i][2] * unit_z;

		/* compute the range and azimuth in pixel space
		           xt[0] = (xt[0] -
		   prm.near_range)/dr-(prm.rshift+prm.sub_int_r)+prm.chirp_ext; xt[1] =
		   prm.prf*(xt[1]-t1)-(prm.ashift+prm.sub_int_a); fprintf(stderr,"xt[0] =
		   %f\n",xt[0]); if(xt[0] < r0 || xt[0] > rf || xt[1] < a0 || xt[1] > af)
		   continue; fprintf(stderr,"%f %f %f %f %f
		   \n",xt[0],xt[1],rp[2],rp[1],rp[0]);         */

		if (otype == 1) {
			/*	     fprintf(stdout,"%f %f %f %f %f
			 * %f\n",rp[1],rp[0],rp[2],unit_x,unit_y,unit_z);  */
			fprintf(stdout, "%f %f %f %f %f %f\n", rp[1], rp[0], rp[2], local_look[0], local_look[1], local_look[2]);
		}
		else if (otype == 2) {
			ds[0] = (float)rp[1];
			ds[1] = (float)rp[0];
			ds[2] = (float)rp[2];
			ds[3] = (float)local_look[0];
			ds[4] = (float)local_look[1];
			ds[5] = (float)local_look[2];
			fwrite(ds, sizeof(float), 6, stdout);
		}
		else if (otype == 3) {
			dd[0] = rp[1];
			dd[1] = rp[0];
			dd[2] = rp[2];
			dd[3] = local_look[0];
			dd[4] = local_look[1];
			dd[5] = local_look[2];
			fwrite(dd, sizeof(double), 6, stdout);
		}
	}

	/* free the orb_pos array  */
	for (j = 0; j < 4; j++) {
		free(orb_pos[j]);
	}
	free(orb_pos);
	free(orb);
	return (0);
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

	while ((x3 - x0) > TOL && fabs(f1-f2) > 1e-6) {
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
