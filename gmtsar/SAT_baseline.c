/*	$Id$	*/
/*******************************************************************************
 * Compute the interferometric baseline between reference and repeat           *
 * images.  Also estimates the crude offset for image matching.                *
 *******************************************************************************/
/********************************************************************************
 * Creator:  Matt Wei, 04/26/10 *
 *										*
 * Based on the program ALOS_baseline by * Creator:  Sandwell and Rob Mellors
 ** (San Diego State University, Scripps Institution of Oceanography)  * Date :
 *06/07/2007                                                         *
 * Modifications: * 29/26/2018 by Xiaohua Xu * Adding the along track component
 *for Baseline called B_offset      * This further ensures phase closure *
 ********************************************************************************/

#include "gmtsar.h"
#include "orbit.h"

char *USAGE = "Usage: (two modes)\n"
              "mode 1:\n\n"
              "SAT_baseline PRM_master PRM_master\n\n"
              "This is used to compute height information\n"
              "(writes out height information for appending to PRM file)\n"
              "\nmode 2:\n\n"
              "SAT_baseline PRM_master PRM_aligned \n\n"
              "PRM_master 	   PRM file for reference image\n"
              "PRM_aligned 	   PRM file of secondary image\n"
              "Please make sure the orbit file data is in PRM \n"
              "Program runs through repeat orbit to find nearest point \n"
              "to the start, center and end on the reference orbit\n"
              "(writes out parameters for appending to PRM file)\n";
/*
"\nor mode 2:\n\n"
"SAT_baseline -input file\n\n"
"file: list of PRM files\n"
"first file is assumed to be master\n"
"following are aligneds\n"
"(writes out decimal year, Bperp, and PRM name)\n";
 */

/* function prototypes */
double find_dist(double, double, double, double, double, double);
double find_alpha_degrees(double, double);
void find_unit_vectors(double, double, double, double *, double *, double *, double *);
void endpoint_distance(int, double, double, double, double, double *, double *, double *, double *, int *);
void find_parallel_perp_baseline(struct PRM, struct PRM, double, double *, double *);
void write_prm_baseline(struct PRM);
void get_sign(struct PRM, double, double, double, double, int *);
void write_bperp(struct PRM, char *);
void baseline_parse_command_line(char **, int *, int *);
void read_input_file(char *, int, char **);
void baseline(struct PRM *, struct SAT_ORB *, int, int, char **, double);
void read_all_ldr(struct PRM *, struct SAT_ORB *, int);
void read_orb(FILE *, struct PRM *, struct SAT_ORB *);
void llt2rat_sub(char *, double *, double *);
void interpolate_SAT_orbit_slow(struct SAT_ORB *orb, double time, double *, double *, double *, int *);
void interpolate_SAT_orbit(struct SAT_ORB *, double *, double *, double *, double, double *, double *, double *, int *);
void calc_height_velocity(struct SAT_ORB *, struct PRM *, double, double, double *, double *, double *, double *, double *);
void polyfit(double *, double *, double *, int *, int *);

int main(int argc, char **argv) {
	int i;
	int input_flag, nfiles;
	double fs0 = 0.;

	struct PRM *r;       /* reference orbit is 0; repeats are > 0*/
	struct SAT_ORB *orb; /* reference orbit is 0; repeats are > 0*/
	char **filename;
	FILE *prmfile;

	if (argc < 3)
		die(USAGE, "");

	verbose = 0;
	input_flag = 0;

	/* read command line */
	baseline_parse_command_line(argv, &nfiles, &input_flag);

	filename = malloc(nfiles * sizeof(char *));
	for (i = 0; i < nfiles; i++)
		filename[i] = malloc(128 * sizeof(char));

	r = malloc(nfiles * sizeof(struct PRM));

	if (input_flag == 1)
		read_input_file(argv[2], nfiles, filename);

	for (i = 0; i < nfiles; i++) {
		if (input_flag == 0)
			strcpy(filename[i], argv[i + 1]);
		if ((prmfile = fopen(filename[i], "r")) == NULL)
			die("Can't open prmfile ", filename[i]);

        null_sio_struct(&r[i]);
		get_sio_struct(prmfile, &r[i]);

		/* check whether range sampling rate is same. If not, warning */
		if (i == 0)
			fs0 = r[0].fs;
		if (r[i].fs != fs0) {
			fprintf(stderr, "\nWARNING:\nRange_sampling_rate is not consistant.\nYou "
			                "need to do FBD/FBS conversion.\n\n");
		}

		fclose(prmfile);
	}

	printf("SC_identity = %d \n", r[0].SC_identity);
	orb = malloc(nfiles * sizeof(struct SAT_ORB));
	read_all_ldr(r, orb, nfiles);
	baseline(r, orb, nfiles, input_flag, filename, fs0);

	return (EXIT_SUCCESS);
}
/*---------------------------------------------------------------------------*/
/* Get time info and find the orbit */
void read_all_ldr(struct PRM *r, struct SAT_ORB *orb, int nfiles) {
	int i;
	FILE *ldrfile;

	for (i = 0; i < nfiles; i++) {
		if (i == 0)
			printf("......master LED file %s \n", r[0].led_file);
		if (i != 0)
			printf(".........aligned LED file %s \n", r[i].led_file);

		/* open each ldrfile and read into structure r */
		if ((ldrfile = fopen(r[i].led_file, "r")) == NULL)
			die("Can't open ldrfile %s", r[i].led_file);

		read_orb(ldrfile, &r[i], &orb[i]);

		fclose(ldrfile);
	}
}

/* interpolation with polynomial refinement */
void poly_interp(struct SAT_ORB *orb, double *t, double *b, double t_ref, double x0, double y0, double z0, double shift) {

	int ntt = 100, k, ir;
	double *time, *bs;
	double ddt, xs, ys, zs;
	double d3[3];
	int nc = 3;

	ddt = 0.01 / ntt;
	time = malloc(ntt * sizeof(double));
	bs = malloc(ntt * sizeof(double));
	for (k = 0; k < ntt; k++) {
		time[k] = (k - ntt / 2 + 0.5) * ddt;
		interpolate_SAT_orbit_slow(orb, t_ref + time[k] + shift, &xs, &ys, &zs, &ir);
		bs[k] = find_dist(xs, ys, zs, x0, y0, z0);
		bs[k] = bs[k] * bs[k];
	}
	polyfit(time, bs, d3, &ntt, &nc);
	*t = t_ref + shift - d3[1] / (2.0 * d3[2]);
	*b = sqrt(d3[0] - d3[1] * d3[1] / 4.0 / d3[2]);

	free(time);
	free(bs);
}

/*---------------------------------------------------------------------------*/
void baseline(struct PRM *r, struct SAT_ORB *orb, int nfiles, int input_flag, char **filename, double fs0) {

	int ii, nd, ir;
	int k, ns, ns2, m1, m2, m3, sign1, sign2, sign3;
	double dr, dt, ds;
	double bv1, bv2, bv3, bh1, bh2, bh3;
	double t11, t12, t13, t21, t22, t23;
	double x11, y11, z11, x12, y12, z12, x13, y13, z13;
	double x21, y21, z21, x22, y22, z22, x23, y23, z23;
	double ru1, ru2, ru3, xu1, yu1, zu1, xu2, yu2, zu2, xu3, yu3, zu3;
	double ts, xs, ys, zs;
	double b1, b2, b3, bpara, bperp; //, b_tmp;
	double *pt, *p, *pv, pt0;
	double height, re_c, vg, vtot, rdot;
	double radar_look[3] = {0, 0, 0};
	double re, rho, b, theta;
	double t_vel, t_sta, x_vel, x_sta, y_vel, y_sta, z_vel, z_sta, ru33, xu33, yu33, zu33, o1x[3], o1y[3], o1z[3], ru4, xu4, yu4,
	    zu4;
	double glob_look[3], target[3], target_llt[3], fll;
	double target_rat_ref[3] = {0, 0, 0};
	double target_rat_rep[3] = {0, 0, 0};
	double far_range;
	double dstart;
	// double d1=0.0,d2=0.0,d3=0.0;

	/* work on the reference orbit */
	get_seconds(r[0], &t11, &t12);
	t13 = (t11 + t12) / 2.;
	dr = 0.5 * SOL / fs0;
	dt = 0.5 / r[0].prf;
	ns = (int)((t12 - t11) / dt); /* seconds of frame */
	dt = (t12 - t11) / (ns - 1);

	/* set the extension to 50% unless ERS or Envisat and then set to 200% */
	ns2 = ns * 0.5;
	if (r[0].SC_identity < 3 || r[0].SC_identity == 4 || r[0].SC_identity == 10)
		ns2 = ns * 2;

	nd = orb[0].nd;
	for (ii = 1; ii < nfiles; ii++) {
		if (nd < orb[ii].nd)
			nd = orb[ii].nd;
	}
	/* compute reference start, center and end point */
	interpolate_SAT_orbit_slow(&orb[0], t11, &x11, &y11, &z11, &ir);
	interpolate_SAT_orbit_slow(&orb[0], t12, &x12, &y12, &z12, &ir);
	interpolate_SAT_orbit_slow(&orb[0], t13, &x13, &y13, &z13, &ir);

	/* allocate memory */
	pt = malloc(nd * sizeof(double));
	p = malloc(nd * sizeof(double));
	pv = malloc(nd * sizeof(double));

	/* work on the repeat orbit */
	ii = 1;

	get_seconds(r[ii], &t21, &t22);
	t23 = (t21 + t22) / 2.;
	pt0 = (24.0 * 60.0 * 60.0) * orb[ii].id + orb[ii].sec;
	for (k = 0; k < nd; k++)
		pt[k] = pt0 + k * orb[ii].dsec;
	verbose = 0;

	/* look at other orbit information and recalculate the height using the ashift
	 */
	calc_height_velocity(&orb[ii], &r[ii], t22, t22, &height, &re_c, &vg, &vtot, &rdot);
	r[ii].ht_end = height + re_c - r[ii].RE;
	calc_height_velocity(&orb[ii], &r[ii], t23, t23, &height, &re_c, &vg, &vtot, &rdot);
	r[ii].ht = height + re_c - r[ii].RE;
	calc_height_velocity(&orb[ii], &r[ii], t21, t21, &height, &re_c, &vg, &vtot, &rdot);
	r[ii].ht_start = height + re_c - r[ii].RE;

	/* loop over repeat orbit 				*/
	/* add 2 scene lengths (ns2) for buffer at each end 	*/
	b1 = b2 = b3 = -1.0;

	pt0 = (24.0 * 60.0 * 60.0) * orb[ii].id + orb[ii].sec;
	for (k = 0; k < nd; k++)
		pt[k] = pt0 + k * orb[ii].dsec;

	/* set some defeault values 				*/

	m1 = -99999;
	x21 = y21 = z21 = -99999.0;
	x22 = y22 = z22 = -99999.0;
	x23 = y23 = z23 = -99999.0;

	/* roughly compute baseline */
	for (k = -ns2; k < ns + ns2; k++) {
		ts = t21 + k * dt;
		interpolate_SAT_orbit(&orb[ii], pt, p, pv, ts, &xs, &ys, &zs, &ir);

		ds = find_dist(xs, ys, zs, x11, y11, z11);
		if (b1 < 0.0 || ds < b1)
			endpoint_distance(k, ds, xs, ys, zs, &b1, &x21, &y21, &z21, &m1);

		ds = find_dist(xs, ys, zs, x12, y12, z12);
		if (b2 < 0.0 || ds < b2)
			endpoint_distance(k, ds, xs, ys, zs, &b2, &x22, &y22, &z22, &m2);

		ds = find_dist(xs, ys, zs, x13, y13, z13);
		if (b3 < 0.0 || ds < b3)
			endpoint_distance(k, ds, xs, ys, zs, &b3, &x23, &y23, &z23, &m3);
	}

	// refine the baseline computation with polynomial fit

	dstart = fabs(t11 - t21);
	// check whether it's the computation of same orbit.
	if (dstart > 10.) {
		poly_interp(&orb[ii], &ts, &b1, t21, x11, y11, z11, m1 * dt);
		interpolate_SAT_orbit_slow(&orb[ii], ts, &xs, &ys, &zs, &ir);
		x21 = xs;
		y21 = ys;
		z21 = zs;
		r[ii].B_offset_start = (ts - t21) * r[ii].vel;
		// seems this does not make much difference
		// if (r[0].SC_identity == 10) {
		/* approximate a secondary shift along orbit for alpha and baseline change
		 */
		// poly_interp(&orb[0],&ts,&b_tmp,t11,x21,y21,z21,m1*dt);
		// d1 = + (t11-ts)*r[0].vel;
		// r[ii].B_offset_start = r[ii].B_offset_start + d1;
		//}
		/* compute more orbital information at the min baseline based on m1 */
		// calc_height_velocity(&orb[ii], &r[ii], ts, ts, &height, &re_c, &vg,
		// &vtot, &rdot);

		poly_interp(&orb[ii], &ts, &b2, t21, x12, y12, z12, m2 * dt);
		interpolate_SAT_orbit_slow(&orb[ii], ts, &xs, &ys, &zs, &ir);
		x22 = xs;
		y22 = ys;
		z22 = zs;
		r[ii].B_offset_end = (ts - t22) * r[ii].vel;
		// if (r[0].SC_identity == 10) {
		// poly_interp(&orb[0],&ts,&b_tmp,t11,x22,y22,z22,m2*dt);
		// d2 = (t12-ts)*r[0].vel;
		// r[ii].B_offset_end = r[ii].B_offset_end + d2;
		//}

		poly_interp(&orb[ii], &ts, &b3, t21, x13, y13, z13, m3 * dt);
		interpolate_SAT_orbit_slow(&orb[ii], ts, &xs, &ys, &zs, &ir);
		x23 = xs;
		y23 = ys;
		z23 = zs;
		r[ii].B_offset_center = (ts - t23) * r[ii].vel;
		// if (r[0].SC_identity == 10) {
		// poly_interp(&orb[0],&ts,&b_tmp,t11,x23,y23,z23,m3*dt);
		// d3 = (t13-ts)*r[0].vel;
		// r[ii].B_offset_center = r[ii].B_offset_center + d3;
		//}
	}

	/* Not sure what these codes are used for
	// change back the dt settings
	dt = 0.5/r[0].prf;
	ns = (int) ((t12 - t11)/dt);    // seconds of frame
	dt = (t12 - t11)/(ns - 1);

	// set the extension to 50% unless ERS or Envisat and then set to 200%
	ns2 = ns*0.5;
	if(r[0].SC_identity < 3 || r[0].SC_identity == 4 || r[0].SC_identity == 10)
	ns2=ns*2;
  */

	/* fd_orbit = -2.0*rdot/r[0].lambda; */

	/* shouldn't happen ..					*/
	if (x21 == -99999.0)
		die("x11 not initialized", "");
	if (x22 == -99999.0)
		die("x12 not initialized", "");
	if (m1 == -99999)
		die("m1 not initialized", "");

	/* put in structure 					*/
	r[ii].baseline_start = b1;
	r[ii].baseline_end = b2;
	r[ii].baseline_center = b3;

	/* compute unit vectors 				*/
	find_unit_vectors(x11, y11, z11, &ru1, &xu1, &yu1, &zu1);
	find_unit_vectors(x12, y12, z12, &ru2, &xu2, &yu2, &zu2);
	find_unit_vectors(x13, y13, z13, &ru3, &xu3, &yu3, &zu3);

	/* compute sign of horizontal baseline 			*/
	get_sign(r[ii], x11, y11, x21, y21, &sign1);
	get_sign(r[ii], x12, y12, x22, y22, &sign2);
	get_sign(r[ii], x13, y13, x23, y23, &sign3);

	/* compute baseline components (horizontal and vertical) */
	bv1 = (x21 - x11) * xu1 + (y21 - y11) * yu1 + (z21 - z11) * zu1;
	bh1 = sign1 * sqrt(b1 * b1 - bv1 * bv1); //+ d1*d1);

	bv2 = (x22 - x12) * xu2 + (y22 - y12) * yu2 + (z22 - z12) * zu2;
	bh2 = sign2 * sqrt(b2 * b2 - bv2 * bv2); // + d2*d2);

	bv3 = (x23 - x13) * xu3 + (y23 - y13) * yu3 + (z23 - z13) * zu3;
	bh3 = sign3 * sqrt(b3 * b3 - bv3 * bv3); // + d3*d3);

	/* angle from horizontal 				*/
	r[ii].alpha_start = find_alpha_degrees(bv1, bh1);
	r[ii].alpha_end = find_alpha_degrees(bv2, bh2);
	r[ii].alpha_center = find_alpha_degrees(bv3, bh3);

	/* calculate parallel baseline 				*/
	find_parallel_perp_baseline(r[0], r[ii], dr, &bpara, &bperp);
	r[ii].bpara = bpara;
	r[ii].bperp = bperp;

	/* find expected offset in pixels (rshift and yshift) 	*/
	/*
	r[ii].ashift = -1*m1;
	r[ii].rshift = -1*(int) (r[ii].bpara/dr + (r[ii].near_range -
	r[0].near_range)/dr); fprintf(stderr,"ashift     =  %d\nrshift    =
	%d\n",r[ii].ashift,r[ii].rshift);
	*/
	/* a more accurate way to estimate offset in pixels    */

	/* since t11 point is out of scene, I need add a fraction of orbit time to it
	 */

	pt0 = (24.0 * 60.0 * 60.0) * orb[0].id + orb[0].sec;
	for (k = 0; k < nd; k++)
		pt[k] = pt0 + k * orb[0].dsec;
	t_sta = t11 + 2.0;
	interpolate_SAT_orbit(&orb[0], pt, p, pv, t_sta, &x_sta, &y_sta, &z_sta, &ir);

	re = r[0].RE;
	far_range = r[0].near_range + dr * r[0].num_rng_bins;
	rho = (r[0].near_range + far_range) / 2.0;
	b = r[0].ht + r[0].RE;
	theta = acos((b * b + rho * rho - re * re) / 2 / b / rho);

	/* find a unit look vector point to middle range of first line of the scene */
	/* in radar coordinate  				*/
	radar_look[0] = 0;
	radar_look[1] = cos(theta);
	radar_look[2] = -sin(theta);
	if (strncmp(r[0].lookdir, "L", 1) == 0)
		radar_look[2] = sin(theta);

	/* form a rotation matrix to convert the look vector from radar
	               coordinate to global cartesian coordinate           */

	find_unit_vectors(x_sta, y_sta, z_sta, &ru33, &xu33, &yu33, &zu33);

	o1y[0] = -xu33;
	o1y[1] = -yu33;
	o1y[2] = -zu33;

	/* get the velocity of the satellite   			*/
	pt0 = (24.0 * 60.0 * 60.0) * orb[0].id + orb[0].sec;
	for (k = 0; k < nd; k++)
		pt[k] = pt0 + k * orb[0].dsec;
	t_vel = t11 + 2.1;
	interpolate_SAT_orbit(&orb[0], pt, p, pv, t_vel, &x_vel, &y_vel, &z_vel, &ir);

	find_unit_vectors(x_vel - x_sta, y_vel - y_sta, z_vel - z_sta, &ru4, &xu4, &yu4, &zu4);
	o1x[0] = xu4;
	o1x[1] = yu4;
	o1x[2] = zu4;

	cross3(o1x, o1y, o1z);

	/* look vector in global cartesian cooridate is     	*/
	glob_look[0] = radar_look[0] * o1x[0] + radar_look[1] * o1y[0] + radar_look[2] * o1z[0];
	glob_look[1] = radar_look[0] * o1x[1] + radar_look[1] * o1y[1] + radar_look[2] * o1z[1];
	glob_look[2] = radar_look[0] * o1x[2] + radar_look[1] * o1y[2] + radar_look[2] * o1z[2];

	/* shoot the vector on the ground and get the tie point */

	target[0] = x_sta + glob_look[0] * rho;
	target[1] = y_sta + glob_look[1] * rho;
	target[2] = z_sta + glob_look[2] * rho;

	fll = (r[0].ra - r[0].rc) / r[0].ra;
	xyz2plh(target, target_llt, r[0].ra, fll);

	printf("lon_tie_point =  %f\n", (target_llt[1] > 180.0) ? target_llt[1] - 360.0 : target_llt[1]);
	printf("lat_tie_point =  %f\n", target_llt[0]);

	llt2rat_sub(filename[0], target_llt, target_rat_ref);
	llt2rat_sub(filename[ii], target_llt, target_rat_rep);

	/* find expected offset in pixels (rshift and yshift)   */
	r[ii].ashift = target_rat_rep[1] - target_rat_ref[1];
	r[ii].rshift = target_rat_rep[0] - target_rat_ref[0];

	/* write out prm format 				*/
	if (input_flag == 0)
		write_prm_baseline(r[ii]);

	/* write out in x,y format bperp baseline only 		*/
	if (input_flag == 1)
		write_bperp(r[ii], filename[ii]);
}
/*---------------------------------------------------------------------------*/
void read_input_file(char *inputfilename, int nfiles, char **filename) {
	int i;
	FILE *inputfile;

	if ((inputfile = fopen(inputfilename, "r")) == NULL)
		die("Can't open ", inputfilename);

	for (i = 0; i < nfiles; i++)
		fscanf(inputfile, " %s ", filename[i]);

	fclose(inputfile);
}
/*---------------------------------------------------------------------------*/
void baseline_parse_command_line(char **argv, int *nfiles, int *input_flag) {
	char tmp[132];
	FILE *inputfile;

	*nfiles = 0;

	if (strncmp(argv[1], "-input", 6) != 0) {

		printf("using command line\n");
		*nfiles = 2;
	}
	else {

		printf("using input file \n");
		*input_flag = 1;

		if ((inputfile = fopen(argv[2], "r")) == NULL)
			die("Can't open ", argv[2]);

		/* count the number of files in inputfile */
		while (fscanf(inputfile, " %s ", tmp) != EOF)
			(*nfiles)++;

		fclose(inputfile);
	}
}
/*---------------------------------------------------------------------------*/
void write_bperp(struct PRM r1, char *string) {
	FILE *GMTscript;
	float year, dyear;

	GMTscript = fopen("GMTscript", "w");

	year = floor(r1.SC_clock_start / 1000.0);
	dyear = (r1.SC_clock_start - 1000.0 * year) / 365.0;
	year = year + dyear;

	/* write date, baseline */
	fprintf(stdout, " %f %f %s\n", year, r1.bperp, string);

	/*
	        fprintf(GMTscript,"psxy Bperp.dat -R%f/%f/%f/%f -JX8 -P >
	   Bperp.ps\n");
	*/

	fclose(GMTscript);
}
/*---------------------------------------------------------------------------*/
void get_sign(struct PRM r, double x11, double y11, double x21, double y21, int *sign) {
	double rlnref, rlnrep;

	rlnref = atan2(y11, x11);
	rlnrep = atan2(y21, x21);

	*sign = 1.0;

	if (strncmp(r.orbdir, "D", 1) == 0)
		*sign = -1 * (*sign);
	if (strncmp(r.lookdir, "L", 1) == 0)
		*sign = -1 * (*sign);

	if (rlnrep < rlnref)
		*sign = -1 * (*sign);
}
/*---------------------------------------------------------------------------*/
void write_prm_baseline(struct PRM rep) {
	printf("SC_vel              = %.12f \n", rep.vel);
	printf("SC_height           = %.12f \n", rep.ht);
	printf("SC_height_start     = %.12f \n", rep.ht_start);
	printf("SC_height_end       = %.12f \n", rep.ht_end);
	printf("earth_radius        = %.12f \n", rep.RE);
	printf("rshift              = %d \n", rep.rshift);
	printf("sub_int_r           = 0.0 \n");
	printf("ashift              = %d\n", rep.ashift);
	printf("sub_int_a           = 0.0 \n");
	printf("B_parallel          = %.12f \n", rep.bpara);
	printf("B_perpendicular     = %.12f \n", rep.bperp);
	printf("baseline_start      = %.12f \n", rep.baseline_start);
	printf("baseline_center     = %.12f \n", rep.baseline_center);
	printf("baseline_end        = %.12f \n", rep.baseline_end);
	printf("alpha_start         = %.12f \n", rep.alpha_start);
	printf("alpha_center        = %.12f \n", rep.alpha_center);
	printf("alpha_end           = %.12f \n", rep.alpha_end);
	printf("B_offset_start      = %.12f \n", rep.B_offset_start);
	printf("B_offset_center     = %.12f \n", rep.B_offset_center);
	printf("B_offset_end        = %.12f \n", rep.B_offset_end);
}

/*---------------------------------------------------------------------------*/
void find_parallel_perp_baseline(struct PRM ref, struct PRM rep, double dr, double *bpara, double *bperp) {
	double rc, ra, rad, rlook, far_range, arg1, arg2;

	rad = M_PI / 180.0;

	/* this is not ellipsoid axis */

	/* use reference ?? */
	rc = ref.RE + ref.ht;
	ra = ref.RE;

	far_range = rep.near_range + dr * rep.num_rng_bins;

	/* calculate the look angle 1/2 way between the near and far range */

	arg1 = (rep.near_range * rep.near_range + rc * rc - ra * ra) / (2. * rep.near_range * rc);
	arg2 = (far_range * far_range + rc * rc - ra * ra) / (2. * far_range * rc);
	rlook = acos((arg1 + arg2) / 2.0);

	/* add the incidence angle correction to get the incidence angle  */

	arg1 = (-rep.near_range * rep.near_range + rc * rc + ra * ra) / (2. * ra * rc);
	arg2 = (-far_range * far_range + rc * rc + ra * ra) / (2. * ra * rc);
	rlook = rlook + acos((arg1 + arg2) / 2.0);

	*bpara = rep.baseline_start * sin(rlook - rep.alpha_start * rad);
	*bperp = rep.baseline_start * cos(rlook - rep.alpha_start * rad);
}
/*---------------------------------------------------------------------------*/
double find_alpha_degrees(double bv, double bh) {
	double a, rad;

	rad = M_PI / 180.0;

	a = atan2(bv, bh);
	a = a / rad;

	return (a);
}
/*---------------------------------------------------------------------------*/
double find_dist(double xs, double ys, double zs, double x, double y, double z) {
	double ds;
	double dx, dy, dz;

	dx = xs - x;
	dy = ys - y;
	dz = zs - z;

	ds = sqrt(dx * dx + dy * dy + dz * dz);

	return (ds);
}
/*---------------------------------------------------------------------------*/
void endpoint_distance(int k, double ds, double xs, double ys, double zs, double *b, double *x, double *y, double *z, int *m) {
	*b = ds;
	*x = xs;
	*y = ys;
	*z = zs;
	*m = k;
}
/*---------------------------------------------------------------------------*/
void find_unit_vectors(double x, double y, double z, double *ru, double *xu, double *yu, double *zu) {
	*ru = sqrt(x * x + y * y + z * z);
	*xu = x / (*ru);
	*yu = y / (*ru);
	*zu = z / (*ru);
}
/*---------------------------------------------------------------------------*/
