/*	$Id: phasediff.c 109 2015-01-19 23:01:24Z sandwell $	*/
/***************************************************************************
 * phasediff reads two complex SAR_SLC images and computes the phase       *
 * difference of the images removing the effects of the curved earth       *
 * and optionally the topography. Orbit information is provided in the     *
 * second PRM-file.                                                        *
 * The baseline changes linearly between the start and end of the image.   *
 * This linear baseline model is appropriate for a single frame but higher *
 * order terms are needed for long swaths of data.                         *
 * Model phase can also be removed as an option.                           *
 * The images should be matched to a level that will                       *
 * produce interference fringes and must be the same size.                 *
 ***************************************************************************/

/***************************************************************************
 * Creator:  David T. Sandwell  and Evelyn J. Price                        *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  7/7/95                                                        *
 * Rewritten: Rob Mellors (San Diego State University)                     *
 * Date   :  11/7/09                                                       *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * DATE                                                                    *
 *                                                                         *
 * 11/18/96     Changed sign of earth flattening phase to fit with signs   *
 *              of baseline estimates.                                     *
 * 11/18/96     Changed to read GIPS headers to gather parameters.         *
 * 11/21/96     Changed to update baseline estimate.                       *
 * 06/19/97     Changed to NOT update baseline parameters                  *
 * 06/25/97     Changed to change baseline along the frame                 *
 * 03/10/99     Changed to remove reference phase due to topography        *
 * 04/06/01     Changed to remove model or other known phase               *
 * 06/21/02     Changed to scale the interferogram by 1/sqrt(amp)          *
 * 08/16/06     Changed to use the topo_ra to shift the range              *
 *              coordinate of the repeat image for long baselines.         *
 * 07/02/08     provide ability to shift topo_ra                           *
 * 07/11/08     Added spacecraft height start and end to make phase        *
 *              continuous across patches                                  *
 * 09/18/08     Added TEC_start and TEC_end to account for the ionosphere  *
 * 08/18/09     Added small correction to the phase for the elevation-     *
 *              dependent range shift, important for long baselines        *
 * 11/07/09     Code rewritten by RJM to read PRM-files and write grd-files*
 * 04/24/10     changed calc_phase subroutine based on pdiff.c             *
 *              changed the phase calculation part to match pdiff.c        *
 *              this change is to account for slight range change due to   *
 *              topography, which is accounted for in pdiff.c              *
 * 10/23/10     changed calc_phase to calc_drho and completely rewrote     *
 *              the topographic phase correction to use all the nonlinear  *
 *              terms.
 ***************************************************************************/
#include "gmtsar.h"

char *USAGE = "phasediff [GMTSAR] - Compute phase difference of two images\n\n"
              "\nUsage: "
              "phasediff ref.PRM rep.PRM [-topo topo_ra.grd] [-model "
              "modelphase.grd]\n (topo_ra and model in GMT grd format)\n";

/*--------------------------------------------------------------*/
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

/*--------------------------------------------------------------*/
void calc_average_topo(double *avet, int xdimt, int ydimt, float *topo) {
	double sumt;
	int k, nsum;

	sumt = 0.0;
	nsum = 0;

	/* compute the average topography, save the value, and remove it from the
	 * topography */
	for (k = 0; k < xdimt * ydimt; k++) {
		sumt += topo[k];
		nsum++;
	}
	*avet = sumt / nsum;
	for (k = 0; k < xdimt * ydimt; k++) {
		topo[k] = topo[k] - (float)*avet;
	}
	if (verbose)
		fprintf(stderr, " mean topo %lf\n", *avet);
}

void print_prm_params(struct PRM p1, struct PRM p2) {
	fprintf(stderr, " SLC 1: num_rng_bins %d num_lines %d \n", p1.num_rng_bins, p1.num_lines);
	fprintf(stderr, " SLC 2: num_rng_bins %d num_lines %d \n", p2.num_rng_bins, p2.num_lines);
	fprintf(stderr, " lambda %f \n", p2.lambda);
	fprintf(stderr, " baseline_start %f \n", p2.baseline_start);
	fprintf(stderr, " baseline_end %f \n", p2.baseline_end);
	fprintf(stderr, " B_offset_start %f \n", p2.B_offset_start);
	fprintf(stderr, " B_offset_end %f \n", p2.B_offset_end);
	fprintf(stderr, " alpha_start %f \n", p2.alpha_start);
	fprintf(stderr, " alpha_end %f \n", p2.alpha_end);
	fprintf(stderr, " near_range %f \n", p2.near_range);
	fprintf(stderr, " rng_samp_rate %f \n", p2.fs);
	fprintf(stderr, " sc_clock_start %f \n", p2.SC_clock_start);
	fprintf(stderr, " sc_clock_stop %f \n", p2.SC_clock_stop);
	fprintf(stderr, " clock_start %f \n", p2.clock_start);
	fprintf(stderr, " clock_stop %f \n", p2.clock_stop);
	fprintf(stderr, " prf %f \n", p2.prf);
}

void fix_prm_params(struct PRM *p, char *s) {
	double delr;

	delr = SOL / p->fs / 2.0;

	/* this is the correction for the range and azimuth shifts of the re-aligned
	 * SLC images */
	//        if(p->sub_int_r < 0.) p->sub_int_r = 0.;
	//        if(p->sub_int_a < 0.) p->sub_int_a = 0.;
	p->near_range = p->near_range + (p->st_rng_bin - p->chirp_ext + p->rshift + p->sub_int_r - 1) * delr;
	p->SC_clock_start = p->SC_clock_start + (p->ashift + p->sub_int_a) / (p->prf * 86400.0) +
	                    (p->nrows - p->num_valid_az) / (2.0 * p->prf * 86400);
	p->SC_clock_stop = p->SC_clock_start + (p->num_valid_az * p->num_patches) / (p->prf * 86400.0);
}

/*--------------------------------------------------------------*/
/* read topo_ra and model files if provided			*/
/* must be in GMT binary grd format					*/
void read_optional_args(void *API, int argc, char **argv, struct PRM *tp, int *topoflag, struct PRM *mp, int *modelflag) {
	int i;
	struct GMT_GRID *M = NULL, *T = NULL; /* Grid structures containing ->header and ->data */

	for (i = 3; i < argc; i++) {
		if (strcmp(argv[i], "-topo") == 0) {
			fprintf(stderr, "reading topo %s\n", argv[i + 1]);
			strcpy(tp->input_file, argv[i + 1]);
			if ((T = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, tp->input_file,
			                       NULL)) == NULL)
				die("cannot open topofile", tp->input_file);
			tp->num_rng_bins = T->header->n_columns;
			tp->num_lines = T->header->n_rows;
			*topoflag = 1;
			i++;
			GMT_Destroy_Data(API, &T);
		}

		if (strcmp(argv[i], "-model") == 0) {
			fprintf(stderr, "reading model %s\n", argv[i + 1]);
			strcpy(mp->input_file, argv[i + 1]);
			if ((M = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, mp->input_file,
			                       NULL)) == NULL)
				die("cannot open model file", mp->input_file);
			mp->num_rng_bins = M->header->n_columns;
			mp->num_lines = M->header->n_rows;
			*modelflag = 1;
			i++;
			GMT_Destroy_Data(API, &M);
		}
	}
}

int main(int argc, char **argv) {
	int j, k, istart;
	int topoflag, modelflag;
	int xdim = 0, ydim = 0; /* size of SLC file */
	int ydim_start;         /* start of SLC filesize */
	int xt, yt;             /* size of topo file, increment */
	int xm, ym;             /* size of model file, increment */
	uint64_t left_node;
	short *d1 = NULL, *d2 = NULL; /* pointers to input data files */
	double *xs = NULL, *shft = NULL, *ss = NULL, *as = NULL, *topo2 = NULL;
	double *real = NULL, *imag = NULL, *range = NULL, *drho = NULL, *drho0 = NULL;
	double drange, dt, tspan, time, time2;
	double ht0, htc, htf, dht, ddht, height;
	double alpha, cnst, pha, avet;
	double B, Bh, Bv, dBh, dBv, ddBh, ddBv;

	double Bx, dBx, ddBx, Bx0, Bxc, Bxf;

	double Bh0, Bhc, Bhf, Bv0, Bvc, Bvf;
	double ys, test, inc[2], wesn[4];
	double xdect, ydect, xdecm, ydecm, rdumt;
	FILE *SLCfile1 = NULL, *SLCfile2 = NULL;
	fcomplex *intfp = NULL, *iptr1 = NULL, *iptr2 = NULL, pshif;
	struct PRM p1, p2, tp, mp;
	void *API = NULL;                       /* GMT control structure */
	struct GMT_GRID *M = NULL, *T = NULL;   /* Grid structures containing ->header and ->data */
	struct GMT_GRID *RE = NULL, *IM = NULL; /* For the real and imaginary grids */

	double *range2 = NULL;

	/* Begin: Initializing new GMT session */
	if ((API = GMT_Create_Session(argv[0], 0U, 0U, NULL)) == NULL)
		return EXIT_FAILURE;

	verbose = 0;
	topoflag = modelflag = 0;
	xdect = ydect = 1.0;
	xdecm = ydecm = 1.0;
	avet = 0.0;
	ydim_start = 0;

	if (argc < 3)
		die(USAGE, "");

	/* read prm file into two pointers */
	get_prm(&p1, argv[1]);
	get_prm(&p2, argv[2]);

	if (verbose)
		fprintf(stderr, "near range: %lf %lf \n", p1.near_range, p2.near_range);

	if (argc > 3)
		read_optional_args(API, argc, argv, &tp, &topoflag, &mp, &modelflag);

	if (debug)
		print_prm_params(p1, p2);
	if (p2.baseline_start < -9000)
		die("baseline < -9000 not set ?", "");

	/* near_range, SC_clock_start, and SC_clock_stop need to be changed */
	fix_prm_params(&p1, argv[1]);
	fix_prm_params(&p2, argv[2]);

	if (verbose)
		fprintf(stderr, "near range: %lf %lf \n", p1.near_range, p2.near_range);

	/* open SLC files */
	if ((SLCfile1 = fopen(p1.SLC_file, "r")) == NULL)
		die("Can't open SLCfile", p1.SLC_file);
	if ((SLCfile2 = fopen(p2.SLC_file, "r")) == NULL)
		die("Can't open SLCfile", p2.SLC_file);

	/* set width and length */
	/* check dimensions of the two SLC files */
	if (p1.num_rng_bins == p2.num_rng_bins) {
		xdim = p1.num_rng_bins;
	}
	else {
		die("The dimensions of range do not match", "");
	}

	if (p1.num_patches * p1.num_valid_az == p2.num_patches * p2.num_valid_az) {
		ydim = p1.num_patches * p1.num_valid_az;
	}
	else {
		die("The dimensions of azimuth do not match", "");
	}
	fprintf(stderr, " xdim %d, ydim %d \n", xdim, ydim);

	/* set heights */
	htc = p1.ht;
	ht0 = p1.ht_start;
	htf = p1.ht_end;

	inc[GMT_X] = inc[GMT_Y] = 1.0; /* Pixels */
	wesn[GMT_XLO] = 0.0;
	wesn[GMT_XHI] = inc[GMT_X] * xdim;
	wesn[GMT_YLO] = 0.0;
	wesn[GMT_YHI] = inc[GMT_Y] * ydim;
	if ((RE = GMT_Create_Data(API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, wesn, inc, GMT_GRID_PIXEL_REG, 0, NULL)) ==
	    NULL)
		die("could not allocate grid header", "");
	if (GMT_Set_Comment(API, GMT_IS_GRID, GMT_COMMENT_IS_REMARK, "real", RE))
		return EXIT_FAILURE;
	if ((IM = GMT_Create_Data(API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, wesn, inc, GMT_GRID_PIXEL_REG, 0, NULL)) ==
	    NULL)
		die("could not allocate grid header", "");
	if (GMT_Set_Comment(API, GMT_IS_GRID, GMT_COMMENT_IS_REMARK, "imag", RE))
		return EXIT_FAILURE;

	/* allocate memory */
	drho = (double *)malloc(xdim * sizeof(double));
	drho0 = (double *)malloc(xdim * sizeof(double));
	range = (double *)malloc(xdim * sizeof(double));
	range2 = (double *)malloc(xdim * sizeof(double));

	intfp = (fcomplex *)malloc(xdim * sizeof(fcomplex));
	iptr1 = (fcomplex *)malloc(xdim * sizeof(fcomplex));
	iptr2 = (fcomplex *)malloc(xdim * sizeof(fcomplex));

	d1 = (short *)malloc(2 * xdim * sizeof(short));
	d2 = (short *)malloc(2 * xdim * sizeof(short));

	shft = (double *)malloc(xdim * sizeof(double));
	xs = (double *)malloc(xdim * sizeof(double));
	topo2 = (double *)malloc(xdim * sizeof(double));
	imag = (double *)malloc(xdim * sizeof(double));
	real = (double *)malloc(xdim * sizeof(double));
	ss = (double *)malloc(xdim * sizeof(double));
	as = (double *)malloc(xdim * sizeof(double));

	/* open and read topo file and allocate memory */
	if (topoflag) {
		if ((T = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, tp.input_file, NULL)) ==
		    NULL)
			return EXIT_FAILURE;
		if (xdim % T->header->n_columns != 0 || ydim % T->header->n_rows != 0)
			die("The dimension SLC must be multiplication factor of the topo_ra", tp.input_file);
		if (GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, tp.input_file, T) == NULL)
			return EXIT_FAILURE;
		rdumt = floor(T->header->z_max + 1.0);

		if (verbose)
			fprintf(stderr, "\n%s %s %d %d\n", T->header->title, tp.input_file, T->header->n_columns, T->header->n_rows);
		if (verbose)
			fprintf(stderr, "\n%f %f %f %f %f\n", T->header->wesn[GMT_XLO], T->header->wesn[GMT_YLO], T->header->inc[GMT_X],
			        T->header->inc[GMT_Y], rdumt);

		/* T->header->inc[GMT_X] or T->header->inc[GMT_Y] may be negative */
		xdect = fabs(T->header->inc[GMT_X]);
		ydect = fabs(T->header->inc[GMT_Y]);

		/* calculate the average and remove the average from the topography */

		calc_average_topo(&avet, T->header->n_columns, T->header->n_rows, T->data);
		if (verbose)
			fprintf(stderr, " read topo file: average %f \n", avet);
	}

	/* open and read the model file and allocate the memory */

	if (modelflag) {
		if ((M = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, mp.input_file, NULL)) ==
		    NULL)
			return EXIT_FAILURE;
		if (xdim % M->header->n_columns != 0 || ydim % M->header->n_rows != 0)
			die("The dimension SLC must be multiplication factor of the modelphase", mp.input_file);
		if (GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, mp.input_file, M) == NULL)
			return EXIT_FAILURE;
		xdecm = fabs(M->header->inc[GMT_X]);
		ydecm = fabs(M->header->inc[GMT_Y]);
	}

	/*   compute the time span and the time spacing */

	tspan = 86400. * fabs(p2.SC_clock_stop - p2.SC_clock_start);
	dt = tspan / (ydim - 1);
	if (tspan < 0.01 || p2.prf < 0.01)
		die("check sc_clock_start, _end, or prf", "");

	/* setup the default parameters */

	drange = SOL / (2.0 * p2.fs);
	alpha = p2.alpha_start * PI / 180.0;
	cnst = -4.0 * PI / p2.lambda;

	for (k = 0; k < xdim; k++) {
		// range[k]=p1.near_range+k*drange;
		range[k] = p1.near_range + k * (1 + p1.stretch_r) * drange;
		topo2[k] = 0.;
		xs[k] = k;
	}

	/* calculate initial baselines */
	/*
	if (p2.SC_identity == 100) {

	    double tmpB = 0.0,alpha = 0.0;
	    tmpB = p2.baseline_start;
	    alpha = p2.alpha_start;
	    p2.baseline_start = sqrt(pow((p1.ashift+p1.sub_int_a -
	p2.ashift-p2.sub_int_a)*p2.vel/p2.prf,2)+pow(tmpB,2)); p2.alpha_start =
	acos(tmpB*sin(alpha*PI/180.0)/p2.baseline_start)/PI*180.0;

	    tmpB = p2.baseline_end;
	    alpha = p2.alpha_end;
	    p2.baseline_end = sqrt(pow((p1.ashift+p1.sub_int_a -
	p2.ashift-p2.sub_int_a)*p2.vel/p2.prf,2)+pow(tmpB,2)); p2.alpha_end =
	acos(tmpB*sin(alpha*PI/180.0)/p2.baseline_end)/PI*180.0;

	    tmpB = p2.baseline_center;
	    alpha = p2.alpha_center;
	    p2.baseline_center = sqrt(pow((p1.ashift+p1.sub_int_a -
	p2.ashift-p2.sub_int_a)*p2.vel/p2.prf,2)+pow(tmpB,2)); p2.alpha_center =
	acos(tmpB*sin(alpha*PI/180.0)/p2.baseline_center)/PI*180.0;
	}
	*/

	Bh0 = p2.baseline_start * cos(p2.alpha_start * PI / 180.0);
	Bv0 = p2.baseline_start * sin(p2.alpha_start * PI / 180.0);
	Bhf = p2.baseline_end * cos(p2.alpha_end * PI / 180.0);
	Bvf = p2.baseline_end * sin(p2.alpha_end * PI / 180.0);
	Bx0 = p2.B_offset_start;
	Bxf = p2.B_offset_end;

	/* first case is quadratic baseline model, second case is default linear model
	 */

	if (p2.baseline_center != NULL_DOUBLE || p2.alpha_center != NULL_DOUBLE || p2.B_offset_center != NULL_DOUBLE) {

		Bhc = p2.baseline_center * cos(p2.alpha_center * PI / 180.0);
		Bvc = p2.baseline_center * sin(p2.alpha_center * PI / 180.0);
		Bxc = p2.B_offset_center;

		dBh = (-3. * Bh0 + 4 * Bhc - Bhf) / tspan;
		dBv = (-3. * Bv0 + 4 * Bvc - Bvf) / tspan;
		ddBh = (2. * Bh0 - 4 * Bhc + 2 * Bhf) / (tspan * tspan);
		ddBv = (2. * Bv0 - 4 * Bvc + 2 * Bvf) / (tspan * tspan);

		dBx = (-3. * Bx0 + 4 * Bxc - Bxf) / tspan;
		ddBx = (2. * Bx0 - 4 * Bxc + 2 * Bxf) / (tspan * tspan);
	}
	else {
		dBh = (Bhf - Bh0) / tspan;
		dBv = (Bvf - Bv0) / tspan;
		dBx = (Bxf - Bx0) / tspan;
		ddBh = ddBv = ddBx = 0.0;
	}

	/* calculate height increment */
	dht = (-3. * ht0 + 4 * htc - htf) / tspan;
	ddht = (2. * ht0 - 4 * htc + 2 * htf) / (tspan * tspan);

	/* revise params in accordance with first_line 	*/
	if ((p1.first_line > 0) && (p2.first_line > 0))
		ydim_start = p1.first_line - 1;

	/* now go through all the rows 		*/
	for (j = ydim_start; j < (ydim + ydim_start); j++) {

		for (k = 0; k < xdim; k++) {
			// range[k]=p1.near_range+k*drange;
			range2[k] = range[k] + j * p1.a_stretch_r * drange;
		}

		/* read data from complex i2 SLC 	*/
		read_SLC_short2float(SLCfile1, p1.SLC_file, d1, &iptr1[0], xdim, 1, DFACT);
		read_SLC_short2float(SLCfile2, p2.SLC_file, d2, &iptr2[0], xdim, 1, DFACT);

		yt = j / ydect; /* for topo_ra */
		ym = j / ydecm; /* for modelphase */

		/* calculate the change in baseline and height along the frame */
		time = j * dt;
		time2 = time * time;
		Bh = Bh0 + dBh * time + ddBh * time2;
		Bv = Bv0 + dBv * time + ddBv * time2;

		Bx = Bx0 + dBx * time + ddBx * time2;

		B = sqrt(Bh * Bh + Bv * Bv);
		alpha = atan2(Bv, Bh);
		height = ht0 + dht * time + ddht * time2;

		for (k = 0; k < xdim; k++) {
			shft[k] = 0.;
			if (topoflag) {
				xt = k / xdect;
				topo2[k] = T->data[xt + T->header->n_columns * yt];
			}
		}

		/* calculate the combined earth curvature and topography correction */
		calc_drho(xdim, range2, topo2, avet, p1.RE, height, B, alpha, Bx, drho);
//        for (k = 0; k < xdim; k++) range2[k] = 0.0;
//        calc_drho(xdim, range2, topo2, avet, p1.RE, height, B, alpha, Bx, drho0);

		// if (j == 50) printf("drho = %.12f\n",drho[50]);

		/* loop over range to make topographic and model phase corrections */
		for (k = 0; k < xdim; k++) {
			//
			// iptr1[k].r = 1.0;
			// iptr1[k].i = 1.0;
			//

			intfp[k] = iptr1[k];
			//pha = cnst * (drho[k] - drho0[k]);
			pha = cnst * drho[k];

			if (modelflag) {
				xm = k / xdecm; /* xm is increment for model phase. note they are all
				                   integers */
				pha = pha - M->data[xm + M->header->n_columns * ym];
			}
			pshif = Cexp(pha);
			//intfp[k] = Cmul(intfp[k], pshif);
            intfp[k].r = pha;
            intfp[k].i = 0.0;
		}

		/* shift the range of the repeat image to improve image matching for very
		 * long baselines > 1000 m */
		if ((topoflag > 0) && (p2.baseline_start > 1000.0)) {

			/* compute the range change with no topography so the range shift can be
			 * determined for the spline */
			calc_drho(xdim, range, shft, avet, p1.RE, height, B, alpha, Bx, drho0);

			for (k = 0; k < xdim; k++) {
				shft[k] = (drho0[k] - drho[k]) / drange;
				real[k] = intfp[k].r;
				imag[k] = intfp[k].i;
			}
			/* shift the real part */
			spline_(&istart, &xdim, xs, real, ss, as);
			for (k = 0; k < xdim; k++) {
				ys = xs[k] + shft[k];
				evals_(&istart, &ys, &xdim, xs, real, ss, &test);
				intfp[k].r = (float)test;
			}
			/* shift imaginary part  */
			spline_(&istart, &xdim, xs, imag, ss, as);
			for (k = 0; k < xdim; k++) {
				ys = xs[k] + shft[k];
				evals_(&istart, &ys, &xdim, xs, real, ss, &test);
				intfp[k].i = (float)test;
			}
		}

		/* make interferogram */
		left_node = GMT_Get_Index(API, RE->header, j, 0);

		for (k = 0; k < xdim; k++) {

			
			//iptr2[k].r = 1.0;
			//iptr2[k].i = 1.0;
			

			//iptr2[k] = Conjg(iptr2[k]);
			//intfp[k] = Cmul(intfp[k], iptr2[k]);
			RE->data[left_node + k] = intfp[k].r;
			IM->data[left_node + k] = intfp[k].i;
		}
	}

	if (GMT_Write_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, "real.grd=bf", RE)) {
		die("Failed to update real.grd grid header", "");
	}
//
//	if (GMT_Write_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, "imag.grd=bf", IM)) {
//		die("Failed to update imag.grd grid header", "");
//	}

	if (GMT_Destroy_Session(API))
		return EXIT_FAILURE; /* Remove the GMT machinery */

	return (EXIT_SUCCESS);
}
