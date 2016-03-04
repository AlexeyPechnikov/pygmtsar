/*	$Id: xcorr.c 73 2013-04-19 17:59:45Z pwessel $	*/
/***************************************************************************/
/* xcorr does a 2-D cross correlation on complex or real images            */
/* either using a time convolution or wavenumber multiplication.           */
/*                                                                         */
/***************************************************************************/

/***************************************************************************
 * Creator:  Rob J. Mellors                                                *
 *           (San Diego State University)                                  *
 * Date   :  November 7, 2009                                              *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * DATE				        	                           *
 *							                   *
 * 011810       Testing and very monor cosmetic modfifications DTS         *
 * 061520	Problem with sub-pixel interpolation RJM		   *
 * 		- fixed bug in 2D interpolationa			   *
 * 		- revised read_xcorr_data to read in all x position	   *
 * 		- reads directly into float rather than int	           *
 * 		- add range interpolation	                           *
 * 		- eliminated obsolete options and code	                   *
 * 		- renamed xcorr_utils.c print_results.c	                   *
 * 		- further testing....					   *
 ***************************************************************************/

/*-------------------------------------------------------*/
#include "gmtsar.h"

char	*USAGE = "xcorr [GMT5SAR] - Compute 2-D cross-correlation of two images\n\n"
"\nUsage: xcorr master.PRM slave.PRM [-time] [-real] [-freq] [-nx n] [-ny n] [-xsearch xs] [-ysearch ys]\n"
"master.PRM     	PRM file for reference image\n"
"slave.PRM     	 	PRM file of secondary image\n"
"-time      		use time cross-correlation\n"
"-freq      		use frequency cross-correlation (default)\n"
"-real      		read float numbers instead of complex numbers\n"
"-noshift  		ignore ashift and rshift in prm file (set to 0)\n"
"-nx  nx    		number of locations in x (range) direction (int)\n"
"-ny  ny    		number of locations in y (azimuth) direction (int)\n"
"-nointerp     		do not interpolate correlation function\n"
"-range_interp ri  	interpolate range by ri (power of two) [default: 2]\n"
"-norange     		do not range interpolate \n"
"-xsearch xs		search window size in x (range) direction (int power of 2 [32 64 128 256])\n"
"-ysearch ys		search window size in y (azimuth) direction (int power of 2 [32 64 128 256])\n"
"-interp  factor    	interpolate correlation function by factor (int) [default, 16]\n"
"-v			verbose\n"
"output: \n freq_xcorr.dat (default) \n time_xcorr.dat (if -time option))\n"
"\nuse fitoffset.csh to convert output to PRM format\n"
"\nExample:\n"
"xcorr IMG-HH-ALPSRP075880660-H1.0__A.PRM IMG-HH-ALPSRP129560660-H1.0__A.PRM -nx 20 -ny 50 \n";

/*-------------------------------------------------------------------------------*/
int do_range_interpolate(void *API, struct FCOMPLEX *c, int nx, int ri, struct FCOMPLEX *work)
{
int	i;

	/* interpolate c and put into work */
	fft_interpolate_1d(API, c, nx, work, ri);

	/* replace original with interpolated (only half) */
	for (i=0; i<nx; i++) {
		c[i].r = work[i+nx/2].r;
		c[i].i = work[i+nx/2].i;
		}	

	return(EXIT_SUCCESS);
}
/*-------------------------------------------------------------------------------*/
/* complex arrays used in fft correlation					*/
/* load complex arrays and mask out slave					*/
/* c1 is master									*/
/* c2 is slave									*/
/* c3 used in fft complex correlation						*/
/* c1, c2, and c3 are npy by npx						*/
/* d1, d2 are npy by nx (length of line in SLC)					*/
/*-------------------------------------------------------------------------------*/
void assign_values(void *API, struct xcorr xc, int iloc)
{
int	i, j, k, sx, mx;
double	mean1, mean2;

	/* master and slave x offsets */
	mx = xc.loc[iloc].x - xc.npx/2;
	sx = xc.loc[iloc].x + xc.x_offset - xc.npx/2;

	for (i=0; i<xc.npy; i++){
		for (j=0; j<xc.npx; j++){
			k = i*xc.npx + j;

			xc.c3[k].i = xc.c3[k].r = 0.0f;

			xc.c1[k].r = xc.d1[i*xc.m_nx + mx + j].r;
			xc.c1[k].i = xc.d1[i*xc.m_nx + mx + j].i;

			xc.c2[k].r = xc.d2[i*xc.s_nx + sx + j].r;
			xc.c2[k].i = xc.d2[i*xc.s_nx + sx + j].i;
			}
		}

	/* range interpolate */
	if (xc.ri > 1) {
		for (i=0; i<xc.npy; i++) {
			do_range_interpolate(API, &xc.c1[i*xc.npx], xc.npx, xc.ri, xc.ritmp);
			do_range_interpolate(API, &xc.c2[i*xc.npx], xc.npx, xc.ri, xc.ritmp);
			}
		}

	/* convert to amplitude and demean */
	mean1 = mean2 = 0.0;
	for (i=0; i<xc.npy*xc.npx; i++) {
		xc.c1[i].r = Cabs(xc.c1[i]);
		xc.c1[i].i = 0.0f;

		xc.c2[i].r = Cabs(xc.c2[i]);
		xc.c2[i].i = 0.0f;

		mean1 += xc.c1[i].r;
		mean2 += xc.c2[i].r;
		}

	mean1 /= (double) (xc.npy * xc.npx);
	mean2 /= (double) (xc.npy * xc.npx);

	for (i=0; i<xc.npy*xc.npx; i++) {
		xc.c1[i].r = xc.c1[i].r - (float)mean1;
		xc.c2[i].r = xc.c2[i].r - (float)mean2;
		}

	/* apply mask */
	for (i=0; i<xc.npy*xc.npx; i++) {
			xc.c1[i].i = xc.c2[i].i = 0.0f;
			xc.c2[i].r = xc.c2[i].r * (float) xc.mask[i];

			xc.i1[i] = (int) (xc.c1[i].r);
			xc.i2[i] = (int) (xc.c2[i].r);
		}

	if (debug) fprintf(stderr," mean %lf\n", mean1);
	if (debug) fprintf(stderr," mean %lf\n", mean2);
}
/*-------------------------------------------------------------------------------*/
void do_correlation(void *API, struct xcorr xc)
{
	int	i, j, iloc, istep;

 	/* opportunity for multiple processors */
	istep = 1;

	/* allocate arrays   			*/
	allocate_arrays(&xc);

	/* make mask 				*/
	make_mask(xc);
	
	iloc = 0;
	for (i=0; i<xc.nyl; i+=istep){

		/* read in data for each row */
		read_xcorr_data(xc, iloc);

		for (j=0; j<xc.nxl; j++){

			if (debug) fprintf(stderr," initial: iloc %d (%d,%d)\n",iloc, xc.loc[iloc].x,xc.loc[iloc].y);

			/* copy values from d1,d2 (real) to c1,c2 (complex) */
			assign_values(API, xc, iloc);

			if (debug) print_complex(xc.c1, xc.npy, xc.npx, 1);
			if (debug) print_complex(xc.c2, xc.npy, xc.npx, 1);

			/* correlate patch with data over offsets in time domain */
			if (xc.corr_flag < 2) do_time_corr(xc, iloc); 

			/* correlate patch with data over offsets in freq domain */
			if (xc.corr_flag == 2) do_freq_corr(API, xc, iloc); 

			/* oversample correlation surface  to obtain sub-pixel resolution */
			if (xc.interp_flag == 1) do_highres_corr(API, xc, iloc); 

			/* write out results */
			print_results(xc, iloc);

			iloc++;
		} /* end of x iloc loop */
	} /* end of y iloc loop */
}
/*-------------------------------------------------------------------------------*/
/* want to avoid circular correlation so mask out most of b			*/
/* could adjust shape for different geometries					*/
/*-------------------------------------------------------------------------------*/
void make_mask(struct xcorr xc)
{
int i,j,imask;
imask = 0;

	for (i=0; i<xc.npy; i++){
		for (j=0; j<xc.npx; j++){
			 xc.mask[i*xc.npx + j] = 1;
			if ((i < xc.ysearch) || (i >= (xc.npy - xc.ysearch))) {
				xc.mask[i*xc.npx +j] = imask;
				}
			if ((j < xc.xsearch) || (j >= (xc.npx - xc.xsearch))) {
				xc.mask[i*xc.npx + j] = imask;
				}
			}
		}
}
/*-------------------------------------------------------------------------------*/
void allocate_arrays(struct xcorr *xc)
{
int	nx, ny, nx_exp, ny_exp;

	xc->d1 = (struct FCOMPLEX *) malloc(xc->m_nx*xc->npy*sizeof(struct FCOMPLEX));
	xc->d2 = (struct FCOMPLEX *) malloc(xc->s_nx*xc->npy*sizeof(struct FCOMPLEX));

	xc->i1 = (int *) malloc(xc->npx*xc->npy*sizeof(int));
	xc->i2 = (int *) malloc(xc->npx*xc->npy*sizeof(int));

	xc->c1 = (struct FCOMPLEX *) malloc(xc->npx*xc->npy*sizeof(struct FCOMPLEX));
	xc->c2 = (struct FCOMPLEX *) malloc(xc->npx*xc->npy*sizeof(struct FCOMPLEX));
	xc->c3 = (struct FCOMPLEX *) malloc(xc->npx*xc->npy*sizeof(struct FCOMPLEX));

	xc->ritmp = (struct FCOMPLEX *) malloc(xc->ri*xc->npx*sizeof(struct FCOMPLEX));
	xc->mask = (short *) malloc(xc->npx*xc->npy*sizeof(short));

	/* this is size of correlation patch */
	xc->corr = (double *) malloc(2*xc->ri*(xc->nxc)*(xc->nyc)*sizeof(double));

	if (xc->interp_flag == 1){
		nx = 2*xc->n2x;
		ny = 2*xc->n2y;
		nx_exp = nx * (xc->interp_factor);
		ny_exp = ny * (xc->interp_factor);
		xc->md = (struct FCOMPLEX *) malloc(nx * ny * sizeof(struct FCOMPLEX));
		xc->cd_exp = (struct FCOMPLEX *) malloc(nx_exp * ny_exp * sizeof(struct FCOMPLEX));
		}
}

/*-------------------------------------------------------*/
int main(int argc,char **argv)
{
	int	input_flag, nfiles;
	struct	xcorr xc;
	clock_t start, end;
	double	cpu_time;
	void	*API = NULL; /* GMT API control structure */

	verbose = 0;
	debug = 0;
	input_flag = 0;
	nfiles = 2;
	xc.interp_flag = 0;
	xc.corr_flag = 2;

	/* Begin: Initializing new GMT5 session */
	if ((API = GMT_Create_Session (argv[0], 0U, 0U, NULL)) == NULL) return EXIT_FAILURE;

	if (argc < 3) die (USAGE, "");

	set_defaults(&xc);
	
	parse_command_line(argc, argv, &xc, &nfiles, &input_flag, USAGE);

	/* read prm files */
	if (input_flag == 0) handle_prm(argv, &xc, nfiles);

	if (debug) print_params(&xc);

	/* output file */
	if (xc.corr_flag == 0) strcpy(xc.filename,"time_xcorr.dat");
	if (xc.corr_flag == 1) strcpy(xc.filename,"time_xcorr_Gatelli.dat");
	if (xc.corr_flag == 2) strcpy(xc.filename,"freq_xcorr.dat");

	xc.file = fopen(xc.filename,"w");
	if (xc.file == NULL) die("Can't open output file",xc.filename); 

	/* x locations, y locations */
	get_locations(&xc);

	/* calculate correlation at all points */
	start = clock();

	do_correlation (API, xc);

	/* write the a_stretch_a based on the PRF differences */
/*
        fprintf(xc.file,"a_stretch_a  =  %f \n",xc.astretcha);
*/

	end = clock();
	cpu_time  = ((double) (end-start)) / CLOCKS_PER_SEC;	
	fprintf(stdout, " elapsed time: %lf \n", cpu_time);

	fclose(xc.data1); 
	fclose(xc.data2);

	if (GMT_Destroy_Session (API)) return EXIT_FAILURE;	/* Remove the GMT machinery */

	return(EXIT_SUCCESS);
}
