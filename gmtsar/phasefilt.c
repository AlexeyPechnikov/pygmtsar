/*	$Id: phasefilt.c 80 2014-08-11 17:04:28Z sandwell $	*/
#include "gmtsar.h"

/*-------------------------------------------------------------------------------------	*/
/*
 *	implementation of adaptive non-linear phase filter based on:				
 * 	Goldstein, R. M. and C. L. Werner, 1998, Radar interferogram filtering 
 *	for geophysical applications, (25),21, 4035-4038.
 *
 *	Baran, I., M. P. Stewart, B. M. Kampes, Z. Perski, and P. Lilly, 2003
 *	A modification of the Golstein radar interferogram filter,
 *	IEE trans. geoscience and remote sensing, (41), 9, 2114-2118.
 *
 *	rjm, SDSU, may 2010
 *	comments: 
 *	- handled edges in a kludgey way
 *	- added memory allocation checks
 *
 *	future:
 *	- added option to calculate power spectra with a taper
 *	- (synthetics suggested a slight improvement from reduced sidelobes)
 *	
 *      Dec 2010 - bug in weighting fixed [apply to complex]
 *     	default psize changed to 32 
 *      '-complex_out' option added to write out filtered real and imag
 *
 *-------------------------------------------------------------------------------------	*/
int verbose;
int debug;

int print_cpatch(int nx, int ny, struct FCOMPLEX *m)
{
int i;
	for (i=0; i<nx*ny; i++){
		fprintf(stderr," (%4.2f %4.2f) ", m[i].r, m[i].i);
		if ((i+1)/nx  == (i+1)/(float) nx) fprintf(stderr,"\n");
		}
	return(EXIT_SUCCESS);
}

int print_patch(int nx, int ny, float *m)
{
int i;
	for (i=0; i<nx*ny; i++){
		fprintf(stderr," %4.2f", m[i]);
		if ((i+1)/nx  == (i+1)/(float) nx) fprintf(stderr,"\n");
		}
	return(EXIT_SUCCESS);
}

int write_grdfile (void *API, struct GMT_GRID *T, char *fname, char *prog, char *type, float *data, int verbose)
{
	if (verbose) fprintf (stderr," writing %s \n", fname);

	strcpy (T->header->command, prog);
	strcpy (T->header->remark, type);
	T->data = data;
	if (GMT_Write_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, fname, T)) die("cannot create ",fname);
	T->data = NULL;
	return(EXIT_SUCCESS);
}

int read_file_hdr (void *API, char *sim, struct GMT_GRID **IM, char *sre, struct GMT_GRID **RE, unsigned int *xdim, unsigned int *ydim)
{
	struct GMT_GRID *I = NULL, *R = NULL;
	if ((I = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, sim, NULL)) == NULL) die("error reading file",sim);
	if ((R = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, sre, NULL)) == NULL) die("error reading file",sre);

	if (!(R->header->nx == I->header->nx && R->header->ny == I->header->ny)) die("dimensions not equal!","");

	*xdim = I->header->nx;
	*ydim = I->header->ny;
	*RE = R;	*IM = I;
	return (EXIT_SUCCESS);
}

int phasefilt_read_data (void *API, char *imname, struct GMT_GRID *IM, char *rename, struct GMT_GRID *RE, struct FCOMPLEX *cdata, float *amp)
{
	long	n, i;

	n = IM->header->nx * IM->header->ny;

	if (GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, imname, IM) == NULL) die("error reading file",imname);
	if (GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, imname, RE) == NULL) die("error reading file",rename);

	for (i=0; i<n; i++){
		cdata[i].r = RE->data[i];
		cdata[i].i = IM->data[i];
		amp[i] = sqrtf(RE->data[i]*RE->data[i] + IM->data[i]*IM->data[i]);
	}

	if (GMT_Destroy_Data (API, &IM)) die("error freeing data",imname);
	if (GMT_Destroy_Data (API, &RE)) die("error freeing data",rename);

	return (EXIT_SUCCESS);
}

int make_wgt(float *wgt, int nxp, int nyp)
{
int	i, j;
float	wx,wy;

	for (i=0; i<nyp/2; i++){
		wy = 1.0f - fabsf((float) (i) - (nyp/2.0f - 1.0f)) / (nyp/2.0f - 1.0f);
		for (j=0; j<nxp/2; j++){
			wx = 1.0f - fabsf((float) (j) - (nxp/2.0f - 1)) / (nxp/2.0f - 1.0f);
			wgt[i*nxp+j] = wgt[(i+1)*nxp-j-1] = wx*wy;
			wgt[(nyp-i-1)*nxp+j] = wgt[(nyp-i)*nxp-j-1] = wx*wy;
			}
		}

	if (debug) print_patch(nxp, nyp, wgt);

	return(EXIT_SUCCESS);
}

/* the classic Goldstein filter							*/
int apply_pspec(int m, int n, float alpha, struct FCOMPLEX *in, struct FCOMPLEX *out)
{
int 	i;
float   wgt;

	if (alpha < 0.0f) die("alpha < 0; something is rotten in Denmark","");
	/* pow(x,a/2) == pow(sqrt(x),a) */  
	for (i=0; i<m*n; i++) {
		wgt = powf((in[i].r*in[i].r + in[i].i*in[i].i),alpha/2.0f);
		out[i].r = wgt*in[i].r;	
		out[i].i = wgt*in[i].i;	
		}

	return(EXIT_SUCCESS);
}

int calc_corr(void *API, char *amp1,  char *amp2, unsigned int xdim, unsigned int ydim, float *amp, float *corr)
{
	unsigned int i, n, xdim2, ydim2;
	float a;
	struct GMT_GRID *A1 = NULL, *A2 = NULL;	/* Grid structures containing ->header and ->data */

	n = xdim * ydim;

	read_file_hdr (API, amp1, &A1, amp2, &A2, &xdim2, &ydim2);

	if ((xdim != xdim2) || (ydim != ydim2)) die("amp files are different size than real and imag files","");

	if (GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, amp1, A1) == NULL) die("error reading data",amp1);
	if (GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, amp2, A2) == NULL) die("error reading data",amp2);

	for (i=0; i<n; i++) {
		a = A1->data[i] * A2->data[i];
		if (a > 0.0f) {
			corr[i] = amp[i] / sqrtf(a);
		} else { 
			corr[i] = 0.0f;
		}
		if (corr[i] < 0.0f) corr[i] = 0.0f;
		if (corr[i] > 1.0f) corr[i] = 1.0f;
	}
	if (GMT_Destroy_Data (API, &A1)) die("error freeing data",amp1);
	if (GMT_Destroy_Data (API, &A2)) die("error freeing data",amp2);
	return(EXIT_SUCCESS);
}

int phasefilt_parse_command_line(char **a, int na, char *USAGE, char *sre, char *sim, float *alp, int *ps, char *amp1, char *amp2, int *dflag, int *comflag)
{
int	n;
int	flag[4];

for (n=0; n<4; n++) flag[n] = 0;

for(n = 1; n < na; n++) {
	if (!strcmp(a[n], "-v")) {
		verbose = 1;
		/* only args after verbose (if set) will be echoed */
	} else if (!strcmp(a[n], "-real")) {
		n++;
		if (n == na) die (" no option after -real!\n","");
		strcpy(sre, a[n]);
		flag[0] = 1;
		if (verbose) fprintf(stderr, "real %s \n", sre);
	} else if (!strcmp(a[n], "-imag")) {
		n++;
		if (n == na) die (" no option after -imag!\n","");
		strcpy(sim, a[n]);
		flag[1] = 1;
		if (verbose) fprintf(stderr, "imag %s \n", sim);
	} else if (!strcmp(a[n], "-alpha")) {
		n++;
		if (n == na) die (" no option after -alpha!\n","");
		*alp = (float)strtod(a[n],NULL);
		if (verbose) fprintf(stderr, "alpha %f \n", *alp);
	} else if (!strcmp(a[n], "-amp1")) {
		n++;
		if (n == na) die (" no option after -amp1!\n","");
		strcpy(amp1, a[n]);
		flag[2] = 1;
		if (verbose) fprintf(stderr, "amp1 %s \n", amp1);
	} else if (!strcmp(a[n], "-amp2")) {
		n++;
		if (n == na) die (" no option after -amp2!\n","");
		strcpy(amp2, a[n]);
		flag[3] = 1;
		if (verbose) fprintf(stderr, "amp2 %s \n", amp2);
	} else if (!strcmp(a[n], "-psize")) {
		n++;
		if (n == na) die (" no option after -psize!\n","");
		*ps = atoi(a[n]);
		if (verbose) fprintf(stderr, "patch size %d \n", *ps);
	} else if (!strcmp(a[n], "-diff")) {
		*dflag = 1;
		if (verbose) fprintf(stderr, "calculate difference \n");
	} else if (!strcmp(a[n], "-complex_out")) {
		*comflag = 1;
	} else if (!strcmp(a[n], "-debug")) {
		verbose = 1;
		debug = 1;
	} else {
		fprintf(stderr," %s *** option not recognized ***\n\n",a[n]);
		fprintf(stderr,"%s",USAGE);
		exit(1);
		}
	}

	/*
	flag[0]	real 
	flag[1]	imag 
	flag[2]	amp1 	only needed if modified filter
	flag[3]	amp2	only needed if modified filter 
	*/

	/* check for real and imag files are set */
	if (flag[0] == 0) die("real file not set","");
	if (flag[1] == 0) die("imag file not set","");

	/* if amp files are defined set alpha = -1 to denote coherence-based alpha */
	if ((flag[2] == 1) && (flag[3] == 1)) *alp = -1.0f; 

	if ((flag[2] == 1) && (flag[3] == 0)) die("amp1 set but not amp2 (needed for modified filter)","");
	if ((flag[2] == 0) && (flag[3] == 1)) die("amp2 set but not amp1 (needed for modified filter)","");

	/* if alpha < 0 check that both amp files are set */
	if ((*alp < 0) && (flag[2] == 0)) die("amp1 file not set (needed for modified filter)","");
	if ((*alp < 0) && (flag[3] == 0)) die("amp2 file not set (needed for modified filter)","");

return(EXIT_SUCCESS);
}

char *USAGE  = "phasefilt [GMT5SAR] - Apply adaptive non-linear phase filter\n\n"
"\n USAGE:\nphasefilt -imag imag.grd -real real.grd [-alpha alpha][-psize size][-amp1 amp1.grd -amp2 amp2.grd][-diff][-v]\n"
" applies Goldstein adaptive filter to phase [output: filtphase.grd]\n"
" or applies modified Goldstein adaptive filter to phase [output: filtphase.grd, corrfilt.grd]\n"
"-imag [required] GMT format file of imaginary component\n"
"-real [required] GMT format file of real component\n"
"-alpha 	exponent for filter - usually between 0.0 and 1.5 (0.0 should not filter).\n"
"		default: 0.5	[Goldstein filter] (anything above 1.0 may be excessive)\n"
"		alpha < 0 will set alpha = (1 - coherence) [modified Goldstein]\n"
"-psize 	patch size for filtering. Must be power of two.\n"
"		default: 32\n"
"-amp1 		GMT format file of amplitude image of image 1. Needed (and applies) modified filter.\n"
"-amp2 		GMT format file of amplitude image of image 2. Needed (and applies) modified filter.\n"
"-diff 		Calculate difference between input phase and output phase.\n"
"-complex_out	Write out filtered real and imaginary (filtphase_real.grd and filtphase_imag.grd)\n"
"-v 		Verbose.\n"
"-debug 	Debug.\n"
"\n"
"example:\n"
"phasefilt -imag imag.grd -real real.grd -alpha 0.5\n";

int main(int argc, char **argv){
	unsigned int	i, j, ii, jj, k1, k2;
	unsigned int	xdim, ydim;
	unsigned int	nxp, nyp; 
	int	psize, corrflag, dflag, comflag;
	float	alpha, pcorr, swgt;
	float	*outphase = NULL, *wgt = NULL, *amp = NULL, *corr = NULL, *diff = NULL, *ftmp = NULL;
	struct 	FCOMPLEX *data = NULL, *fdata = NULL, *patch0 = NULL, *patch1 = NULL;
	char	sre[256], sim[256];
	char	amp1[256], amp2[256];
	
	void *API = NULL; /* GMT API control structure */
	struct	GMT_GRID *RE = NULL, *IM = NULL, *T = NULL;	/* Grid structure containing ->header and ->data */

	/* Begin: Initializing new GMT5 session */
	if ((API = GMT_Create_Session (argv[0], 0U, 0U, NULL)) == NULL) return EXIT_FAILURE;

	debug = verbose = corrflag = comflag = 0;
	if (argc < 3) die(USAGE,"");

	/* defaults */
	psize = 32;		/* size of patch  # changed from 64 */
	alpha = 0.5f;		/* exponent */
	dflag = 0;		/* write out difference */
	phasefilt_parse_command_line(argv, argc, USAGE, sre, sim, &alpha, &psize, amp1, amp2, &dflag, &comflag);

	/* patch size 		*/
	/* currently square 	*/
	nxp = psize;
	nyp = psize;

	/* read dimensions 	*/
	read_file_hdr (API, sim, &IM, sre, &RE, &xdim, &ydim);
	if ((T = GMT_Create_Data (API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, IM->header->wesn, IM->header->inc, \
		IM->header->registration, 0, NULL)) == NULL) die("could not allocate grid header","");

	if ((data   = malloc (T->header->nm * sizeof (struct FCOMPLEX))) == NULL)  die("error allocating memory","");
	if ((fdata  = malloc (T->header->nm * sizeof (struct FCOMPLEX))) == NULL) die("error allocating memory","");
	if ((patch0 = malloc (nxp * nyp * sizeof (struct FCOMPLEX))) == NULL)     die("error allocating memory","");
	if ((patch1 = malloc (nxp * nyp * sizeof (struct FCOMPLEX))) == NULL)     die("error allocating memory","");
	if ((amp = malloc (T->header->nm * sizeof (float))) == NULL) die("error allocating memory","");
	if ((wgt = malloc (nxp * nyp * sizeof (float))) == NULL) die("error allocating memory","");
	if ((outphase = malloc (xdim * (ydim+nyp) * sizeof (float))) == NULL) die("error allocating memory","");
	
	/* array for difference */
	if (dflag) {if ((diff = malloc (T->header->nm * sizeof (float))) == NULL) die("error allocating memory","");}

	/* array for writing filtered real and imag  */
	if (comflag) {if ((ftmp = malloc (T->header->nm * sizeof (float))) == NULL) die("error allocating memory","");}

	for (i=0; i<T->header->nm; i++) outphase[i] = fdata[i].r = fdata[i].i = 0.0f;

	/* read data 	*/
	phasefilt_read_data (API, sim, IM, sre, RE, data, amp);

	/* calculate correlation 	*/
	if (alpha < 0.0) { 
		corr = (float *) malloc(T->header->nm * sizeof(float));
		calc_corr (API, amp1, amp2, xdim, ydim, amp, corr); 
		corrflag = 1;
	}

	if (verbose) {
		if (corrflag == 1) fprintf(stderr,"phasefile: using coherence-dependent alpha (1 - gamma)\n");
		if (corrflag == 0) fprintf(stderr,"phasefilt: constant alpha (%6.2f)\n", alpha);
	}

	/* create weights for each patch 				*/
	/* each patch overlaps each other by half and summed 		*/
	/* ideally, total wgt for each pixl = 1				*/
	/* except at edges...						*/
	make_wgt(wgt, nxp, nyp);

	for (ii=0; ii<(ydim-nyp); ii+=nyp/2){
		for (jj=0; jj<(xdim-nxp); jj+=nxp/2){

			pcorr = 0.0f;
			swgt = 0.0f;
			for (i=0; i<nyp; i++){
				for (j=0; j<nxp; j++) {
					k1 = (ii+i)*xdim + jj + j;
					k2 = i*nxp+j;
					patch0[k2].r = data[k1].r;
					patch0[k2].i = data[k1].i;
					if (corrflag) {
						pcorr += wgt[k2]*corr[k1];
						swgt  += wgt[k2];
						}
					}
				}

			/* set alpha to 1.0 - coherence 	*/
			/* Baran et al., 2003			*/
			if (corrflag) alpha = 1.0f - pcorr/swgt;

			GMT_FFT_2D (API, (float *)patch0, nxp, nyp, GMT_FFT_FWD, GMT_FFT_COMPLEX);

			apply_pspec(nxp, nyp, alpha, patch0, patch1);

			GMT_FFT_2D (API, (float *)patch1, nxp, nyp, GMT_FFT_INV, GMT_FFT_COMPLEX);

			for (i=0; i<nyp; i++){
				for (j=0; j<nxp; j++) {
					k1 = (ii+i)*xdim +jj+j;
					k2 = i*nxp+j;
					fdata[k1].r = fdata[k1].r + wgt[k2]*patch1[k2].r;
					fdata[k1].i = fdata[k1].i + wgt[k2]*patch1[k2].i;
					}
				}
			}
		}

	for (i=0; i<T->header->nm; i++) outphase[i] = atan2f(fdata[i].i,fdata[i].r);
	
	write_grdfile (API, T, "filtphase.grd", "phasefilt", "phase", outphase, verbose);

	if (corrflag) write_grdfile (API, T, "filtcorr.grd", "phasefilt", "corr", corr, verbose);

	if (comflag) {
		for (i=0; i<T->header->nm; i++)  ftmp[i] = fdata[i].r; 
		write_grdfile (API, T, "filtphase_real.grd", "phasefilt", "real", ftmp, verbose);
		for (i=0; i<T->header->nm; i++)  ftmp[i] = fdata[i].i; 
		write_grdfile (API, T, "filtphase_imag.grd", "phasefilt", "imag", ftmp, verbose);
		free ((void *)ftmp);
	}

	if (dflag) {
		for (i=0; i<T->header->nm; i++) diff[i] = atan2f(data[i].i,data[i].r) - outphase[i];
		write_grdfile (API, T, "filtdiff.grd", "phasefilt", "diff", diff, verbose);
		free ((void *)diff);
	}

	if (GMT_Destroy_Session (API)) return EXIT_FAILURE;	/* Remove the GMT machinery */

	free ((void *)patch0);
	free ((void *)patch1);
	free ((void *)data);
	free ((void *)fdata);
	free ((void *)amp);
	free ((void *)outphase);
	free ((void *)wgt);
	
	return(EXIT_SUCCESS);
}
