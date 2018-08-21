/*	$Id: conv.c 109 2015-01-19 23:01:24Z sandwell $	*/
/***************************************************************************/
/* conv convolves a 2-D filter with an array and outputs the results       */
/* at a sub-sampled interval.  Basically it does the same thing as         */
/* the GIPS program ihbox but it runs a lot faster.                        */
/*                                                                         */
/***************************************************************************/

/***************************************************************************
 * Creator:  David T. Sandwell                                             *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  04/17/98                                                      *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * DATE					                                   *
 * 04/17/98     Started hacking at the code                                *
 * 06/07/98     Removed edge effects of derivative filter                  *
 * 01/14/10     Modified to read and write grd files - RJM                 *
 ***************************************************************************/

#include "gmtsar.h"

char *USAGE = "conv [GMTSAR] - 2-D image convolution\n\n"
"Usage: conv idec jdec filter_file input output \n"
"   idec           - row decimation factor \n"
"   jdec           - column decimation factor \n"
"   filter_file    - eg. filters/gauss17x5 \n"
"   input          - name of file to be filtered (I*2 or R*4) \n"
"   output         - name of filtered output file (R*4 only) \n\n"
"   examples:\n"
"   conv 4 2 filters/gauss9x5 IMG-HH-ALPSRP109430660-H1.0__A.PRM test.grd \n"
"   (makes and filters amplitude file from an SLC-file) \n\n"
"   conv 4 2 filters/gauss5x5 infile.grd outfile.grd \n"
"   (filters a float file) \n ";

int	input_file_type, format_flag;

/*-------------------------------------------------------------*/
int determine_file_type(char *name, int *input_file_type)
{
int	n, m;
char	tail[8];

	*input_file_type = 1;

	n = (int)strlen(name);
	m = n - 3;
	strncpy(&tail[0], &name[m], 4);
	if (verbose) fprintf(stderr," name %s tail %s \n",name, tail);

	if ((strncmp(tail, "PRM", 3) == 0) || (strncmp(tail, "prm", 3) == 0)){
		if (verbose) fprintf(stderr," input: PRM file\n");
		*input_file_type = 2;
		}

	if (*input_file_type == 1) if (verbose) fprintf(stderr," input: GMT binary\n");

	return(EXIT_SUCCESS);
}

/*-------------------------------------------------------------*/
FILE	*read_PRM_file(char *prmfilename, char *input_file_name, struct PRM p, int *xdim, int *ydim)
{
FILE *f_input_prm, *f_input;
void change_name(char * );

	if (verbose) fprintf(stderr," reading PRM file %s\n",prmfilename);
	if ((f_input_prm = fopen(prmfilename,"r")) == NULL) die("Can't open input header",prmfilename);
	get_sio_struct(f_input_prm, &p);
	strcpy(input_file_name, p.SLC_file); 
	format_flag = 2;
	if(strncmp(p.dtype,"c",1) == 0) format_flag = 3;
	if (verbose) fprintf(stderr," reading PRM file %s\n",input_file_name);
	if ((f_input = fopen(input_file_name,"r")) == NULL) die("Can't open input data ",input_file_name);
	*xdim = p.num_rng_bins;
	*ydim = p.num_valid_az * p.num_patches;

	return(f_input);
}

/*-------------------------------------------------------------*/
int read_float(float *indat, int xdim, FILE *f_input, int yarr, float *buffer, int ibuff)
{
int i, j;

	for (i=0;i<ibuff;i++){
		fread(indat, sizeof(float), xdim, f_input);
		for (j=0; j<xdim; j++) buffer[j+xdim*(i+yarr)] = indat[j];
		}

	return(EXIT_SUCCESS);
}
/*-------------------------------------------------------------*/
int read_SLC_int(short *ci2, int xdim, FILE *f_input, int yarr, float *buffer, double dfact, int ibuff)
{
int i, j;
double	df2;

	df2 = dfact*dfact;

	/* read i2 complex and calculate amplitude */
	/* use square of amplitude to match gips ihconv */
	for (i=0;i<ibuff;i++){
		fread(ci2, 2*sizeof(short), xdim, f_input);
		for (j=0; j<xdim; j++) buffer[j+xdim*(i+yarr)]
			= (float) (df2*ci2[2*j]*ci2[2*j] + df2*ci2[2*j+1]*ci2[2*j+1]);
	}

	return(EXIT_SUCCESS);
}

/*-------------------------------------------------------------*/
int read_SLC_float(float *cf2, int xdim, FILE *f_input, int yarr, float *buffer, double dfact, int ibuff)
{
	int i, j;
	double	df2;
	
	df2 = dfact*dfact;
	
	/* read r4 complex and calculate amplitude */
	/* use square of amplitude to match gips ihconv */
	for (i=0;i<ibuff;i++){
		fread(cf2, 2*sizeof(float), xdim, f_input);
		for (j=0; j<xdim; j++) buffer[j+xdim*(i+yarr)]
			= (float) (df2*cf2[2*j]*cf2[2*j] + df2*cf2[2*j+1]*cf2[2*j+1]);
	}
	return(EXIT_SUCCESS);
}

int main(int argc, char **argv)
{
	int	idec, jdec;
	int	iout, jout;
	int	i, j, ic, jc, norm, ic0, ic1; 
	int	ydim=0, xdim=0;		/* size of input file */
	int	xarr, yarr, narr, yarr2;
	int	nbuff, ibuff, imove;
	int	iend, ylen, iread;
	uint64_t	left_node;
	unsigned int row;
	char	input_name[128], output_name[128], prmfilename[128], *c = NULL;
	short	*cindat = NULL;
	float	*cfdat = NULL;
	double	inc[2], wesn[4], xmax = 0.0, ymax = 0.0;
	float	*filter = NULL,*buffer = NULL,*indat = NULL;
	float 	filtin, filtdat,rnorm,rnormax,anormax;
	FILE	*f_filter = NULL, *f_input = NULL;
	struct	PRM p;
	void	*API = NULL; /* GMT control structure */
	struct	GMT_GRID *Out = NULL;	/* For the output grid */
	struct	GMT_GRID *In = NULL;	/* Grid structure containing ->header and ->data */

	if (argc < 6) die("\n",USAGE);

	/* Begin: Initializing new GMT session */
	if ((API = GMT_Create_Session (argv[0], 0U, 0U, NULL)) == NULL) return EXIT_FAILURE;

	ibuff = 512;
	verbose = 0;

	null_sio_struct(&p);
	input_file_type = 1;	/* default; GMT binary */

	/* format_flag = 1 => float		*/
	/* format_flag = 2 => i*2 complex	*/
	/* format_flag = 3 => r*4 complex	*/
	format_flag = 1;

	idec = atoi(argv[1]);		/* y decimation factor */
	jdec = atoi(argv[2]);		/* x decimation factor */
	
	if (idec <= 0 || jdec <= 0) die("idec and jdec should be positive integers.", "");
	
	/* open and filter file */
	if ((f_filter = fopen(argv[3],"r")) == NULL) die("Can't open filter","");

	/* get input and output file names */
	strcpy(input_name,argv[4]);
	strcpy(output_name,argv[5]);

	/* read input file 	*/
	/* GMT file or PRM file ? look at suffix (grd or PRM)*/
	determine_file_type(input_name, &input_file_type);

	if (verbose) fprintf(stderr," input file type: %d \n", input_file_type);

	switch (input_file_type)
	{
		case 1:
			if (verbose) fprintf(stderr," reading GMT binary\n");
			if ((In = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, input_name, NULL)) == NULL) die("Can't open ",input_name);
			if ((c = strstr (input_name, "=bf"))) c[0] = '\0';	/* Chop off any trailing =bf flag */
			if ((f_input = fopen(input_name,"r")) == NULL) die("Can't open ",input_name);
			fseek (f_input, 892L, SEEK_SET);	/* Skip past the header */
			xdim=In->header->nx;
			ydim=In->header->ny;
			xmax=In->header->wesn[GMT_XHI];
			ymax=In->header->wesn[GMT_YHI];
			if (verbose) fprintf(stderr,"%d %d \n",In->header->nx,In->header->ny);
			format_flag = 1;
			break;
		case 2:
			strcpy(prmfilename, input_name);
			f_input = read_PRM_file(prmfilename, input_name, p, &xdim, &ydim);
			if (verbose) fprintf(stderr," reading PRM file format %d\n", format_flag);
			xmax=xdim;
			ymax=ydim;
			break;
		default:
			die("confused about input file type","quitting");
	}

	if (verbose) fprintf(stderr," read file %s %d %d \n",input_name, xdim, ydim);

	/* read size of filter and make sure dimensions are odd */
	if (fscanf(f_filter,"%d%d", &xarr, &yarr) != 2 || xarr < 1 || yarr < 1 || (xarr & 1) == 0 || (yarr & 1) == 0) 
			die("filter incomplete","");
	if (ibuff < yarr) die("increase dimension of ibuff","");

	/* size of output file */
	iout = jout = 0;
	for (ic = 0; ic<ydim; ic = ic+idec) iout = iout + 1;
	for (jc = 0; jc<xdim; jc = jc+jdec) jout = jout + 1;
	if (verbose) fprintf(stderr," original: ydim %d xdim %d new: %d %d decimation: %d %d\n",ydim,xdim,iout,jout,idec,jdec); 
	inc[GMT_X] = round(xmax/(double)jout);
	inc[GMT_Y] = round(ymax/(double)iout);
    jout = floor(xmax/inc[GMT_X]);
    iout = floor(ymax/inc[GMT_Y]);
	if (verbose) fprintf(stderr," creating GMT grid %s \n",output_name);
	wesn[GMT_XLO] = 0.0;	wesn[GMT_XHI] = inc[GMT_X] * jout;
	wesn[GMT_YLO] = 0.0;	wesn[GMT_YHI] = inc[GMT_Y] * iout;

	if ((Out = GMT_Create_Data (API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, wesn, inc, \
		GMT_GRID_PIXEL_REG, 0, NULL)) == NULL) die("could not allocate output grid","");
	if (GMT_Set_Comment (API, GMT_IS_GRID, GMT_COMMENT_IS_TITLE, "conv", Out)) die("could not set title","");

	/* parameters for convolution buffer */
	narr  = xarr * yarr;
	yarr2 = (int) (yarr/2.0);
	if (ydim < ibuff) ibuff = ydim;
	imove = ibuff - yarr;
	nbuff = xdim * ibuff;

	if (( filter = (float *) malloc(sizeof(float) * narr)) == NULL) die("memory allocation","");
	if (( buffer = (float *) malloc(2 * sizeof(float) * nbuff)) == NULL) die("memory allocation","");

	if (format_flag == 1) if (( indat = (float *) malloc(4*xdim)) == NULL) die("memory allocation",""); 
	if (format_flag == 2) if (( cindat = (short *) malloc(4*xdim)) == NULL) die("memory allocation",""); 
	if (format_flag == 3) if (( cfdat = (float *) malloc(8*xdim)) == NULL) die("memory allocation","");

	/* read the filter and calculate normalization constants*/
	anormax = rnormax = 0.0f;
	for(i=0;i<narr;i++) {
		if (fscanf(f_filter,"%f",&filtin) == EOF) die("filter incomplete","");
		filter[i] = filtin;
		anormax = anormax + (float) fabs(filter[i]);
		rnormax = rnormax + filter[i];
		}

	norm = 0.0f;
	if(fabs(rnormax) > 0.05*anormax) norm = 1.0f;
	ic0 = 0;
	iend = ylen = ibuff;

	if (verbose) fprintf(stderr," format_flag %d \n", format_flag);
	
	/* read the data (512 lines) */
	if (format_flag == 1) read_float(indat, xdim, f_input, 0, buffer, ibuff); 
	if (format_flag == 2) read_SLC_int(cindat, xdim, f_input, 0, buffer, DFACT, ibuff); 
	if (format_flag == 3) read_SLC_float(cfdat, xdim, f_input, 0, buffer, DFACT, ibuff);

	for (ic=0, row = 0; ic<iout*idec; ic=ic+idec){

	    if (ic/2000.0 == ic/1000) if (verbose) fprintf(stderr," line %d\n",ic);

	    /* check buffer and shift data up if necessary */
	    if ((ic+yarr2) >= iend && (ic+yarr2) < (ydim-1)) {

	        for (i=0;i<yarr;i++) for (j=0;j<xdim;j++) buffer[j+xdim*i] = buffer[j+xdim*(i+(ylen-yarr))];

                iread = MIN(imove,(ydim - iend));
		ic0 = iend - yarr;
                iend = iend + iread;
                ylen = iread + yarr;

		/* now read in more data into end of buffer */
			if (format_flag == 1) read_float(indat, xdim, f_input, yarr, buffer, iread);
			if (format_flag == 2) read_SLC_int(cindat, xdim, f_input, yarr, buffer, DFACT, iread);
			if (format_flag == 3) read_SLC_float(cfdat, xdim, f_input, yarr, buffer, DFACT, iread);

		} /* end of ic loop */  
		left_node = GMT_Get_Index (API, Out->header, row, 0);
	        jout = 0;
	        ic1 = ic - ic0;

			/* now do the 2d convolution */
	        for (jc=0;jc<floor(xmax/inc[GMT_X])*jdec;jc=jc+jdec) {
	            	conv2d (buffer, &ylen, &xdim, filter, &yarr, &xarr, &filtdat, &ic1, &jc, &rnorm);
					/* use a zero or null value if there is not enough data in the filter */
	            	Out->data[left_node+jout] = 0.0f;
	            	if (norm > 0 ) {
	    	       	 	if (fabs(rnorm) > (0.01*rnormax)) Out->data[left_node+jout] = filtdat/rnorm;
		       	 	} else { 
				if (fabs(rnorm) < 0.0001*anormax) Out->data[left_node+jout] = filtdat; 
			}

			jout++;

	        } /* end of jc loop */
		row++;
	} /* end of data loop */
	fclose (f_input);

	if (GMT_Write_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, output_name, Out)) {
		die("Failed to write output grid","");
	}
	
	if (GMT_Destroy_Session (API)) return EXIT_FAILURE;	/* Remove the GMT machinery */

	return(EXIT_SUCCESS);
}
