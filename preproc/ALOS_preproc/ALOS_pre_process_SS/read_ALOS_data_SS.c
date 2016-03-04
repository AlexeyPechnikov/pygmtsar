/***************************************************************************/
/* read_ALOS_data_SS reads an ALOS IMG file containing raw signal data     */
/* and creates a raw-file and PRM-file suitable for our esarp processor.   */
/* The program skips the first 720 bytes of the IMG file but copies the    */
/* remaining data to the IMG.raw file after checking and fixing problems.  */
/* The first record is read to determine the linelength, starting PRF,     */
/* and near_range.  If the line length or PRF change then the program      */
/* halts.  If the near_range changes then the lines are shifted and        */
/* unconstrained values at the ends are set to NULL_DATA (15 or 16).       */
/* (random sequence of 15 and 16's)                                        */
/* During this processing the available parameters are added to the        */
/* PRM-file.                                                               */
/***************************************************************************/

/***************************************************************************
 * Creator:  David T. Sandwell and Rob Mellors                             *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  10/03/2008                                                    *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * DATE                                                                    *
 *                                                                         *
 * 06/29/2006   added the near_range as an optional command-line argument  *
 * 02/19/2007   added the ability to remove duplicate lines                *
 * 03/07/2007   more robust to bad data in file                            *
 * 03/26/2007   ability to swap bytes (read on PC) RJM                     *
 * 03/28/2007	part of subroutin  RJM                                     *
 * removed set n_azimuth to 9000 rather than default                       *
 * 07/17/08     creates new file when prf changes RJM                      *
 * 07/17/08     reformatted; added functions      RJM                      *
 * 09/29/08     added the ability to dump a subswath of WB1 scansar        *
 ***************************************************************************/

/* the data header information is read into the structure dfd
the line prefix information is read into sdr
Values read here (and or generated) are:
num_rng_bins bytes_per_line good_bytes_per_line 
PRF pulse_dur near_range
num_lines num_patches 
SC_clock_start SC_clock_stop  - this includes year YYYDDD.DDDDDDD
clock_start clock_stop        - this has not year for better precision DDD.DDDDDDDDD
*/

#include "image_sio.h"
#include "lib_functions.h"

/* fast random number generator */
#define znew   (int) (z=36969*(z&65535)+(z>>16))
typedef unsigned long UL;
static UL z=362436069, t[256];
void settable(UL i1)
{ int i; z=i1;
for(i=0;i<256;i=i+1)  t[i]=znew;
}

/* burst information from JAXA manual for WB1 */
static int n_data_burst[5]={4976,4720,5376,4432,4688};
static int n_burst[5]={247,356,274,355,327};
static int PRF[5]={0,0,0,0,0};

void swap_ALOS_data_info(struct sardata_info *sdr);
long read_sardata_info(FILE *, struct PRM *, int *, int *);
void print_params(struct PRM *prm);
int assign_sardata_params(struct PRM *, int, int *, int *);
int check_shift(struct PRM *, int *, int *, int *, int);
int set_file_position(FILE *, long *, int);
int fill_shift_data(int, int, int, int, int, char *, char *, FILE *);
int handle_prf_change(struct PRM *, FILE *, long *, int);

struct 	sardata_record r1;
struct	sardata_descriptor dfd;
struct	sardata_info sdr;

long read_ALOS_data_SS (FILE *imagefile, FILE *outfile, struct PRM *prm, long *byte_offset, 
                        int *nsub, int *burst_skip, int *num_burst) {

	char *data, *shift_data, *gap_data;
	int	header_size;		/* file header size          720 bytes */
	int	line_prefix_size; 	/* line header size           412 bytes*/
	int     record_length0;		/* data record size start  10788 bytes */
	int	totrecl;		/* total record length	   11200 bytes = line_prefix_size + record_length */
	int	record_length1;		/* data record size curr.  10788 bytes */
	int	line_suffix_size;	/* number of bytes after data 36 bytes */
	int	kswath0=0;		/* subswath at start of file */
	int	skip_swath;             /* number of subswaths to skip to get to desired subswath */
	int	skiprecl;               /* number of lines to skip to get to desired subswath */
	int	data_length;		/* bytes of data			*/
	int     n = 1, ishift, shift, shift0;
	int	k, kburst = 0, totburst = 0;
	int	j, ngap, nlines = 0, ntot;
	int	nprfchange, nburstchange;

	double 	tbias = 0.0, get_clock();
	double  ttot=0.,dt=0.,tgap=0.;		/* total time, burst interval, burst gap, fractional gap */
	settable(12345);

	if (debug) fprintf(stderr,".... reading header \n");

	/* read header information */
	read_sardata_info(imagefile, prm, &header_size, &line_prefix_size);
	assign_sardata_params(prm, line_prefix_size, &line_suffix_size, &record_length0);
	totrecl = line_prefix_size + record_length0;
	set_file_position(imagefile, byte_offset, header_size);

	/* get the sarting burst number and the PRF information for all the bursts and rewind the file */
	for (k=0;k<5;k++) totburst = totburst + n_burst[k];
	for(j = 0; j < totburst; j++) {
		if ( fread((void *) &sdr,sizeof(struct sardata_info), 1, imagefile) !=1 ) break;
		fseek(imagefile, record_length0, SEEK_CUR);
		if (swap) swap_ALOS_data_info(&sdr);
		for (k=0;k<5;k++) {
			if(sdr.n_data_pixels == n_data_burst[k]) {
				PRF[k]=sdr.PRF;
				if(j == 0) kswath0 = k;
			}
		}
	}

        /* get the time interval of a burst  and rewind the file */
	dt = 0.;
	for (k=0;k<5;k++) if(PRF[k] > 0.)  dt = dt + (double)n_burst[k]/(0.001*PRF[k]);
 	tgap = dt - (double)n_burst[*nsub]/(0.001*PRF[*nsub]);
        ngap = (int)(tgap*(0.001*PRF[*nsub])+0.5);
	prm->num_valid_az = 6*(n_burst[*nsub] + ngap);    /* use 6 bursts per aperture */
	rewind(imagefile);
	fseek(imagefile, header_size, SEEK_SET);
	if (verbose) fprintf(stderr," totburst start_busrt# burstcycle_time %d %d %f \n", totburst, kswath0, dt);

	/* set the total number of lines output based on the number of patches requested or default 1000 */
	ntot = *num_burst * (n_burst[*nsub] + ngap);
	if (verbose) fprintf(stderr," dt, tgap, n_burst, ngap ,ntot %f %f %d %d %d %d \n", dt, tgap, n_burst[*nsub], ngap, prm->num_valid_az,ntot);

	/* allocate the memory for data */
	if ((data = (char *) malloc(record_length0)) == NULL) die("couldn't allocate memory for input indata.\n","");
	if ((shift_data = (char *) malloc(record_length0)) == NULL) die("couldn't allocate memory for input indata.\n","");
	if ((gap_data = (char *) malloc(record_length0)) == NULL) die("couldn't allocate memory for input indata.\n","");

	/* seek to beginning of next burst of sub swath nsub, read the line header, reset the parameters, and set the counters */
	skip_swath = ((*nsub +5) - kswath0)%5;
	skiprecl = 0;
	for( k = kswath0; k < kswath0 + skip_swath; k++) skiprecl = skiprecl + n_burst[k%5];
	skiprecl = skiprecl + *burst_skip * totburst;
	if (verbose) fprintf(stderr, " skip_swath skiprecl %d %d \n",skip_swath, skiprecl);
	fseek(imagefile,skiprecl*totrecl, SEEK_CUR);

        /* recalculate the PRF at the new skip location */
	for(j = 0; j < totburst; j++) {
		if ( fread((void *) &sdr,sizeof(struct sardata_info), 1, imagefile) !=1 ) break;
		fseek(imagefile, record_length0, SEEK_CUR);
		if (swap) swap_ALOS_data_info(&sdr);
		for (k=0;k<5;k++) {
			if(sdr.n_data_pixels == n_data_burst[k]) {
				PRF[k]=sdr.PRF;
			}
		}
	}
        /* get the time interval of a burst at the new skip location and rewind the file */
	dt = 0.;
	for (k=0;k<5;k++) if(PRF[k] > 0.)  dt = dt + (double)n_burst[k]/(0.001*PRF[k]);
 	tgap = dt - (double)n_burst[*nsub]/(0.001*PRF[*nsub]);
        ngap = (int)(tgap*(0.001*PRF[*nsub])+0.5);
	prm->num_valid_az = 6*(n_burst[*nsub] + ngap);    /* use 6 bursts per aperture */
	if (verbose) fprintf(stderr," totburst start_busrt# burstcycle_time %d %d %f \n", totburst, kswath0, dt);

	/* rewind and seek again to start in the correct location */
	rewind(imagefile);
	fseek(imagefile, header_size, SEEK_SET);
	fseek(imagefile,skiprecl*totrecl, SEEK_CUR);
	
        /* read the line header and set the parameters for this subswath */
	fread((void *) &sdr,line_prefix_size, 1, imagefile);
	if (swap) swap_ALOS_data_info(&sdr);
        prm->clock_start =  get_clock(sdr, tbias);
        prm->SC_clock_start = ((double) sdr.sensor_acquisition_year)*1000 + prm->clock_start;
	prm->prf = sdr.PRF;
	if(prm->near_range < 0) prm->near_range = sdr.slant_range;
	fseek(imagefile, -1*line_prefix_size, SEEK_CUR);
	n = sdr.sequence_number -1;
	//m = sdr.sequence_number;
	*byte_offset = ftell(imagefile);
	if (verbose) fprintf(stderr," n skiprecl byte_offset %d %d %ld \n", n, skiprecl, *byte_offset);

	/* now start at the beginning but with the intervals known */
	if (verbose) fprintf(stderr,".... reading data (byte %ld) \n",ftell(imagefile));
	shift0 = 0;
  	//m = 0;

	/* read the rest of the file */
	while ( (fread((void *) &sdr,sizeof(struct sardata_info), 1, imagefile)) == 1 ) {
        	n++;
		if (swap) swap_ALOS_data_info(&sdr);

	/* check to make sure there is no prf-change in any of the subswaths */
		for (k=0;k<5;k++) {
			if(sdr.n_data_pixels == n_data_burst[k]) {
          		if ((sdr.PRF) != PRF[k]) {
				PRF[*nsub] = 0.;
			}
			}
		}

		/* detect a different subswath and skip intil the next subswath */
		if(sdr.n_data_pixels == n_data_burst[*nsub]) {

			kburst++;

          		if (sdr.sequence_number != n) printf(" missing line: n, seq# %d %d \n", n, sdr.sequence_number);

			/* check for changes in record_length and PRF */
          		record_length1 = sdr.record_length - line_prefix_size;
          		if (record_length0  != record_length1)  die("record_length changed",""); 

			/* if prf changes exit */
          		if ((sdr.PRF) != PRF[*nsub]) {
				fprintf(stderr," ERROR  PRF changed, oldPRF, newPRF %f %f \n",PRF[*nsub]*.001,sdr.PRF*.001);
				*byte_offset=ftell(imagefile);
				nprfchange = (*byte_offset - header_size)/totrecl;
				nburstchange = nprfchange/totburst;
				fprintf(stderr," rec# burst# %d %d \n",nprfchange,nburstchange);
                                break;
			}

			/* check shift to see if it varies from beginning or from command line value */
			check_shift(prm, &shift, &ishift, &shift0, record_length1);
		
			if ((verbose) && (n/2000.0 == n/2000)) {
				fprintf(stderr," Working on line %d prf %f record length %d slant_range %d \n"
				,sdr.sequence_number, 0.001*sdr.PRF, record_length1, sdr.slant_range);
			}

			/* read data (and trailing bytes) */
          		if ( fread ((char *) data, record_length1, (size_t) 1, imagefile) != 1 ) break;

			data_length = 2*n_data_burst[*nsub];

			/* write line header to output data  */
          		fwrite((void *) &sdr, line_prefix_size, 1, outfile);
			nlines++;

			/* check to see if this is enough data */
			if(nlines >= ntot) break;

			/* Shimada says the first 13 lines are bad.  I only saw 11 bad. shift these lines outside the window */
			if(kburst < 12) ishift = data_length;
			/* ishift and write the data */
			fill_shift_data(shift, ishift, data_length, line_suffix_size, record_length1, data, shift_data, outfile); 

			/* if kburst it equal to the length of this burst then write the appropriate number of zero lines */
			if(kburst == n_burst[*nsub]) {
				ttot = ttot + dt;
				if(verbose) fprintf(stderr," %d %d %f %f %d %f \n",*nsub+1,kburst,dt,tgap,ngap,ttot);
				kburst = 0;
			/* write the appropriate number of zero lines  */
				for(j=0;j<ngap;j++) {
          				fwrite((void *) &sdr, line_prefix_size, 1, outfile);
					nlines++;
				/* set the gap data to 15.5 */
					for (k=0;k<record_length0;k++) gap_data[k]=NULL_DATA+znew%2;
					fwrite((char *) gap_data, record_length1, 1, outfile); 
				}
			}
		}
		else {
			record_length1 = sdr.record_length - line_prefix_size;
			if ( fread ((char *) data, record_length1, (size_t) 1, imagefile) != 1 ) break;
		}
	}
      
	/* calculate end time and fix prf */
	prm->prf = 0.001*PRF[*nsub];

        prm->clock_stop =  get_clock(sdr, tbias);
        prm->SC_clock_stop = ((double) sdr.sensor_acquisition_year)*1000 + prm->clock_stop;

	/* m is non-zero only in the event of a prf change */
	prm->num_lines = nlines-1;
	prm->num_patches = (int)((1.0*prm->num_lines)/(1.0*prm->num_valid_az));
	if (prm->num_lines == 0) prm->num_lines = 1;

	if (verbose) print_params(prm); 

	free(data);
	free(shift_data);
	fclose (outfile);

	return(*byte_offset);
}
/***************************************************************************/
double get_clock(struct sardata_info sdr, double tbias)
{
double	time;

	//nsd = 24.0*60.0*60.0;	/* seconds in a day */

	time =  (double) sdr.sensor_acquisition_DOY +
		(double) sdr.sensor_acquisition_msecs_day/1000.0/86400.0 +
		tbias/86400.0;

	return(time);
}
/***************************************************************************/
void print_params(struct PRM *prm)
{
	fprintf(stdout,"input_file		= %s \n",prm->input_file);
	fprintf(stdout,"num_rng_bins		= %d \n",prm->num_rng_bins);
	fprintf(stdout,"bytes_per_line		= %d \n",prm->bytes_per_line);
	fprintf(stdout,"good_bytes_per_line	= %d \n",prm->good_bytes);
	fprintf(stdout,"first_sample		= %d \n",prm->first_sample);
	fprintf(stdout,"PRF			= %f \n",prm->prf);
	fprintf(stdout,"pulse_dur		= %e \n",prm->pulsedur);
	fprintf(stdout,"near_range		= %f \n",prm->near_range);
	fprintf(stdout,"num_lines		= %d \n",prm->num_lines);
	fprintf(stdout,"num_patches		= %d \n",prm->num_patches);
       	fprintf(stdout,"SC_clock_start		= %16.10lf \n",prm->SC_clock_start);
       	fprintf(stdout,"SC_clock_stop		= %16.10lf \n",prm->SC_clock_stop);
       	fprintf(stdout,"clock_start		= %16.12lf \n",prm->clock_start);
       	fprintf(stdout,"clock_stop		= %16.12lf \n",prm->clock_stop);
}
/***************************************************************************/
long read_sardata_info(FILE *imagefile, struct PRM *prm, int *header_size, int *line_prefix_size)
{
long nitems;

	*header_size = sizeof(struct sardata_record) + sizeof(struct sardata_descriptor);
	*line_prefix_size = sizeof(struct sardata_info);

	if (*header_size != 720) die("header size is not 720 bytes\n","");
	if (*line_prefix_size != 412) die("header size is not 720 bytes\n","");

	if (debug) fprintf(stderr," header_size %d line_prefix_size %d swap data %d\n", *header_size, *line_prefix_size, swap);

	/* make sure that we are at the beginning */
	/* re-read header even if resetting after a PRF change */
	 rewind(imagefile);

	if (verbose) fprintf(stderr,".... reading header (byte %ld) \n",ftell(imagefile));

	/* data processed before Sept 15, 2006 have a timing bias of 0.9 s */
	/* data processed after this data have a smaller bias 0.0 s */

	nitems = fread((void *) &r1, sizeof(struct sardata_record), 1, imagefile);
	if (debug) { 
		fprintf(stderr,SARDATA_RECORD_WCS,SARDATA_RECORD_RVL(&r1));
		fprintf(stderr," read %ld bytes at position %ld\n", (sizeof(struct sardata_record)), ftell(imagefile));
		}

	nitems = fread((void *) &dfd, sizeof(struct sardata_descriptor), 1, imagefile);
	if (debug) {
		fprintf(stderr,SARDATA_DESCRIPTOR_WCS,SARDATA_DESCRIPTOR_RVL(&dfd));
		fprintf(stderr," read %ld bytes at position %ld\n", (sizeof(struct sardata_descriptor)), ftell(imagefile));
		}

	nitems = fread((void *) &sdr, sizeof(struct sardata_info), 1, imagefile);
	if (debug) fprintf(stderr," read %ld bytes at position %ld\n", (sizeof(struct sardata_info)), ftell(imagefile));

	/* swap data little end/ big end if needed */
	if (swap) swap_ALOS_data_info(&sdr);

	if (debug) fprintf(stderr,SARDATA__WCS,SARDATA_RVL(sdr));

	return(nitems);
}
/***************************************************************************/
int assign_sardata_params(struct PRM *prm, int line_prefix_size, int *line_suffix_size, int *record_length0)
{
double tbias = 0.0, get_clock();

	prm->prf = sdr.PRF;
	prm->pulsedur = (1e-9)*sdr.chirp_length;

	*record_length0 = sdr.record_length - line_prefix_size;

        prm->clock_start =  get_clock(sdr, tbias);
        prm->SC_clock_start = ((double) sdr.sensor_acquisition_year)*1000 + prm->clock_start;

	/* record_length is 21100 */
	/* beginning of line has a 412 byte prefix */
	/* end of line has a 80 byte (40 pixels) suffix (right-fill pixels)*/
	/* record_length0 (data length) is (20688 - 412) = 20276 */
	/* n_data_pixels  10304 */
	/* 2 bytes per pixel */
	/* 412 bytes + (2*10304) bytes + (40*2) bytes  = 21100 bytes*/

	prm->bytes_per_line = sdr.record_length;
   	prm->good_bytes = sdr.record_length; /* make these match the standard FBD data */
	prm->num_rng_bins = 5652; /* make these match the standard FBD data */
	
	*line_suffix_size = sdr.record_length - prm->good_bytes;

//	if (prm->near_range < 0) prm->near_range = sdr.slant_range; 

	if (*record_length0 > 50000) {
		fprintf(stderr, "**** record_length is %d !\n", *record_length0);
		die("expect something like 21100 .... try -swap option?\n","exiting");
		}

	return(EXIT_SUCCESS);
}
/***************************************************************************/
int check_shift(struct PRM *prm, int *shift, int *ishift, int *shift0, int record_length1)
{
        *shift = 2*floor(0.5 + (sdr.slant_range - prm->near_range)/(0.5*SOL/prm->fs));
        *ishift = abs(*shift);

         if (*ishift > record_length1) { 
          	printf(" end: shift exceeds data window %d \n", *shift);
		die("exitting","");
          	}

          if(*shift != *shift0) {
	  	printf(" near_range, shift = %d %d \n", sdr.slant_range, *shift);
            	*shift0 = *shift;
	  	}

	return(EXIT_SUCCESS);
}
/***************************************************************************/
int set_file_position(FILE *imagefile, long *byte_offset, int header_size)
{
	if (*byte_offset < 0) {
		*byte_offset = 0;
		rewind(imagefile);
		fseek(imagefile, header_size, SEEK_SET);
		} else {
		fseek(imagefile, *byte_offset, SEEK_SET);
		}

	return(EXIT_SUCCESS);
}
/***************************************************************************/
int fill_shift_data(int shift, int ishift, int data_length, 
	int line_suffix_size, int record_length1, char *data, char *shift_data, FILE *outfile)
{
int	k;

	/* NULL_DATA = 15; znew randonly is 0 or 1			      */
       	if ( shift > 0) {					
         	for (k = 0; k < ishift; k++) shift_data[k] = NULL_DATA+znew%2;
            	for (k = 0; k < data_length - ishift; k++) shift_data[k + ishift] = data[k];
            	for (k = data_length - ishift; k < record_length1 - ishift;  k++) shift_data[k + ishift] = NULL_DATA+znew%2;
		}

	/* if data is shifted, fill in with data vlues of NULL_DATA at end */
	  if ( shift <= 0) {
            	for (k = 0; k < data_length - ishift - line_suffix_size; k++) shift_data[k] = data[k+ishift];
            	for (k = data_length - ishift - line_suffix_size; k < record_length1; k++ ) shift_data[k] = NULL_DATA+znew%2;
          	}

	/* write the shifted data out */
        fwrite((char *) shift_data, record_length1, 1, outfile);

	return(EXIT_SUCCESS);
}
/***************************************************************************/
int handle_prf_change(struct PRM *prm, FILE *imagefile, long *byte_offset, int n)
{
        prm->num_lines = n;

        /* skip back to beginning of the line */
        fseek(imagefile, -1*sizeof(struct sardata_info), SEEK_CUR);

        /* define byte_offset */
        *byte_offset = ftell(imagefile);

        /* tell the world */
        printf(" *** PRF changed from %lf to  %lf  at line %d (byte %ld)\n", (0.001*prm->prf),(0.001*sdr.PRF), n, *byte_offset);
        printf(" end: PRF changed from %lf to  %lf  at line %d \n", (0.001*prm->prf),(0.001*sdr.PRF), n);

        return(EXIT_SUCCESS);
}
/***************************************************************************/
