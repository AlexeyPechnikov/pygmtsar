/***************************************************************************/
/* read_ALOS_data reads an ALOS IMG file containing SLC data.              */
/* The program skips the first 720 bytes of the IMG file but copies the    */
/* remaining data to the IMG.SLC file after checking and fixing problems.  */
/* The first record is read to determine the linelength, starting PRF,     */
/* and near_range.  If the line length or PRF change then the program      */
/* halts.  If the near_range changes then the lines are shifted and        */
/* unconstrained values at the ends are set to zero.                       */
/* During this processing the available parameters are added to the        */
/* PRM-file.                                                               */
/***************************************************************************/

/***************************************************************************
 * Creator:  David T. Sandwell                                             *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  03/18/2013                                                     *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * DATE                                                                    *
 *                                                                         *
 * 06/29/2006   added the near_range as an optional command-line argument  *
 * 02/19/2007   added the ability to remove duplicate lines		   *
 * 03/07/2007   more robust to bad data in file				   *
 * 03/26/2007   ability to swap bytes (read on PC) RJM			   *
 * 03/28/2007	part of subroutine		   RJM		           *
 * removed set n_azimuth to 9000 rather than default			   *
 * 07/17/08     creates new file when prf changes RJM			   *
 * 07/17/08     reformatted; added functions      RJM                      *
 ***************************************************************************/

/*
the data header information is read into the structure dfd
the line prefix information is read into sdr
Values read here (and or generated) are:
num_rng_bins bytes_per_line good_bytes_per_line
PRF pulse_dur near_range
num_lines num_patches
SC_clock_start SC_clock_stop
clock_start clock_stop
*/

#include "image_sio.h"
#include "lib_functions.h"

void swap32(float *, float *, int);
void swap_ALOS_data_info(struct sardata_info *sdr);
long read_sardata_info(FILE *, struct PRM *, int *, int *);
void print_params(struct PRM *prm);
int assign_sardata_params(struct PRM *, int, int *, int *);
int check_shift(struct PRM *, int *, int *, int *, int);
int set_file_position(FILE *, long *, int);
int reset_params(struct PRM *prm, long *, int *, int *);
int fill_shift_data(int, int, int, int, int, char *, char *, FILE *);
int handle_prf_change(struct PRM *, FILE *, long *, int);

struct sardata_record r1;
struct sardata_descriptor dfd;
struct sardata_info sdr;

long read_ALOS_data_SLC(FILE *imagefile, FILE *outfile, struct PRM *prm, long *byte_offset) {

	float *rdata, *rdata_swap;
	short *i2data;
	int record_length0;        /* length of record read at start of file */
	int record_length1;        /* length of record read in file 	*/
	int start_sdr_rec_len = 0; /* sdr record length for fisrt record */
	int slant_range_old = 0;   /* slant range of previous record */
	int line_suffix_size;      /* number of bytes after data 		*/
	int data_length;           /* bytes of data			*/
	int n, m, ishift, shift, shift0, jj;
	int header_size, line_prefix_size;
	int nclip = 0, nsum = 0;
	double rtest, sgn, rsum = 0., rmad = 0., tfac = 1.;

	double get_clock();

	if (debug)
		fprintf(stderr, ".... reading header \n");

	/* read header information */
	read_sardata_info(imagefile, prm, &header_size, &line_prefix_size);

	/* calculate parameters (data length, range bins, etc) */
	assign_sardata_params(prm, line_prefix_size, &line_suffix_size, &record_length0);

	/* allocate memory */
	if ((rdata = (float *)malloc(record_length0)) == NULL)
		die("couldn't allocate memory for input indata.\n", "");
	if ((rdata_swap = (float *)malloc(record_length0)) == NULL)
		die("couldn't allocate memory for input indata.\n", "");
	if ((i2data = (short *)malloc(record_length0)) == NULL)
		die("couldn't allocate memory for output indata.\n", "");

	// fprintf(stderr,"after allocate data\n");

	/* if byte_offset < 0 this is the first time through 	*/
	/* if prf change has occurred, set file to byte_offset  */
	set_file_position(imagefile, byte_offset, header_size);

	if (verbose)
		fprintf(stderr, ".... reading data (byte %ld) \n", ftell(imagefile));

	shift0 = 0;
	n = 1;
	m = 0;

	/* read the rest of the file */
	while ((fread((void *)&sdr, sizeof(struct sardata_info), 1, imagefile)) == 1) {
		fseek(imagefile, prefix_off, SEEK_CUR); /* skip extra bytes in ALOS-2 prefix */
		n++;

		/* checks for little endian/ big endian */
		if (swap)
			swap_ALOS_data_info(&sdr);

		/* if this is partway through the file due to prf change, reset sequence,
		 * PRF, and near_range */
		if ((*byte_offset > 0) && (n == 2))
			reset_params(prm, byte_offset, &n, &m);

		if (n == 2)
			start_sdr_rec_len = sdr.record_length;
		if (sdr.record_length != start_sdr_rec_len) {
			printf(" ***** warning sdr.record_length error %d \n", sdr.record_length);
			sdr.record_length = start_sdr_rec_len;
			sdr.PRF = prm->prf;
			sdr.slant_range = slant_range_old;
		}
		if (sdr.sequence_number != n)
			printf(" missing line: n, seq# %d %d \n", n, sdr.sequence_number);

		/* check for changes in record_length and PRF */
		record_length1 = sdr.record_length - line_prefix_size - line_suffix_size;
		if (record_length0 != record_length1)
			die("record_length changed", "");

		/* if prf changes, close file and set byte_offset */
		if ((sdr.PRF) != prm->prf) {
			handle_prf_change(prm, imagefile, byte_offset, n);
			break;
		}

		/* check shift to see if it varies from beginning or from command line value
		 */
		check_shift(prm, &shift, &ishift, &shift0, record_length1);

		if ((verbose) && (n / 2000.0 == n / 2000)) {
			fprintf(stderr, " Working on line %d prf %f record length %d slant_range %d \n", sdr.sequence_number, 0.001 * sdr.PRF,
			        record_length1, sdr.slant_range);
		}

		/* read data */
		data_length = record_length1 / 4;
		slant_range_old = sdr.slant_range;
		if (fread((float *)rdata, 4, data_length, imagefile) != data_length)
			break;

		/* skip over suffix */
		fseek(imagefile, line_suffix_size, SEEK_CUR);

		/* write line header to output data  */
		/*fwrite((void *) &sdr, line_prefix_size, 1, outfile);*/

		/* swap the floats if needed and compute some statistics */
		if (swap)
			swap32(rdata, rdata_swap, data_length);
		for (jj = 0; jj < data_length; jj++) {
			sgn = 1.;
			// if(jj%2 != 0) sgn = -1.;  experiment to switch sign of imaginary
			// component
			if (swap) {
				rtest = sgn * slc_fact * rdata_swap[jj];
			}
			else {
				rtest = sgn * slc_fact * rdata[jj];
			}
			/* compute the average value of output */
			rsum = rsum + fabs(rtest);
			nsum = nsum + 1;
			i2data[jj] = (short)clipi2(rtest);
			if ((int)fabs(rtest) > I2MAX)
				nclip = nclip + 1;
		}

		/* write data */
		if (shift == 0) {
			fwrite((float *)i2data, 2, data_length, outfile);
			/* if data is shifted, fill in with data values of NULL_DATA at start or
			 * end*/
		}
		else if (shift != 0) {
			fprintf(stderr, "near range change");
			break;
		}
	}

//printf("chirp_length = %.12d\n",sdr.chirp_length);
//printf("chirp_linear_coeff = %.12d\n",sdr.chirp_linear_coeff);

	/* calculate end time and fix prf */
	prm->prf = 0.001 * prm->prf;

	prm->clock_stop = get_clock(sdr, tbias);
	prm->SC_clock_stop = ((double)sdr.sensor_acquisition_year) * 1000 + prm->clock_stop;

	/* m is non-zero only in the event of a prf change */
	prm->num_lines = n - m - 1;

	/* make sure num lines is divisible by 32 */
	prm->num_lines = prm->num_lines - prm->num_lines % 32;

	/* for SLC data the prm->nrows = prm->num_valid_az = prm->num_lines and
	 * prm->num_patches = 1 */
	prm->nrows = prm->num_lines;
	prm->num_valid_az = prm->num_lines;
	prm->num_patches = 1;

	if (verbose)
		print_params(prm);

	fprintf(stderr, " %d integers were clipped \n", nclip);
	rmad = rsum / nsum;
	tfac = 2000. / rmad;
	if (tfac < 0.333 || tfac > 3.0) {
		fprintf(stderr, " %f median absolute deviation after scaling is \n", rmad);
		fprintf(stderr, " ERROR *** reset SCL_factor to something closer to %f \n", tfac * slc_fact);
	}
	free(rdata);
	free(rdata_swap);
	free(i2data);
	fclose(imagefile);
	fclose(outfile);
	return (*byte_offset);
}
/***************************************************************************/
double get_clock(struct sardata_info sdr, double tbias) {
	double time;

	// nsd = 24.0*60.0*60.0;	/* seconds in a day */

	time = (double)sdr.sensor_acquisition_DOY + (double)sdr.sensor_acquisition_msecs_day / 1000.0 / 86400.0 + tbias / 86400.0;

	return (time);
}
/***************************************************************************/
void print_params(struct PRM *prm) {
	fprintf(stdout, "input_file		= %s \n", prm->input_file);
	fprintf(stdout, "num_rng_bins		= %d \n", prm->num_rng_bins);
	fprintf(stdout, "bytes_per_line		= %d \n", prm->bytes_per_line);
	fprintf(stdout, "good_bytes_per_line	= %d \n", prm->good_bytes);
	fprintf(stdout, "first_sample		= %d \n", prm->first_sample);
	fprintf(stdout, "PRF			= %f \n", prm->prf);
	fprintf(stdout, "pulse_dur		= %e \n", prm->pulsedur);
	fprintf(stdout, "near_range		= %f \n", prm->near_range);
	fprintf(stdout, "num_lines		= %d \n", prm->num_lines);
	fprintf(stdout, "num_patches		= %d \n", prm->num_patches);
	fprintf(stdout, "SC_clock_start		= %16.10lf \n", prm->SC_clock_start);
	fprintf(stdout, "SC_clock_stop		= %16.10lf \n", prm->SC_clock_stop);
	fprintf(stdout, "clock_start             = %16.12lf \n", prm->clock_start);
	fprintf(stdout, "clock_stop              = %16.12lf \n", prm->clock_stop);
}
/***************************************************************************/
long read_sardata_info(FILE *imagefile, struct PRM *prm, int *header_size, int *line_prefix_size) {
	long nitems;

	*header_size = sizeof(struct sardata_record) + sizeof(struct sardata_descriptor);
	*line_prefix_size = sizeof(struct sardata_info) + prefix_off;

	if (*header_size != 720)
		die("header size is not 720 bytes\n", "");
	if (*line_prefix_size != (412 + prefix_off))
		die("prefix size is incorrect \n", "");

	if (debug)
		fprintf(stderr, " header_size %d line_prefix_size %d swap data %d\n", *header_size, *line_prefix_size, swap);

	/* make sure that we are at the beginning */
	/* re-read header even if resetting after a PRF change */
	rewind(imagefile);

	if (verbose)
		fprintf(stderr, ".... reading header (byte %ld) \n", ftell(imagefile));

	/* data processed before Sept 15, 2006 have a timing bias of 0.9 s */
	/* data processed after this data have a smaller bias 0.0 s */

	nitems = fread((void *)&r1, sizeof(struct sardata_record), 1, imagefile);
	if (debug) {
		fprintf(stderr, SARDATA_RECORD_WCS, SARDATA_RECORD_RVL(&r1));
		fprintf(stderr, " read %ld bytes at position %ld\n", (sizeof(struct sardata_record)), ftell(imagefile));
	}

	nitems = fread((void *)&dfd, sizeof(struct sardata_descriptor), 1, imagefile);
	if (debug) {
		fprintf(stderr, SARDATA_DESCRIPTOR_WCS, SARDATA_DESCRIPTOR_RVL(&dfd));
		fprintf(stderr, " read %ld bytes at position %ld\n", (sizeof(struct sardata_descriptor)), ftell(imagefile));
	}

	nitems = fread((void *)&sdr, sizeof(struct sardata_info), 1, imagefile);
	fseek(imagefile, prefix_off, SEEK_CUR); /* skip extra bytes in ALOS-2 prefix */
	if (debug)
		fprintf(stderr, " read %ld bytes at position %ld\n", (sizeof(struct sardata_info) + prefix_off), ftell(imagefile));

	/* swap data little end/ big end if needed */
	if (swap)
		swap_ALOS_data_info(&sdr);

	if (debug)
		fprintf(stderr, SARDATA__WCS, SARDATA_RVL(sdr));

	return (nitems);
}
/***************************************************************************/
int assign_sardata_params(struct PRM *prm, int line_prefix_size, int *line_suffix_size, int *record_length0) {
	double get_clock();

	prm->prf = sdr.PRF;
	prm->pulsedur = (1e-9) * sdr.chirp_length;
    prm->chirp_slope = sdr.chirp_linear_coeff/(1e-6);

	prm->clock_start = get_clock(sdr, tbias);
	prm->SC_clock_start = ((double)sdr.sensor_acquisition_year) * 1000 + prm->clock_start;

	/* record_length is 21100 */
	/* beginning of line has a 412 byte prefix */
	/* end of line has a 80 byte (40 pixels) suffix (right-fill pixels)*/
	/* record_length0 (data length) is (20688 - 412) = 20276 */
	/* n_data_pixels  10304 */
	/* 2 bytes per pixel */
	/* 412 bytes + (2*10304) bytes + (40*2) bytes  = 21100 bytes*/

	prm->good_bytes = 8 * sdr.n_data_pixels + line_prefix_size;
	prm->chirp_ext = 0;
	prm->num_rng_bins = sdr.n_data_pixels + prm->chirp_ext; /* chirp_ext formerly nextend */
	prm->bytes_per_line = sdr.record_length;

	*line_suffix_size = sdr.record_length - prm->good_bytes;
	*record_length0 = sdr.record_length - line_prefix_size - *line_suffix_size;

	/* set to 4 times the number of range bins for PRM file */
	prm->good_bytes = 4 * prm->num_rng_bins;
	prm->bytes_per_line = 4 * prm->num_rng_bins;
	/* set the mean values to zero */
	prm->xmi = 0.;
	prm->xmq = 0.;
	prm->first_sample = 0;

	if (prm->near_range < 0)
		prm->near_range = sdr.slant_range;

	if (*record_length0 > 340000) {
		fprintf(stderr, "**** record_length is %d !\n", *record_length0);
		die("expect something like 340000 .... try -swap option?\n", "exiting");
	}

	return (EXIT_SUCCESS);
}
/***************************************************************************/
int check_shift(struct PRM *prm, int *shift, int *ishift, int *shift0, int record_length1) {
	*shift = 2 * floor(0.5 + (sdr.slant_range - prm->near_range) / (0.5 * SOL / prm->fs));
	*ishift = abs(*shift);

	if (*ishift > record_length1) {
		printf(" end: shift exceeds data window %d \n", *shift);
		die("exitting", "");
	}

	if (*shift != *shift0) {
		printf(" near_range, shift = %d %d \n", sdr.slant_range, *shift);
		*shift0 = *shift;
	}

	return (EXIT_SUCCESS);
}
/***************************************************************************/
int set_file_position(FILE *imagefile, long *byte_offset, int header_size) {
	if (*byte_offset < 0) {
		*byte_offset = 0;
		rewind(imagefile);
		fseek(imagefile, header_size, SEEK_SET);
	}
	else {
		fseek(imagefile, *byte_offset, SEEK_SET);
	}

	return (EXIT_SUCCESS);
}
/***************************************************************************/
int reset_params(struct PRM *prm, long *byte_offset, int *n, int *m) {
	double get_clock();

	prm->clock_start = get_clock(sdr, tbias);
	prm->SC_clock_start = ((double)sdr.sensor_acquisition_year) * 1000 + prm->clock_start;
	prm->prf = sdr.PRF;
	prm->near_range = sdr.slant_range;
	*n = sdr.sequence_number;
	*m = *n;
	*byte_offset = 0;
	if (verbose) {
		fprintf(stderr, " new parameters: \n sequence number %d \n PRF  %f\n near_range  %lf\n", *n, 0.001 * prm->prf,
		        prm->near_range);
	}
	return (EXIT_SUCCESS);
}
/***************************************************************************/
int fill_shift_data(int shift, int ishift, int data_length, int line_suffix_size, int record_length1, char *data,
                    char *shift_data, FILE *outfile) {
	int k;

	/* NULL_DATA = 15; znew randonly is 0 or 1			      */
	if (shift > 0) {
		for (k = 0; k < ishift; k++)
			shift_data[k] = NULL_DATA;
		for (k = 0; k < data_length - ishift; k++)
			shift_data[k + ishift] = data[k];
	}

	/* if data is shifted, fill in with data vlues of NULL_DATA at end */
	if (shift < 0) {
		for (k = 0; k < data_length - ishift - line_suffix_size; k++)
			shift_data[k] = data[k + ishift];
		for (k = data_length - ishift - line_suffix_size; k < record_length1; k++)
			shift_data[k] = NULL_DATA;
	}

	/* write the shifted data out */
	fwrite((char *)shift_data, data_length, 1, outfile);

	return (EXIT_SUCCESS);
}
/***************************************************************************/
int handle_prf_change(struct PRM *prm, FILE *imagefile, long *byte_offset, int n) {
	prm->num_lines = n;

	/* skip back to beginning of the line */
	fseek(imagefile, -1 * (sizeof(struct sardata_info) + prefix_off), SEEK_CUR);

	/* define byte_offset */
	*byte_offset = ftell(imagefile);

	/* tell the world */
	printf(" *** PRF changed from %lf to  %lf  at line %d (byte %ld)\n", (0.001 * prm->prf), (0.001 * sdr.PRF), n, *byte_offset);
	printf(" end: PRF changed from %lf to  %lf  at line %d \n", (0.001 * prm->prf), (0.001 * sdr.PRF), n);

	return (EXIT_SUCCESS);
}
/***************************************************************************/
