/*	$Id$	*/
/************************************************************************
 * esarp reads a file of raw radar echos, and an associated parameter    *
 * 	file and outputs a SAR processed image.  It's subroutines       *
 *	include range compression, Doppler Centroid frequency           *
 *       estimation, range migration, and azimuth compression. This      *
 *    	processor follows one obtained from Howard Zebker as part of the*
 *	Stanford interferometry package.                                *
 ************************************************************************/
/************************************************************************
 * Creator: Evelyn J. Price	(Scripps Institution of Oceanography)   *
 * Date   : 11/18/96						        *
 ************************************************************************/
/************************************************************************
 * Modification History: *
 *				                                        *
 * Date   : May 2003  Hannah Chenh                                       *
 *          - change int *count to int count                             *
 *          - change (void *) &count to &count in fread                  *
 * Date   : Oct 2006 Meng Wei					        *
 *	   - add aastrech() 					        *
 *                                                                       *
 ************************************************************************/
/* delete the doppler estimation part */
/* move the doppler units change right after get params */
/* Apr. 24 2006 */
/* delete nrows from declearation, get it from PRM */
/* use SC_identity to enable reading of different data types */

#include "gmtsar.h"
#include "soi.h"

void print_time(float timer) {
	int min;
	float sec;

	min = (int)timer / 60;
	sec = timer - ((float)min) * 60.0f;

	fprintf(stderr, "Processing Elapsed Time: %d min %.2f sec\n", min, sec);
}

int main(int argc, char *argv[]) {
	unsigned char *indata;
	FILE *fph = NULL, *fpq2 = NULL, *fpi = NULL;
	int ranfft;
	int n, i, j, k, ipatch, low_ind, hi_ind;
	long num_to_seek;
	double delr, rtest, itest, shft, dfact;
	fcomplex **fdata = NULL, *fft_vec = NULL, *ref = NULL;
	scomplex *i2data = NULL;
	int i2test = 1, nclip = 0;
	int pcount = 24, count = 0;
	int ineg = 0;
	void *API = NULL; /* GMT API control structure */

	/* Begin: Initializing new GMT session */
	if ((API = GMT_Create_Session(argv[0], 0U, 0U, NULL)) == NULL)
		return EXIT_FAILURE;

	if (argc < 3) {
		fprintf(stderr, "esarp [GMTSAR] - Produce SAR processed image\n\n");
		fprintf(stderr, "\nUsage: %s filein.PRM fileout.SLC [R4]\n\n", argv[0]);
		exit(-1);
	}

	if ((fph = fopen(argv[1], "r")) == NULL) {
		fprintf(stderr, "Can't open file %s\n", argv[1]);
		exit(-1);
	}

	/* check to see if integer output is desired */
	if (argc == 4) {
		if (argv[3][0] == 'r' || argv[3][0] == 'R') {
			i2test = 0;
		}
	}

	/* get input parameters */
	get_params(fph);
	fclose(fph);
	delr = SOL / fs / 2.0;

	/*  if this is ALOS then increase the scale factor by 2  to improve dynamic
	 * range */
	if (SC_identity == 5 && xmi1 < 16.) {
		dfact = 2. * I2SCALE;
	}
	else if (SC_identity == 8) {
		dfact = 2. * I2SCALE;
	}
	else {
		dfact = I2SCALE;
	}
	printf("I2SCALE %f \n", dfact);

	/* change the units of doppler parameters */
	radopp(&fd1, &fdd1, &fddd1, near_range, delr);

	/* compute range parameters */
	if (num_rng_bins == 0)
		num_rng_bins = good_bytes / 2;
	ranfft = fft_bins(num_rng_bins);
	near_range = near_range + (st_rng_bin - nextend + xshift - 1) * delr;
	far_range = near_range + delr * (num_rng_bins - 1);

	/* Open files. */
	if ((fpi = fopen(input_file, "r")) == NULL) {
		fprintf(stderr, "Can't open file %s\n", input_file);
		exit(-1);
	}
	if ((fpq2 = fopen(argv[2], "w")) == NULL) {
		fprintf(stderr, "Bad output file pointer.\n");
		exit(-1);
	}

	/* allocate memory for and seek in data input file */
	if ((indata = (unsigned char *)malloc(good_bytes * sizeof(unsigned char))) == NULL) {
		fprintf(stderr, "Sorry, couldn't allocate memory for input indata.\n");
		exit(-1);
	}
	if ((fdata = (fcomplex **)malloc(nrows * sizeof(fcomplex))) == NULL) {
		fprintf(stderr, "Sorry, couldn't allocate memory for input data.\n");
		exit(-1);
	}
	for (i = 0; i < nrows; i++) {
		if ((fdata[i] = (fcomplex *)malloc(num_rng_bins * sizeof(fcomplex))) == NULL) {
			fprintf(stderr, "sorry, couldn't allocate memory for input data.\n");
			exit(-1);
		}
	}
	if (i2test) {
		if ((i2data = (scomplex *)malloc(num_rng_bins * sizeof(scomplex))) == NULL) {
			fprintf(stderr, "Sorry, couldn't allocate memory for i2data \n");
			exit(-1);
		}
	}
	if ((fft_vec = (fcomplex *)malloc(ranfft * sizeof(fcomplex))) == NULL) {
		fprintf(stderr, "Sorry, couldn't allocate memory for ref fnc.\n");
		exit(-1);
	}

	/* get the fft of the reference function */
	if ((ref = (fcomplex *)malloc(ranfft * sizeof(fcomplex))) == NULL) {
		fprintf(stderr, "Sorry, couldn't allocate memory for ref fnc.\n");
		exit(-1);
	}
	fprintf(stderr, "Computing range reference function.\n");
	rng_ref(API, ranfft, delr, ref);

	/* check for negative yshift and set ineg accordingly */
	if (yshift < 0) {
		ineg = yshift;
	}

	/* read in data, range compress block, estimate doppler params */
	for (ipatch = 1; ipatch <= num_patches; ipatch++) {

		clock();

		/* seek over IMOP file, first_line, data from last patch */
		num_to_seek = ((long)((ipatch - 1) * num_valid_az + first_line + 1)) * ((long)bytes_per_line);

		if ((n = fseek(fpi, num_to_seek, 0)) != 0) {
			perror(argv[0]);
			exit(-1);
		}

		fprintf(stderr, "Processing patch %d\n", ipatch);
		print_time((float)clock() / CLOCKS_PER_SEC);
		fprintf(stderr, "Range Compression\n");

		if (ineg >= 0)
			if ((n = fseek(fpi, ((long)(bytes_per_line) * ((long)yshift)), 1)) != 0) {
				perror(argv[0]);
				exit(-1);
			}

		for (k = 0; k < nrows; k++) {

			if (ineg >= 0) {
				/* this code reads the number of good data from the header */
				if (SC_identity == 1 || SC_identity == 2 || SC_identity == 5) {
					if ((n = fseek(fpi, (long)(pcount), 1)) != 0) {
						perror(argv[0]);
						exit(-1);
					}
					fread(&count, 4 * sizeof(char), 1, fpi);
					if ((n = fseek(fpi, (long)(first_sample * 2 - pcount - 4), 1)) != 0) {
						perror(argv[0]);
						exit(-1);
					}
				}

				/* this code actually reads the data including trailing bytes that may
				 * be bad */
				fread((void *)indata, 2 * sizeof(unsigned char), good_bytes / 2, fpi);
				fseek(fpi, (long)(bytes_per_line - (first_sample * 2 + good_bytes)), 1);

				/* if the count is bad then use the good_bytes. TXP. */
				if (count == 0)
					count = good_bytes / 2;
				/* now fill the row with good data or zero depending on the SC_identity
				 */
				for (i = 0; i < ranfft; i++) {
					if (i < num_rng_bins) {
						/* ERS-1 and ERS-2 */
						if (SC_identity == 1 || SC_identity == 2) {
							if (i < good_bytes / 2 - first_sample) {
								if ((((int)indata[2 * i]) != 35) && (((int)indata[2 * i + 1]) != 35) && (i < (int)count)) {
									fft_vec[i].r = (float)(indata[2 * i] - xmi1);
									fft_vec[i].i = (float)(indata[2 * i + 1] - xmq1);
								}
								else {
									fft_vec[i].r = 0.0f;
									fft_vec[i].i = 0.0f;
								}
							}
							else {
								fft_vec[i].r = 0.0f;
								fft_vec[i].i = 0.0f;
							}
						}

						/* ENVISAT and ALOS - the number 35 is not used */
						else if (SC_identity == 4 || SC_identity == 5 || SC_identity == 8) {
							if (i < good_bytes / 2 - first_sample) {
								fft_vec[i].r = (float)(indata[2 * i] - xmi1);
								fft_vec[i].i = (float)(indata[2 * i + 1] - xmq1);
							}
							else {
								fft_vec[i].r = 0.0f;
								fft_vec[i].i = 0.0f;
							}
						}

						/* other data formats */
						else {
							fft_vec[i].r = (float)(indata[2 * i] - xmi1);
							fft_vec[i].i = (float)(indata[2 * i + 1] - xmq1);
						}
					}
					else {
						fft_vec[i].r = 0.0f;
						fft_vec[i].i = 0.0f;
					}
				}
			}

			/* use zero lines for a negative yshift */
			else {
				for (i = 0; i < ranfft; i++) {
					fft_vec[i].r = 0.0f;
					fft_vec[i].i = 0.0f;
				}
				ineg = ineg + 1;
			}

			/* range compress line of data */
			rng_cmp(API, ranfft, fft_vec, ref);

			if (xshift >= 0) {
				for (i = xshift; i < num_rng_bins + xshift; i++) {
					fdata[k][i - xshift] = fft_vec[i];
				}
			}
			else {
				for (i = 0; i < num_rng_bins; i++) {
					if (i < (-1 * xshift)) {
						fdata[k][i].r = 0.0f;
						fdata[k][i].i = 0.0f;
					}
					else {
						fdata[k][i] = fft_vec[i + xshift];
					}
				}
			}
		}

		/* transform columns */
		print_time((float)clock() / CLOCKS_PER_SEC);
		fprintf(stderr, "Azimuthal Transform\n");
		trans_col(API, num_rng_bins, nrows, fdata);

		/* range migrate patch */
		print_time((float)clock() / CLOCKS_PER_SEC);
		fprintf(stderr, "Range Migration\n");
		rmpatch(fdata, nrows, delr, fd1, fdd1, fddd1);

		/*azimuth compress patch */
		print_time((float)clock() / CLOCKS_PER_SEC);
		fprintf(stderr, "Azimuthal Compression\n");
		acpatch(API, fdata, nrows, delr, fd1, fdd1, fddd1);

		print_time((float)clock() / CLOCKS_PER_SEC);
		/*apply a azimuth-dependent azimuth shift if a_stretch_a !=0  */
		if (a_stretch_a != 0) {
			fprintf(stderr, "apply azimuth-dependent azimuth shift \n");
			aastretch(fdata, ipatch, nrows, num_valid_az, num_rng_bins, a_stretch_a);
			print_time((float)clock() / CLOCKS_PER_SEC);
		}

		fprintf(stderr, "Writing Data\n");

		low_ind = (nrows - num_valid_az) / 2;
		hi_ind = (nrows + num_valid_az) / 2;

		for (i = low_ind; i < hi_ind; i += nlooks) {

			/* move the data from the 2-D array into a vector to prepare for shift and
			 * output */
			for (j = 0; j < num_rng_bins; j++) {
				fft_vec[j].r = fdata[i][j].r;
				fft_vec[j].i = fdata[i][j].i;
			}

			/* apply a azimuth-dependent range shift if a_stretch_r != 0 */
			if (a_stretch_r != 0.0) {
				for (j = num_rng_bins; j < ranfft; j++) {
					fft_vec[j].r = 0.0f;
					fft_vec[j].i = 0.0f;
				}
				shft = -a_stretch_r * (i - low_ind + num_valid_az * (ipatch - 1));
				shift(API, ranfft, fft_vec, shft);
			}

			/*write data as complex float or complex I2 */
			if (i2test) {
				for (j = 0; j < num_rng_bins; j++) {
					rtest = dfact * fft_vec[j].r;
					itest = dfact * fft_vec[j].i;
					i2data[j].r = (short)clipi2(rtest);
					i2data[j].i = (short)clipi2(itest);

					if ((int)(rtest) > I2MAX || (int)(itest) > I2MAX) {
						nclip++;
					}
				}

				if ((n = (int)fwrite((void *)i2data, 2 * sizeof(short), num_rng_bins, fpq2)) != num_rng_bins) {
					fprintf(stderr, "Problem writing integer data.\n");
				}
			}
			else {
				if ((n = (int)fwrite((void *)fft_vec, 2 * sizeof(float), num_rng_bins, fpq2)) != num_rng_bins) {
					fprintf(stderr, "Problem writing float data.\n");
				}
			}
		}
	}

	print_time((float)clock() / CLOCKS_PER_SEC);
	fclose(fpq2);

	/* close input file */
	fclose(fpi);
	free((char *)ref);
	free((char **)fdata);
	free((unsigned char *)indata);
	free((char *)fft_vec);
	if (i2test) {
		free((char *)i2data);
		fprintf(stderr, "number of points clipped to short int %d \n", nclip);
	}

	if (GMT_Destroy_Session(API))
		return EXIT_FAILURE; /* Remove the GMT machinery */
	return (EXIT_SUCCESS);
}
