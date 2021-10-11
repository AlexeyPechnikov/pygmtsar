/*	$Id: parse_xcorr_input.c 109 2015-01-19 23:01:24Z sandwell $	*/
#include "gmtsar.h"
#include "xcorr.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*-------------------------------------------------------*/
void print_params(struct xcorr *xc) {
	/*
	        fprintf(stdout," input %s \n",xc->param_name);
	*/
	fprintf(stdout, " format %d [0 => complex short, 1 => real float, 2 => real float (NetCDF)] \n", xc->format);
	fprintf(stdout, " m_nx %d \n", xc->m_nx);
	fprintf(stdout, " m_ny %d \n", xc->m_ny);
	fprintf(stdout, " s_nx %d \n", xc->s_nx);
	fprintf(stdout, " s_ny %d \n", xc->s_ny);
	fprintf(stdout, " nxl %d \n", xc->nxl);
	fprintf(stdout, " nyl %d \n", xc->nyl);
	fprintf(stdout, " rshift %d \n", xc->x_offset);
	fprintf(stdout, " ashift %d \n", xc->y_offset);
	fprintf(stdout, " nx_corr %d \n", xc->nx_corr);
	fprintf(stdout, " ny_corr %d \n", xc->ny_corr);
	fprintf(stdout, " yoff %d \n", xc->ysearch);
	fprintf(stdout, " xoff %d \n", xc->xsearch);
	fprintf(stdout, " npx %d \n", xc->npx);
	fprintf(stdout, " npy %d \n", xc->npy);
	fprintf(stdout, " npx %d \n", xc->npx);
	fprintf(stdout, " npy %d \n", xc->npy);
	fprintf(stderr, "data file 1 %s \n", xc->data1_name);
	fprintf(stderr, "data file 2 %s \n", xc->data2_name);
}
/*-------------------------------------------------------*/
void set_defaults(struct xcorr *xc) {
	xc->format = 0; /* data format */
	xc->ri = 2;     /* range interpolation factor */

	/* default values for time correlation */
	if ((xc->corr_flag == 0) || (xc->corr_flag == 1)) {
		xc->nxl = 16;     /* number of points in x direction */
		xc->nyl = 32;     /* number of points in y direction */
		xc->x_offset = 0; /* starting x offset point*/
		xc->y_offset = 0; /* starting y offset point*/

		/* 128; 64 */
		xc->nx_corr = 128; /* size of correlation window */
		xc->ny_corr = 128; /* size of correlation window */

		xc->xsearch = 64; /* size of search x offset */
		xc->ysearch = 64; /* size of search y offset */
	}

	/* default values for fft correlation 		*/
	/* xc->npx and xc->npx must be powers of two	*/
	if (xc->corr_flag == 2) {
		xc->nxl = 16;     /* number of points in x direction */
		xc->nyl = 32;     /* number of points in y direction */
		xc->x_offset = 0; /* starting x offset point*/
		xc->y_offset = 0; /* starting y offset point*/

		/* 128; 64 */
		xc->nx_corr = 128; /* size of correlation window */
		xc->ny_corr = 128; /* size of correlation window */

		xc->xsearch = 64; /* size of search x offset */
		xc->ysearch = 64; /* size of search y offset */
	}

	/* default values for fft interpolation 		*/
	xc->n2x = 8;
	xc->n2y = 8;

	xc->npx = xc->nx_corr + 2 * xc->xsearch; /* size of data required*/
	xc->npy = xc->ny_corr + 2 * xc->ysearch; /* size of data required*/

	xc->nxc = 2 * xc->xsearch; /* size of correlation patch*/
	xc->nyc = 2 * xc->ysearch; /* size of correlation patch*/

	xc->astretcha = 0.0;    /* azimuth stretch partameter */
	xc->interp_flag = 1;    /* interpolate or not ? */
	xc->interp_factor = 16; /* interpolatation factor */
}
/*-------------------------------------------------------*/
/* reads options 				*/
/* start with third arguement			*/

void parse_command_line(int na, char **a, struct xcorr *xc, int *nfiles, int *input_flag, char *USAGE) {
	int n;
	FILE *inputfile;
	char tmp[128];

	for (n = 3; n < na; n++) {
		if (!strcmp(a[n], "-freq")) {
			xc->corr_flag = 2;
			fprintf(stderr, " using frequency cross correlation\n");
		}
		else if (!strcmp(a[n], "-time")) {
			xc->corr_flag = 0;
			fprintf(stderr, " using time cross correlation\n");
		}
		else if (!strcmp(a[n], "-time4")) {
			xc->corr_flag = 1;
			fprintf(stderr, " using 4th power time cross correlation\n");
		}
		else if (!strcmp(a[n], "-real")) {
			xc->format = 1;
			fprintf(stderr, " no interpolation\n");
		}
		else if (!strcmp(a[n], "-nointerp")) {
			xc->interp_flag = 0;
			fprintf(stderr, " no interpolation\n");
		}
		else if (!strcmp(a[n], "-noshift")) {
			xc->offset_flag = 1;
			fprintf(stderr, " ignoring shift in prm\n");
		}
		else if (!strcmp(a[n], "-interp")) {
			n++;
			if (n == na)
				die(" no option after -interp!\n", "");
			xc->interp_flag = 1;
			xc->interp_factor = atoi(a[n]);
			fprintf(stderr, " setting interpolation factor to %d\n", xc->interp_factor);
		}
		else if (!strcmp(a[n], "-range_interp")) {
			n++;
			if (n == na)
				die(" no option after -range_interp!\n", "");
			xc->ri = atoi(a[n]);
			fprintf(stderr, " setting range interpolation factor to %d\n", xc->ri);
		}
		else if (!strcmp(a[n], "-nx")) {
			n++;
			if (n == na)
				die(" no option after -nx!\n", "");
			xc->nxl = atoi(a[n]);
			fprintf(stderr, " setting nx to %d\n", xc->nxl);
		}
		else if (!strcmp(a[n], "-ny")) {
			n++;
			if (n == na)
				die(" no option after -ny!\n", "");
			xc->nyl = atoi(a[n]);
			fprintf(stderr, " setting ny to %d\n", xc->nyl);
		}
		else if (!strcmp(a[n], "-xsearch")) {
			n++;
			if (n == na)
				die(" no option after -xsearch!\n", "");
			xc->xsearch = atoi(a[n]);
			xc->nx_corr = 2 * xc->xsearch;
			xc->npx = xc->nx_corr + 2 * xc->xsearch; /* size of data required*/
			xc->nxc = 2 * xc->xsearch;               /* size of correlation patch*/
			if (((xc->xsearch - 1) & xc->xsearch))
				die(" xsearch needs to be power of 2! (32 64 128 256) \n", "");
			fprintf(stderr, " setting xsearch to %d\n", xc->xsearch);
			fprintf(stderr, " setting nx_corr to %d\n", xc->nx_corr);
		}
		else if (!strcmp(a[n], "-ysearch")) {
			n++;
			if (n == na)
				die(" no option after -ysearch!\n", "");
			xc->ysearch = atoi(a[n]);
			xc->ny_corr = 2 * xc->ysearch;
			xc->npy = xc->ny_corr + 2 * xc->ysearch; /* size of data required*/
			xc->nyc = 2 * xc->ysearch;               /* size of correlation patch*/
			if (((xc->ysearch - 1) & xc->ysearch))
				die(" ysearch needs to be power of 2! (32 64 128 256) \n", "");
			fprintf(stderr, " setting ysearch to %d\n", xc->ysearch);
			fprintf(stderr, " setting ny_corr to %d\n", xc->ny_corr);
		}
		else if (!strcmp(a[n], "-v")) {
			verbose = 1;
			fprintf(stderr, " verbose output \n");
		}
		else if (!strcmp(a[n], "-norange")) {
			xc->ri = 1;
			fprintf(stderr, " no range interpolation \n");
		}
		else if (!strncmp(a[n], "-input", 1)) {
			n++;
			if (n == na)
				die(" no option after -input!\n", "");
			fprintf(stderr, "using input file \n");
			*input_flag = 1;

			if ((inputfile = fopen(a[2], "r")) == NULL)
				die("Can't open ", a[2]);

			while (fscanf(inputfile, " %s ", tmp) != EOF)
				(*nfiles)++;
			fclose(inputfile);
		}
		else {
			fprintf(stderr, " %s *** option not recognized ***\n\n", a[n]);
			fprintf(stderr, " %s someone made a mistake!\n\n", a[n]);
			fprintf(stderr, " %s I think it was you! \n\n", a[n]);
			fprintf(stderr, "%s", USAGE);
			exit(1);
		}
	}
}
/*---------------------------------------------------------------------------*/
void handle_prm(void *API, char **argv, struct xcorr *xc, int nfiles) {

	int i;
	char **filename;
	FILE *prmfile;
	struct PRM *r;
	char tmp_c[1024];

	if (debug)
		fprintf(stderr, "handle_prm %d\n", nfiles);

	filename = malloc(nfiles * sizeof(char *));
	for (i = 0; i < nfiles; i++)
		filename[i] = malloc(128 * sizeof(char));

	r = malloc(nfiles * sizeof(struct PRM));

	for (i = 0; i < nfiles; i++) {
		strcpy(filename[i], argv[i + 1]);
		strcpy(tmp_c, &filename[i][strlen(filename[i]) - 4]);
		// fprintf(stderr,"%s\n",tmp_c);
		if (strcmp(tmp_c, ".PRM") == 0) {

			// fprintf(stderr," Reading in PRM file: %s\n",filename[i]);
			if ((prmfile = fopen(filename[i], "r")) == NULL)
				die("Can't open prmfile ", filename[i]);
            null_sio_struct(&r[i]);
			get_sio_struct(prmfile, &r[i]);

			if (i == 0) {
				strcpy(xc->data1_name, r[i].SLC_file);
				if ((xc->data1 = fopen(xc->data1_name, "r")) == NULL)
					die("Cannot open SLC_file", xc->data1_name);
				xc->m_nx = r[i].num_rng_bins;
				xc->m_ny = r[i].num_patches * r[i].num_valid_az;
			}

			if (i == 1) {
				strcpy(xc->data2_name, r[i].SLC_file);
				xc->data2 = fopen(xc->data2_name, "r");
				if ((xc->data2 = fopen(xc->data2_name, "r")) == NULL)
					die("Cannot open SLC_file", xc->data2_name);
				xc->s_nx = r[i].num_rng_bins;
				xc->s_ny = r[i].num_patches * r[i].num_valid_az;
				xc->x_offset = r[i].rshift;
				xc->y_offset = r[i].ashift;
			}
			/* added rjm 5/25/2010 to avoid div by zero and Nan */
			if (r[0].prf > 0) {
				xc->astretcha = (r[1].prf - r[0].prf) / r[0].prf;
			}
			else {
				xc->astretcha = 0.0;
			}

			if (xc->offset_flag == 1) {
				xc->x_offset = xc->y_offset = 0;
				fprintf(stderr, " setting ashift and rshift to zero\n");
			}

			fclose(prmfile);
		}
		else {
			// fprintf(stderr," Reading in netcdf file: %s\n",filename[i]);
			xc->format = 2;

			if (i == 0) {
				if ((xc->D1 = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL,
				                            filename[i], NULL)) == NULL)
					die("cannot open topofile", filename[i]);
				strcpy(xc->data1_name, filename[i]);
				xc->m_nx = xc->D1->header->n_columns;
				xc->m_ny = xc->D1->header->n_rows;
				if ((GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, filename[i],
				                   xc->D1)) == NULL)
					die("cannot open topofile", filename[i]);
			}
			if (i == 1) {
				if ((xc->D2 = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL,
				                            filename[i], NULL)) == NULL)
					die("cannot open topofile", filename[i]);
				strcpy(xc->data2_name, filename[i]);
				xc->s_nx = xc->D2->header->n_columns;
				xc->s_ny = xc->D2->header->n_rows;
				xc->x_offset = 0;
				xc->y_offset = 0;
				if ((GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, filename[i],
				                   xc->D2)) == NULL)
					die("cannot open topofile", filename[i]);
			}
		}
	}
	fprintf(stderr, " %d %d %d %d %d %d %f\n", xc->m_nx, xc->m_ny, xc->s_nx, xc->s_ny, xc->x_offset, xc->y_offset, xc->astretcha);

	free(r);
}
/*---------------------------------------------------------------------------*/
