/**************************************************************************/
/*                                                                        */
/*   Mosaic 5 subswaths of 2 frames of ALOS ScanSAR in radar coordinates  */
/*                                                                        */
/*   This program is written for better unwrapping because                */
/*   there will be phase difference for each subswath if unwrapping       */
/*   individually and the phase difference might not be a constant.       */
/*                                                                        */
/*   This program will read the clock_start and near_range info from      */
/*   5 PRM files of each subswath and read the pixel values from each     */
/*   GMT grids, calculate the coordinates of each subswath in meters,     */
/*   put each subswath in the corresponding coordinates. The boundary in  */
/*   range direction of each subswath is chosen to be at the middle of    */
/*   the overlapping region.                                              */
/*                                                                        */
/*   After this step, should be phase unwrapping.                         */
/*                                                                        */
/*   May 18 2018                                                          */
/*                                                                        */
/*   Xiaopeng Tong                                                        */
/*                                                                        */
/*   Update history:                                                      */
/*                                                                        */
/*     May 25 2018    modified to mosaic 2 frames                         */
/*     May 30 2018    can divide the mosaic into 5 subswaths              */
/*                                                                        */
/**************************************************************************/

char *USAGE = " ALOS_mosaic_ss_2frames prmfile1 grd1 prmfile2 grd2 prmfile3 grd3 "
              "prmfile4 grd4 prmfile5 grd5 prmfile6 grd6 prmfile7 grd7 prmfile8 grd8 "
              "prmfile9 grd9 prmfile10 grd10 mosaic dec dir \n\n"
              " USAGE:  \n"
              "  mosaic 5 subswath of 2 frames of ALOS ScanSAR in radar coordinates\n\n"
              "  prmfile[1-10]    the PRM files \n"
              "  grd[1-10]        the GMT grid files (corr, phase, phasefilt) \n"
              "  dec              decimation factor in radar coordinates (100) \n"
              "  dir              dir is 1 means mosaic 5 subswaths \n"
              "                   dir is 0 means divide the mosaic into 5 subswaths \n"
              "                   the first frame is 1-5 and the second frame is 6-10 \n";

#include "gmt.h"
#include "gmtsar.h"

void set_ALOS_defaults(struct PRM *);
void get_sio_struct(FILE *, struct PRM *);

int main(int argc, char **argv) {

	/* define variables */
	FILE *prmfile;
	struct PRM prm;
	void *API = NULL; /* GMT control structure */
	struct GMT_GRID *grid1 = NULL, *grid2 = NULL, *grid3 = NULL, *grid4 = NULL, *grid5 = NULL, *grid6 = NULL, *grid7 = NULL,
	                *grid8 = NULL, *grid9 = NULL, *grid10 = NULL, *grid11 = NULL;
	double inc[2], wesn[4];
	double clock_start[10], clock_start_rel[10], near_range[10], near_range_rel[10];
	double prf[10], rng_samp_rate[10], vel[10];
	int num_rng[10], num_azi[10];
	double clock_start_min;
	int ra[10][4];
	int range_ub, azimuth_ub;
	double sol = 3.0e8;
	double rsr, dr;
	double dec;
	int row, col;
	int i;
	int range, azimuth;
	int range_boundary1, range_boundary2, range_boundary3, range_boundary4, range_boundary5;
	int azimuth_boundary1, azimuth_boundary2, azimuth_boundary3, azimuth_boundary4, azimuth_boundary5;
	int node1, node2, node3, node4, node5, node6, node7, node8, node9, node10, node11;
	int dir;
	double num_x, num_y;

	/* parse arguments */
	if (argc != 24)
		die("\n", USAGE);

	/* Begin: Initializing new GMT session */
	if ((API = GMT_Create_Session(argv[0], 2, 0, NULL)) == NULL)
		return EXIT_FAILURE;

	/* read in dec */
	dec = atof(argv[22]);

	/* read in dir */
	dir = atoi(argv[23]);

	num_x = 1000;
	num_y = 5000;

	/* open and read PRM files */
	fprintf(stderr, "input:\n");

	/* subswath 1 */
	if ((prmfile = fopen(argv[1], "r")) == NULL) {
		die("\n can't open PRM file:\n", argv[1]);
	}
	set_ALOS_defaults(&prm);
	get_sio_struct(prmfile, &prm);
	clock_start[0] = 86400.0 * prm.clock_start + (prm.nrows - prm.num_valid_az) / (2. * prm.prf);
	near_range[0] = prm.near_range;
	rng_samp_rate[0] = prm.fs;
	num_rng[0] = prm.num_rng_bins;
	num_azi[0] = prm.num_patches * prm.num_valid_az;
	prf[0] = prm.prf;
	vel[0] = prm.vel;
	fclose(prmfile);

	/* subswath 2 */
	if ((prmfile = fopen(argv[3], "r")) == NULL) {
		die("\n can't open PRM file:\n", argv[3]);
	}
	set_ALOS_defaults(&prm);
	get_sio_struct(prmfile, &prm);
	clock_start[1] = 86400.0 * prm.clock_start + (prm.nrows - prm.num_valid_az) / (2. * prm.prf);
	near_range[1] = prm.near_range;
	rng_samp_rate[1] = prm.fs;
	num_rng[1] = prm.num_rng_bins;
	num_azi[1] = prm.num_patches * prm.num_valid_az;
	prf[1] = prm.prf;
	vel[1] = prm.vel;
	fclose(prmfile);

	/* subswath 3 */
	if ((prmfile = fopen(argv[5], "r")) == NULL) {
		die("\n can't open PRM file:\n", argv[5]);
	}
	set_ALOS_defaults(&prm);
	get_sio_struct(prmfile, &prm);
	clock_start[2] = 86400.0 * prm.clock_start + (prm.nrows - prm.num_valid_az) / (2. * prm.prf);
	near_range[2] = prm.near_range;
	rng_samp_rate[2] = prm.fs;
	num_rng[2] = prm.num_rng_bins;
	num_azi[2] = prm.num_patches * prm.num_valid_az;
	prf[2] = prm.prf;
	vel[2] = prm.vel;
	fclose(prmfile);

	/* subswath 4 */
	if ((prmfile = fopen(argv[7], "r")) == NULL) {
		die("\n can't open PRM file:\n", argv[7]);
	}
	set_ALOS_defaults(&prm);
	get_sio_struct(prmfile, &prm);
	clock_start[3] = 86400.0 * prm.clock_start + (prm.nrows - prm.num_valid_az) / (2. * prm.prf);
	near_range[3] = prm.near_range;
	rng_samp_rate[3] = prm.fs;
	num_rng[3] = prm.num_rng_bins;
	num_azi[3] = prm.num_patches * prm.num_valid_az;
	prf[3] = prm.prf;
	vel[3] = prm.vel;
	fclose(prmfile);

	/* subswath 5 */
	if ((prmfile = fopen(argv[9], "r")) == NULL) {
		die("\n can't open PRM file:\n", argv[9]);
	}
	set_ALOS_defaults(&prm);
	get_sio_struct(prmfile, &prm);
	clock_start[4] = 86400.0 * prm.clock_start + (prm.nrows - prm.num_valid_az) / (2. * prm.prf);
	near_range[4] = prm.near_range;
	rng_samp_rate[4] = prm.fs;
	num_rng[4] = prm.num_rng_bins;
	num_azi[4] = prm.num_patches * prm.num_valid_az;
	prf[4] = prm.prf;
	vel[4] = prm.vel;
	fclose(prmfile);

	/* subswath 6 */
	if ((prmfile = fopen(argv[11], "r")) == NULL) {
		die("\n can't open PRM file:\n", argv[11]);
	}
	set_ALOS_defaults(&prm);
	get_sio_struct(prmfile, &prm);
	clock_start[5] = 86400.0 * prm.clock_start + (prm.nrows - prm.num_valid_az) / (2. * prm.prf);
	near_range[5] = prm.near_range;
	rng_samp_rate[5] = prm.fs;
	num_rng[5] = prm.num_rng_bins;
	num_azi[5] = prm.num_patches * prm.num_valid_az;
	prf[5] = prm.prf;
	vel[5] = prm.vel;
	fclose(prmfile);

	/* subswath 7 */
	if ((prmfile = fopen(argv[13], "r")) == NULL) {
		die("\n can't open PRM file:\n", argv[13]);
	}
	set_ALOS_defaults(&prm);
	get_sio_struct(prmfile, &prm);
	clock_start[6] = 86400.0 * prm.clock_start + (prm.nrows - prm.num_valid_az) / (2. * prm.prf);
	near_range[6] = prm.near_range;
	rng_samp_rate[6] = prm.fs;
	num_rng[6] = prm.num_rng_bins;
	num_azi[6] = prm.num_patches * prm.num_valid_az;
	prf[6] = prm.prf;
	vel[6] = prm.vel;
	fclose(prmfile);

	/* subswath 8 */
	if ((prmfile = fopen(argv[15], "r")) == NULL) {
		die("\n can't open PRM file:\n", argv[15]);
	}
	set_ALOS_defaults(&prm);
	get_sio_struct(prmfile, &prm);
	clock_start[7] = 86400.0 * prm.clock_start + (prm.nrows - prm.num_valid_az) / (2. * prm.prf);
	near_range[7] = prm.near_range;
	rng_samp_rate[7] = prm.fs;
	num_rng[7] = prm.num_rng_bins;
	num_azi[7] = prm.num_patches * prm.num_valid_az;
	prf[7] = prm.prf;
	vel[7] = prm.vel;
	fclose(prmfile);

	/* subswath 9 */
	if ((prmfile = fopen(argv[17], "r")) == NULL) {
		die("\n can't open PRM file:\n", argv[17]);
	}
	set_ALOS_defaults(&prm);
	get_sio_struct(prmfile, &prm);
	clock_start[8] = 86400.0 * prm.clock_start + (prm.nrows - prm.num_valid_az) / (2. * prm.prf);
	near_range[8] = prm.near_range;
	rng_samp_rate[8] = prm.fs;
	num_rng[8] = prm.num_rng_bins;
	num_azi[8] = prm.num_patches * prm.num_valid_az;
	prf[8] = prm.prf;
	vel[8] = prm.vel;
	fclose(prmfile);

	/* subswath 10 */
	if ((prmfile = fopen(argv[19], "r")) == NULL) {
		die("\n can't open PRM file:\n", argv[19]);
	}
	set_ALOS_defaults(&prm);
	get_sio_struct(prmfile, &prm);
	clock_start[9] = 86400.0 * prm.clock_start + (prm.nrows - prm.num_valid_az) / (2. * prm.prf);
	near_range[9] = prm.near_range;
	rng_samp_rate[9] = prm.fs;
	num_rng[9] = prm.num_rng_bins;
	num_azi[9] = prm.num_patches * prm.num_valid_az;
	prf[9] = prm.prf;
	vel[9] = prm.vel;
	fclose(prmfile);

	/* open and read GMT grids */
	if (dir == 1) {
		/* subswath 1 */
		if ((grid1 = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, argv[2], NULL)) ==
		    NULL)
			return EXIT_FAILURE;
		/* check if the dimension of the grid and PRM file match */
		fprintf(stderr, "subswath 1: num_rng=%d, num_azi=%d\n", num_rng[0], num_azi[0]);
		fprintf(stderr,
		        "gmt grid nx=%d, gmt grid inc_x=%f, gmt grid ny=%d, gmt grid "
		        "inc_y=%f\n",
		        grid1->header->n_columns, grid1->header->inc[GMT_X], grid1->header->n_rows, grid1->header->inc[GMT_Y]);

		if ((num_rng[0] != round(grid1->header->n_columns * grid1->header->inc[GMT_X])) ||
		    (num_azi[0] != round(grid1->header->n_rows * grid1->header->inc[GMT_Y]))) {
			die("\nPRM file and GMT grid don't match!\n", argv[2]);
		}
		if (GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, argv[2], grid1) == NULL)
			return EXIT_FAILURE;

		/* subswath 2 */
		if ((grid2 = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, argv[4], NULL)) ==
		    NULL)
			return EXIT_FAILURE;
		/* check if the dimension of the grid and PRM file match */
		fprintf(stderr, "subswath 2: num_rng=%d, num_azi=%d\n", num_rng[1], num_azi[1]);
		fprintf(stderr,
		        "gmt grid nx=%d, gmt grid inc_x=%f, gmt grid ny=%d, gmt grid "
		        "inc_y=%f\n",
		        grid2->header->n_columns, grid2->header->inc[GMT_X], grid2->header->n_rows, grid2->header->inc[GMT_Y]);
		if ((num_rng[1] != round(grid2->header->n_columns * grid2->header->inc[GMT_X])) ||
		    (num_azi[1] != round(grid2->header->n_rows * grid2->header->inc[GMT_Y]))) {
			die("\nPRM file and GMT grid don't match!\n", argv[4]);
		}
		if (GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, argv[4], grid2) == NULL)
			return EXIT_FAILURE;

		/* subswath 3 */
		if ((grid3 = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, argv[6], NULL)) ==
		    NULL)
			return EXIT_FAILURE;
		/* check if the dimension of the grid and PRM file match */
		fprintf(stderr, "subswath 3: num_rng=%d, num_azi=%d\n", num_rng[2], num_azi[2]);
		fprintf(stderr,
		        "gmt grid nx=%d, gmt grid inc_x=%f, gmt grid ny=%d, gmt grid "
		        "inc_y=%f\n",
		        grid3->header->n_columns, grid3->header->inc[GMT_X], grid3->header->n_rows, grid3->header->inc[GMT_Y]);
		if ((num_rng[2] != round(grid3->header->n_columns * grid3->header->inc[GMT_X])) ||
		    (num_azi[2] != round(grid3->header->n_rows * grid3->header->inc[GMT_Y]))) {
			die("\nPRM file and GMT grid don't match!\n", argv[6]);
		}
		if (GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, argv[6], grid3) == NULL)
			return EXIT_FAILURE;

		/* subswath 4 */
		if ((grid4 = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, argv[8], NULL)) ==
		    NULL)
			return EXIT_FAILURE;
		/* check if the dimension of the grid and PRM file match */
		fprintf(stderr, "subswath 4: num_rng=%d, num_azi=%d\n", num_rng[3], num_azi[3]);
		fprintf(stderr,
		        "gmt grid nx=%d, gmt grid inc_x=%f, gmt grid ny=%d, gmt grid "
		        "inc_y=%f\n",
		        grid4->header->n_columns, grid4->header->inc[GMT_X], grid4->header->n_rows, grid4->header->inc[GMT_Y]);
		if ((num_rng[3] != round(grid4->header->n_columns * grid4->header->inc[GMT_X])) ||
		    (num_azi[3] != round(grid4->header->n_rows * grid4->header->inc[GMT_Y]))) {
			die("\nPRM file and GMT grid don't match!\n", argv[8]);
		}
		if (GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, argv[8], grid4) == NULL)
			return EXIT_FAILURE;

		/* subswath 5 */
		if ((grid5 = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, argv[10], NULL)) ==
		    NULL)
			return EXIT_FAILURE;
		/* check if the dimension of the grid and PRM file match */
		fprintf(stderr, "subswath 5: num_rng=%d, num_azi=%d\n", num_rng[4], num_azi[4]);
		fprintf(stderr,
		        "gmt grid nx=%d, gmt grid inc_x=%f, gmt grid ny=%d, gmt grid "
		        "inc_y=%f\n",
		        grid5->header->n_columns, grid5->header->inc[GMT_X], grid5->header->n_rows, grid5->header->inc[GMT_Y]);
		if ((num_rng[4] != round(grid5->header->n_columns * grid5->header->inc[GMT_X])) ||
		    (num_azi[4] != round(grid5->header->n_rows * grid5->header->inc[GMT_Y]))) {
			die("\nPRM file and GMT grid don't match!\n", argv[10]);
		}
		if (GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, argv[10], grid5) == NULL)
			return EXIT_FAILURE;

		/* subswath 6 */
		if ((grid6 = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, argv[12], NULL)) ==
		    NULL)
			return EXIT_FAILURE;
		/* check if the dimension of the grid and PRM file match */
		fprintf(stderr, "subswath 6: num_rng=%d, num_azi=%d\n", num_rng[5], num_azi[5]);
		fprintf(stderr,
		        "gmt grid nx=%d, gmt grid inc_x=%f, gmt grid ny=%d, gmt grid "
		        "inc_y=%f\n",
		        grid6->header->n_columns, grid6->header->inc[GMT_X], grid6->header->n_rows, grid6->header->inc[GMT_Y]);
		if ((num_rng[5] != round(grid6->header->n_columns * grid6->header->inc[GMT_X])) ||
		    (num_azi[5] != round(grid6->header->n_rows * grid6->header->inc[GMT_Y]))) {
			die("\nPRM file and GMT grid don't match!\n", argv[12]);
		}
		if (GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, argv[12], grid6) == NULL)
			return EXIT_FAILURE;

		/* subswath 7 */
		if ((grid7 = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, argv[14], NULL)) ==
		    NULL)
			return EXIT_FAILURE;
		/* check if the dimension of the grid and PRM file match */
		fprintf(stderr, "subswath 7: num_rng=%d, num_azi=%d\n", num_rng[6], num_azi[6]);
		fprintf(stderr,
		        "gmt grid nx=%d, gmt grid inc_x=%f, gmt grid ny=%d, gmt grid "
		        "inc_y=%f\n",
		        grid7->header->n_columns, grid7->header->inc[GMT_X], grid7->header->n_rows, grid7->header->inc[GMT_Y]);
		if ((num_rng[6] != round(grid7->header->n_columns * grid7->header->inc[GMT_X])) ||
		    (num_azi[6] != round(grid7->header->n_rows * grid7->header->inc[GMT_Y]))) {
			die("\nPRM file and GMT grid don't match!\n", argv[14]);
		}
		if (GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, argv[14], grid7) == NULL)
			return EXIT_FAILURE;

		/* subswath 8 */
		if ((grid8 = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, argv[16], NULL)) ==
		    NULL)
			return EXIT_FAILURE;
		/* check if the dimension of the grid and PRM file match */
		fprintf(stderr, "subswath 8: num_rng=%d, num_azi=%d\n", num_rng[7], num_azi[7]);
		fprintf(stderr,
		        "gmt grid nx=%d, gmt grid inc_x=%f, gmt grid ny=%d, gmt grid "
		        "inc_y=%f\n",
		        grid8->header->n_columns, grid8->header->inc[GMT_X], grid8->header->n_rows, grid8->header->inc[GMT_Y]);
		if ((num_rng[7] != round(grid8->header->n_columns * grid8->header->inc[GMT_X])) ||
		    (num_azi[7] != round(grid8->header->n_rows * grid8->header->inc[GMT_Y]))) {
			die("\nPRM file and GMT grid don't match!\n", argv[16]);
		}
		if (GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, argv[16], grid8) == NULL)
			return EXIT_FAILURE;

		/* subswath 9 */
		if ((grid9 = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, argv[18], NULL)) ==
		    NULL)
			return EXIT_FAILURE;
		/* check if the dimension of the grid and PRM file match */
		fprintf(stderr, "subswath 9: num_rng=%d, num_azi=%d\n", num_rng[8], num_azi[8]);
		fprintf(stderr,
		        "gmt grid nx=%d, gmt grid inc_x=%f, gmt grid ny=%d, gmt grid "
		        "inc_y=%f\n",
		        grid9->header->n_columns, grid9->header->inc[GMT_X], grid9->header->n_rows, grid9->header->inc[GMT_Y]);
		if ((num_rng[8] != round(grid9->header->n_columns * grid9->header->inc[GMT_X])) ||
		    (num_azi[8] != round(grid9->header->n_rows * grid9->header->inc[GMT_Y]))) {
			die("\nPRM file and GMT grid don't match!\n", argv[18]);
		}
		if (GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, argv[18], grid9) == NULL)
			return EXIT_FAILURE;

		/* subswath 10 */
		if ((grid10 = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, argv[20], NULL)) ==
		    NULL)
			return EXIT_FAILURE;
		/* check if the dimension of the grid and PRM file match */
		fprintf(stderr, "subswath 10: num_rng=%d, num_azi=%d\n", num_rng[9], num_azi[9]);
		fprintf(stderr,
		        "gmt grid nx=%d, gmt grid inc_x=%f, gmt grid ny=%d, gmt grid "
		        "inc_y=%f\n",
		        grid10->header->n_columns, grid10->header->inc[GMT_X], grid10->header->n_rows, grid10->header->inc[GMT_Y]);
		if ((num_rng[9] != round(grid10->header->n_columns * grid10->header->inc[GMT_X])) ||
		    (num_azi[9] != round(grid10->header->n_rows * grid10->header->inc[GMT_Y]))) {
			die("\nPRM file and GMT grid don't match!\n", argv[20]);
		}
		if (GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, argv[20], grid10) == NULL)
			return EXIT_FAILURE;
	}
	else if (dir == 0) {
		/* subswath 1 */
		wesn[GMT_XLO] = 0;
		wesn[GMT_YLO] = 0;
		wesn[GMT_XHI] = num_rng[0];
		wesn[GMT_YHI] = num_azi[0];
		inc[GMT_X] = (wesn[GMT_XHI] - wesn[GMT_XLO]) / num_x;
		inc[GMT_Y] = (wesn[GMT_YHI] - wesn[GMT_YLO]) / num_y;

		fprintf(stderr, "wesn = %f %f %f %f \n", wesn[GMT_XLO], wesn[GMT_YLO], wesn[GMT_XHI], wesn[GMT_YHI]);
		if ((grid1 = GMT_Create_Data(API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, wesn, inc, GMT_GRID_PIXEL_REG, 0,
		                             NULL)) == NULL)
			die("could not allocate grid header", "");
		/* fill the empty grid with nan */
		for (node1 = 0; node1 < grid1->header->size; node1++)
			grid1->data[node1] = NAN;

		/* subswath 2 */
		wesn[GMT_XLO] = 0;
		wesn[GMT_YLO] = 0;
		wesn[GMT_XHI] = num_rng[1];
		wesn[GMT_YHI] = num_azi[1];
		inc[GMT_X] = (wesn[GMT_XHI] - wesn[GMT_XLO]) / num_x;
		inc[GMT_Y] = (wesn[GMT_YHI] - wesn[GMT_YLO]) / num_y;

		fprintf(stderr, "wesn = %f %f %f %f \n", wesn[GMT_XLO], wesn[GMT_YLO], wesn[GMT_XHI], wesn[GMT_YHI]);
		if ((grid2 = GMT_Create_Data(API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, wesn, inc, GMT_GRID_PIXEL_REG, 0,
		                             NULL)) == NULL)
			die("could not allocate grid header", "");
		/* fill the empty grid with nan */
		for (node2 = 0; node2 < grid2->header->size; node2++)
			grid2->data[node2] = NAN;

		/* subswath 3 */
		wesn[GMT_XLO] = 0;
		wesn[GMT_YLO] = 0;
		wesn[GMT_XHI] = num_rng[2];
		wesn[GMT_YHI] = num_azi[2];
		inc[GMT_X] = (wesn[GMT_XHI] - wesn[GMT_XLO]) / num_x;
		inc[GMT_Y] = (wesn[GMT_YHI] - wesn[GMT_YLO]) / num_y;

		fprintf(stderr, "wesn = %f %f %f %f \n", wesn[GMT_XLO], wesn[GMT_YLO], wesn[GMT_XHI], wesn[GMT_YHI]);
		if ((grid3 = GMT_Create_Data(API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, wesn, inc, GMT_GRID_PIXEL_REG, 0,
		                             NULL)) == NULL)
			die("could not allocate grid header", "");
		/* fill the empty grid with nan */
		for (node3 = 0; node3 < grid3->header->size; node3++)
			grid3->data[node3] = NAN;

		/* subswath 4 */
		wesn[GMT_XLO] = 0;
		wesn[GMT_YLO] = 0;
		wesn[GMT_XHI] = num_rng[3];
		wesn[GMT_YHI] = num_azi[3];
		inc[GMT_X] = (wesn[GMT_XHI] - wesn[GMT_XLO]) / num_x;
		inc[GMT_Y] = (wesn[GMT_YHI] - wesn[GMT_YLO]) / num_y;

		fprintf(stderr, "wesn = %f %f %f %f \n", wesn[GMT_XLO], wesn[GMT_YLO], wesn[GMT_XHI], wesn[GMT_YHI]);
		if ((grid4 = GMT_Create_Data(API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, wesn, inc, GMT_GRID_PIXEL_REG, 0,
		                             NULL)) == NULL)
			die("could not allocate grid header", "");
		/* fill the empty grid with nan */
		for (node4 = 0; node4 < grid4->header->size; node4++)
			grid4->data[node4] = NAN;

		/* subswath 5 */
		wesn[GMT_XLO] = 0;
		wesn[GMT_YLO] = 0;
		wesn[GMT_XHI] = num_rng[4];
		wesn[GMT_YHI] = num_azi[4];
		inc[GMT_X] = (wesn[GMT_XHI] - wesn[GMT_XLO]) / num_x;
		inc[GMT_Y] = (wesn[GMT_YHI] - wesn[GMT_YLO]) / num_y;

		fprintf(stderr, "wesn = %f %f %f %f \n", wesn[GMT_XLO], wesn[GMT_YLO], wesn[GMT_XHI], wesn[GMT_YHI]);
		if ((grid5 = GMT_Create_Data(API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, wesn, inc, GMT_GRID_PIXEL_REG, 0,
		                             NULL)) == NULL)
			die("could not allocate grid header", "");
		/* fill the empty grid with nan */
		for (node5 = 0; node5 < grid5->header->size; node5++)
			grid5->data[node5] = NAN;

		/* subswath 6 */
		wesn[GMT_XLO] = 0;
		wesn[GMT_YLO] = 0;
		wesn[GMT_XHI] = num_rng[5];
		wesn[GMT_YHI] = num_azi[5];
		inc[GMT_X] = (wesn[GMT_XHI] - wesn[GMT_XLO]) / num_x;
		inc[GMT_Y] = (wesn[GMT_YHI] - wesn[GMT_YLO]) / num_y;

		fprintf(stderr, "wesn = %f %f %f %f \n", wesn[GMT_XLO], wesn[GMT_YLO], wesn[GMT_XHI], wesn[GMT_YHI]);
		if ((grid6 = GMT_Create_Data(API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, wesn, inc, GMT_GRID_PIXEL_REG, 0,
		                             NULL)) == NULL)
			die("could not allocate grid header", "");
		/* fill the empty grid with nan */
		for (node6 = 0; node6 < grid6->header->size; node6++)
			grid6->data[node6] = NAN;

		/* subswath 7 */
		wesn[GMT_XLO] = 0;
		wesn[GMT_YLO] = 0;
		wesn[GMT_XHI] = num_rng[6];
		wesn[GMT_YHI] = num_azi[6];
		inc[GMT_X] = (wesn[GMT_XHI] - wesn[GMT_XLO]) / num_x;
		inc[GMT_Y] = (wesn[GMT_YHI] - wesn[GMT_YLO]) / num_y;

		fprintf(stderr, "wesn = %f %f %f %f \n", wesn[GMT_XLO], wesn[GMT_YLO], wesn[GMT_XHI], wesn[GMT_YHI]);
		if ((grid7 = GMT_Create_Data(API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, wesn, inc, GMT_GRID_PIXEL_REG, 0,
		                             NULL)) == NULL)
			die("could not allocate grid header", "");
		/* fill the empty grid with nan */
		for (node7 = 0; node7 < grid7->header->size; node7++)
			grid7->data[node7] = NAN;

		/* subswath 8 */
		wesn[GMT_XLO] = 0;
		wesn[GMT_YLO] = 0;
		wesn[GMT_XHI] = num_rng[7];
		wesn[GMT_YHI] = num_azi[7];
		inc[GMT_X] = (wesn[GMT_XHI] - wesn[GMT_XLO]) / num_x;
		inc[GMT_Y] = (wesn[GMT_YHI] - wesn[GMT_YLO]) / num_y;

		fprintf(stderr, "wesn = %f %f %f %f \n", wesn[GMT_XLO], wesn[GMT_YLO], wesn[GMT_XHI], wesn[GMT_YHI]);
		if ((grid8 = GMT_Create_Data(API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, wesn, inc, GMT_GRID_PIXEL_REG, 0,
		                             NULL)) == NULL)
			die("could not allocate grid header", "");
		/* fill the empty grid with nan */
		for (node8 = 0; node8 < grid8->header->size; node8++)
			grid8->data[node8] = NAN;

		/* subswath 9 */
		wesn[GMT_XLO] = 0;
		wesn[GMT_YLO] = 0;
		wesn[GMT_XHI] = num_rng[8];
		wesn[GMT_YHI] = num_azi[8];
		inc[GMT_X] = (wesn[GMT_XHI] - wesn[GMT_XLO]) / num_x;
		inc[GMT_Y] = (wesn[GMT_YHI] - wesn[GMT_YLO]) / num_y;

		fprintf(stderr, "wesn = %f %f %f %f \n", wesn[GMT_XLO], wesn[GMT_YLO], wesn[GMT_XHI], wesn[GMT_YHI]);
		if ((grid9 = GMT_Create_Data(API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, wesn, inc, GMT_GRID_PIXEL_REG, 0,
		                             NULL)) == NULL)
			die("could not allocate grid header", "");
		/* fill the empty grid with nan */
		for (node9 = 0; node9 < grid9->header->size; node9++)
			grid9->data[node9] = NAN;

		/* subswath 10 */
		wesn[GMT_XLO] = 0;
		wesn[GMT_YLO] = 0;
		wesn[GMT_XHI] = num_rng[9];
		wesn[GMT_YHI] = num_azi[9];
		inc[GMT_X] = (wesn[GMT_XHI] - wesn[GMT_XLO]) / num_x;
		inc[GMT_Y] = (wesn[GMT_YHI] - wesn[GMT_YLO]) / num_y;

		fprintf(stderr, "wesn = %f %f %f %f \n", wesn[GMT_XLO], wesn[GMT_YLO], wesn[GMT_XHI], wesn[GMT_YHI]);
		if ((grid10 = GMT_Create_Data(API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, wesn, inc, GMT_GRID_PIXEL_REG, 0,
		                              NULL)) == NULL)
			die("could not allocate grid header", "");
		/* fill the empty grid with nan */ for (node10 = 0; node10 < grid10->header->size; node10++)
			grid10->data[node10] = NAN;
	}

	else {
		die("dir must be either 1 or 0", "");
	}

	/* find the minimum of the clock start */
	clock_start_min = 1e10;
	for (i = 0; i < 5; i++) {
		if (clock_start[i] < clock_start_min) {
			clock_start_min = clock_start[i];
		}
	}

	/* check the range sampling rate should be the same */
	rsr = rng_samp_rate[0];
	for (i = 1; i < 10; i++) {
		if (rsr != rng_samp_rate[i]) {
			die("\nrange sampling rate from PRM files are not the same.\n", " ");
		}
	}
	dr = sol / rsr / 2.0;

	fprintf(stderr, "dr = %f\n", dr);

	for (i = 0; i < 10; i++) {
		fprintf(stderr, "subswath %d, da = %f\n", i + 1, vel[i] / prf[i]);
	}

	/* calculate the relative clock in seconds */
	for (i = 0; i < 10; i++) {
		clock_start_rel[i] = clock_start[i] - clock_start_min;
		near_range_rel[i] = near_range[i] - near_range[0];
		fprintf(stderr, "subswath = %d, clock_start_rel = %f\n", i + 1, clock_start_rel[i]);
		fprintf(stderr, "subswath = %d, near_range_rel = %f\n", i + 1, near_range_rel[i]);
	}

	/* calculate the positions of the grids in dec meters */
	ra[0][0] = round(near_range_rel[0] / dec);
	ra[0][1] = ra[0][0] + round(grid1->header->n_columns * grid1->header->inc[GMT_X] * dr / dec);
	ra[0][2] = round(clock_start_rel[0] * vel[0] / dec);
	ra[0][3] = ra[0][2] + round(grid1->header->n_rows * grid1->header->inc[GMT_Y] * vel[0] / prf[0] / dec);

	ra[1][0] = round(near_range_rel[1] / dec);
	ra[1][1] = ra[1][0] + round(grid2->header->n_columns * grid2->header->inc[GMT_X] * dr / dec);
	ra[1][2] = round(clock_start_rel[1] * vel[1] / dec);
	ra[1][3] = ra[1][2] + round(grid2->header->n_rows * grid2->header->inc[GMT_Y] * vel[1] / prf[1] / dec);

	ra[2][0] = round(near_range_rel[2] / dec);
	ra[2][1] = ra[2][0] + round(grid3->header->n_columns * grid3->header->inc[GMT_X] * dr / dec);
	ra[2][2] = round(clock_start_rel[2] * vel[2] / dec);
	ra[2][3] = ra[2][2] + round(grid3->header->n_rows * grid3->header->inc[GMT_Y] * vel[2] / prf[2] / dec);

	ra[3][0] = round(near_range_rel[3] / dec);
	ra[3][1] = ra[3][0] + round(grid4->header->n_columns * grid4->header->inc[GMT_X] * dr / dec);
	ra[3][2] = round(clock_start_rel[3] * vel[3] / dec);
	ra[3][3] = ra[3][2] + round(grid4->header->n_rows * grid4->header->inc[GMT_Y] * vel[3] / prf[3] / dec);

	ra[4][0] = round(near_range_rel[4] / dec);
	ra[4][1] = ra[4][0] + round(grid5->header->n_columns * grid5->header->inc[GMT_X] * dr / dec);
	ra[4][2] = round(clock_start_rel[4] * vel[4] / dec);
	ra[4][3] = ra[4][2] + round(grid5->header->n_rows * grid5->header->inc[GMT_Y] * vel[4] / prf[4] / dec);

	ra[5][0] = round(near_range_rel[5] / dec);
	ra[5][1] = ra[5][0] + round(grid6->header->n_columns * grid6->header->inc[GMT_X] * dr / dec);
	ra[5][2] = round(clock_start_rel[5] * vel[5] / dec);
	ra[5][3] = ra[5][2] + round(grid6->header->n_rows * grid6->header->inc[GMT_Y] * vel[5] / prf[5] / dec);

	ra[6][0] = round(near_range_rel[6] / dec);
	ra[6][1] = ra[6][0] + round(grid7->header->n_columns * grid7->header->inc[GMT_X] * dr / dec);
	ra[6][2] = round(clock_start_rel[6] * vel[6] / dec);
	ra[6][3] = ra[6][2] + round(grid7->header->n_rows * grid7->header->inc[GMT_Y] * vel[6] / prf[6] / dec);

	ra[7][0] = round(near_range_rel[7] / dec);
	ra[7][1] = ra[7][0] + round(grid8->header->n_columns * grid8->header->inc[GMT_X] * dr / dec);
	ra[7][2] = round(clock_start_rel[7] * vel[7] / dec);
	ra[7][3] = ra[7][2] + round(grid8->header->n_rows * grid8->header->inc[GMT_Y] * vel[7] / prf[7] / dec);

	ra[8][0] = round(near_range_rel[8] / dec);
	ra[8][1] = ra[8][0] + round(grid9->header->n_columns * grid9->header->inc[GMT_X] * dr / dec);
	ra[8][2] = round(clock_start_rel[8] * vel[8] / dec);
	ra[8][3] = ra[8][2] + round(grid9->header->n_rows * grid9->header->inc[GMT_Y] * vel[8] / prf[8] / dec);

	ra[9][0] = round(near_range_rel[9] / dec);
	ra[9][1] = ra[9][0] + round(grid10->header->n_columns * grid10->header->inc[GMT_X] * dr / dec);
	ra[9][2] = round(clock_start_rel[9] * vel[9] / dec);
	ra[9][3] = ra[9][2] + round(grid10->header->n_rows * grid10->header->inc[GMT_Y] * vel[9] / prf[9] / dec);

	/* debug */
	/*for (i=0;i<5;i++) {
	  for (j=0;j<4;j++) {
	    fprintf(stderr, "i=%d j=%d ra=%d\n", i,j,ra[i][j]);
	  }
	}*/

	if (dir == 1) {
		fprintf(stderr, "mosaic the 5 subswath...\n");

		/* the lower boundary in the azimuth and range direction is 0 */
		/* find the upper boundary in the azimuth and range direction */
		azimuth_ub = -999999;
		for (i = 5; i < 10; i++) {
			if (ra[i][3] > azimuth_ub) {
				azimuth_ub = ra[i][3];
			}
		}
		range_ub = ra[4][1];

		/* make an output grid */
		inc[GMT_X] = inc[GMT_Y] = dec; /* in meters */
		wesn[GMT_XLO] = 0;
		wesn[GMT_YLO] = 0;
		wesn[GMT_XHI] = range_ub * inc[GMT_X];
		wesn[GMT_YHI] = azimuth_ub * inc[GMT_Y];

		/* debug */
		/* fprintf(stderr, "inc[GMT_X] = %f , inc[GMT_Y] = %f\n", inc[GMT_X],
		inc[GMT_Y]); fprintf(stderr, "wesn = %f %f %f %f \n", wesn[GMT_XLO],
		wesn[GMT_YLO], wesn[GMT_XHI], wesn[GMT_YHI]); */

		if ((grid11 = GMT_Create_Data(API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, wesn, inc, GMT_GRID_PIXEL_REG, 0,
		                              NULL)) == NULL)
			die("could not allocate grid header", "");

		/* fill the empty grid with nan */
		for (node11 = 0; node11 < grid11->header->size; node11++)
			grid11->data[node11] = NAN;

		/* put in the subswath grids */
		range_boundary1 = round((ra[0][1] + ra[1][0]) / 2);
		azimuth_boundary1 = round((ra[0][3] + ra[5][2]) / 2);
		/* debug */
		// fprintf(stderr,"range_boundary1=%d\n",range_boundary1);

		/* put in subswath 1 */
		for (range = ra[0][0]; range < range_boundary1; range++) {
			for (azimuth = ra[0][2]; azimuth < azimuth_boundary1; azimuth++) {
				node11 = GMT_Get_Index(API, grid11->header, azimuth, range);
				col = round((range * dec - near_range_rel[0]) / dr / grid1->header->inc[GMT_X]);
				row = grid1->header->n_rows -
				      round((azimuth * dec - clock_start_rel[0] * vel[0]) / (vel[0] / prf[0]) / grid1->header->inc[GMT_Y]);
				if (col >= 0 && col <= grid1->header->n_columns && row >= 0 && row <= grid1->header->n_rows) {
					node1 = GMT_Get_Index(API, grid1->header, row, col);
					grid11->data[node11] = grid1->data[node1];
				}
			}
		}

		/* put in subswath 2 */
		range_boundary2 = round((ra[1][1] + ra[2][0]) / 2);
		azimuth_boundary2 = round((ra[1][3] + ra[6][2]) / 2);
		/* debug */
		// fprintf(stderr,"range_boundary2=%d\n",range_boundary2);
		for (range = range_boundary1; range < range_boundary2; range++) {
			for (azimuth = ra[1][2]; azimuth < azimuth_boundary2; azimuth++) {
				node11 = GMT_Get_Index(API, grid11->header, azimuth, range);
				col = round((range * dec - near_range_rel[1]) / dr / grid2->header->inc[GMT_X]);
				row = grid2->header->n_rows -
				      round((azimuth * dec - clock_start_rel[1] * vel[1]) / (vel[1] / prf[1]) / grid2->header->inc[GMT_Y]);
				if (col >= 0 && col <= grid2->header->n_columns && row >= 0 && row <= grid2->header->n_rows) {
					node2 = GMT_Get_Index(API, grid2->header, row, col);
					grid11->data[node11] = grid2->data[node2];
				}
			}
		}

		/* put in subswath 3 */
		range_boundary3 = round((ra[2][1] + ra[3][0]) / 2);
		azimuth_boundary3 = round((ra[2][3] + ra[7][2]) / 2);
		/* debug */
		// fprintf(stderr,"range_boundary3=%d\n",range_boundary3);
		for (range = range_boundary2; range < range_boundary3; range++) {
			for (azimuth = ra[2][2]; azimuth < azimuth_boundary3; azimuth++) {
				node11 = GMT_Get_Index(API, grid11->header, azimuth, range);
				col = round((range * dec - near_range_rel[2]) / dr / grid3->header->inc[GMT_X]);
				row = grid3->header->n_rows -
				      round((azimuth * dec - clock_start_rel[2] * vel[2]) / (vel[2] / prf[2]) / grid3->header->inc[GMT_Y]);
				if (col >= 0 && col <= grid3->header->n_columns && row >= 0 && row <= grid3->header->n_rows) {
					node3 = GMT_Get_Index(API, grid3->header, row, col);
					grid11->data[node11] = grid3->data[node3];
				}
			}
		}

		/* put in subswath 4 */
		range_boundary4 = round((ra[3][1] + ra[4][0]) / 2);
		azimuth_boundary4 = round((ra[3][3] + ra[8][2]) / 2);
		/* debug */
		// fprintf(stderr,"range_boundary4=%d\n",range_boundary4);
		for (range = range_boundary3; range < range_boundary4; range++) {
			for (azimuth = ra[3][2]; azimuth < azimuth_boundary4; azimuth++) {
				node11 = GMT_Get_Index(API, grid11->header, azimuth, range);
				col = round((range * dec - near_range_rel[3]) / dr / grid4->header->inc[GMT_X]);
				row = grid4->header->n_rows -
				      round((azimuth * dec - clock_start_rel[3] * vel[3]) / (vel[3] / prf[3]) / grid4->header->inc[GMT_Y]);
				if (col >= 0 && col <= grid4->header->n_columns && row >= 0 && row <= grid4->header->n_rows) {
					node4 = GMT_Get_Index(API, grid4->header, row, col);
					grid11->data[node11] = grid4->data[node4];
				}
			}
		}

		/* put in subswath 5 */
		range_boundary5 = range_ub;
		azimuth_boundary5 = round((ra[4][3] + ra[9][2]) / 2);
		/* debug */
		// fprintf(stderr,"range_boundary5=%d\n",range_boundary5);
		for (range = range_boundary4; range < range_boundary5; range++) {
			for (azimuth = ra[4][2]; azimuth < azimuth_boundary5; azimuth++) {
				node11 = GMT_Get_Index(API, grid11->header, azimuth, range);
				col = round((range * dec - near_range_rel[4]) / dr / grid5->header->inc[GMT_X]);
				row = grid5->header->n_rows -
				      round((azimuth * dec - clock_start_rel[4] * vel[4]) / (vel[4] / prf[4]) / grid5->header->inc[GMT_Y]);
				if (col >= 0 && col <= grid5->header->n_columns && row >= 0 && row <= grid5->header->n_rows) {
					node5 = GMT_Get_Index(API, grid5->header, row, col);
					grid11->data[node11] = grid5->data[node5];
				}
			}
		}

		/* put in subswath 6 */
		range_boundary1 = round((ra[5][1] + ra[6][0]) / 2);
		for (range = ra[5][0]; range < range_boundary1; range++) {
			for (azimuth = azimuth_boundary1; azimuth < ra[5][3]; azimuth++) {
				node11 = GMT_Get_Index(API, grid11->header, azimuth, range);
				col = round((range * dec - near_range_rel[5]) / dr / grid6->header->inc[GMT_X]);
				row = grid6->header->n_rows -
				      round((azimuth * dec - clock_start_rel[5] * vel[5]) / (vel[5] / prf[5]) / grid6->header->inc[GMT_Y]);
				if (col >= 0 && col <= grid6->header->n_columns && row >= 0 && row <= grid6->header->n_rows) {
					node6 = GMT_Get_Index(API, grid6->header, row, col);
					grid11->data[node11] = grid6->data[node6];
				}
			}
		}

		/* put in subswath 7 */
		range_boundary2 = round((ra[6][1] + ra[7][0]) / 2);
		for (range = range_boundary1; range < range_boundary2; range++) {
			for (azimuth = azimuth_boundary2; azimuth < ra[6][3]; azimuth++) {
				node11 = GMT_Get_Index(API, grid11->header, azimuth, range);
				col = round((range * dec - near_range_rel[6]) / dr / grid7->header->inc[GMT_X]);
				row = grid7->header->n_rows -
				      round((azimuth * dec - clock_start_rel[6] * vel[6]) / (vel[6] / prf[6]) / grid7->header->inc[GMT_Y]);
				if (col >= 0 && col <= grid7->header->n_columns && row >= 0 && row <= grid7->header->n_rows) {
					node7 = GMT_Get_Index(API, grid7->header, row, col);
					grid11->data[node11] = grid7->data[node7];
				}
			}
		}

		/* put in subswath 8 */
		range_boundary3 = round((ra[7][1] + ra[8][0]) / 2);
		for (range = range_boundary2; range < range_boundary3; range++) {
			for (azimuth = azimuth_boundary3; azimuth < ra[7][3]; azimuth++) {
				node11 = GMT_Get_Index(API, grid11->header, azimuth, range);
				col = round((range * dec - near_range_rel[7]) / dr / grid8->header->inc[GMT_X]);
				row = grid8->header->n_rows -
				      round((azimuth * dec - clock_start_rel[7] * vel[7]) / (vel[7] / prf[7]) / grid8->header->inc[GMT_Y]);
				if (col >= 0 && col <= grid8->header->n_columns && row >= 0 && row <= grid8->header->n_rows) {
					node8 = GMT_Get_Index(API, grid8->header, row, col);
					grid11->data[node11] = grid8->data[node8];
				}
			}
		}

		/* put in subswath 9 */
		range_boundary4 = round((ra[8][1] + ra[9][0]) / 2);
		for (range = range_boundary3; range < range_boundary4; range++) {
			for (azimuth = azimuth_boundary4; azimuth < ra[8][3]; azimuth++) {
				node11 = GMT_Get_Index(API, grid11->header, azimuth, range);
				col = round((range * dec - near_range_rel[8]) / dr / grid9->header->inc[GMT_X]);
				row = grid9->header->n_rows -
				      round((azimuth * dec - clock_start_rel[8] * vel[8]) / (vel[8] / prf[8]) / grid9->header->inc[GMT_Y]);
				if (col >= 0 && col <= grid9->header->n_columns && row >= 0 && row <= grid9->header->n_rows) {
					node9 = GMT_Get_Index(API, grid9->header, row, col);
					grid11->data[node11] = grid9->data[node9];
				}
			}
		}

		/* put in subswath 10 */
		range_boundary5 = ra[9][1];
		for (range = range_boundary4; range < range_boundary5; range++) {
			for (azimuth = azimuth_boundary5; azimuth < ra[9][3]; azimuth++) {
				node11 = GMT_Get_Index(API, grid11->header, azimuth, range);
				col = round((range * dec - near_range_rel[9]) / dr / grid10->header->inc[GMT_X]);
				row = grid10->header->n_rows -
				      round((azimuth * dec - clock_start_rel[9] * vel[9]) / (vel[9] / prf[9]) / grid10->header->inc[GMT_Y]);
				if (col >= 0 && col <= grid10->header->n_columns && row >= 0 && row <= grid10->header->n_rows) {
					node10 = GMT_Get_Index(API, grid10->header, row, col);
					grid11->data[node11] = grid10->data[node10];
				}
			}
		}

		/* output and clean up */
		if (GMT_Write_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, argv[21], grid11))
			return EXIT_FAILURE;
	}
	else if (dir == 0) {
		fprintf(stderr, "divide the mosaic into 5 subswath...\n");
		/* open the mosaic grid */
		if ((grid11 = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, argv[21], NULL)) ==
		    NULL)
			return EXIT_FAILURE;
		if (GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, argv[21], grid11) == NULL)
			return EXIT_FAILURE;

		/* subswath 1 */
		for (col = 0; col < grid1->header->n_columns; col++) {
			for (row = 0; row < grid1->header->n_rows; row++) {
				node1 = GMT_Get_Index(API, grid1->header, row, col);
				range = round((near_range_rel[0] + col * grid1->header->inc[GMT_X] * dr) / dec);
				azimuth = round(
				    (clock_start_rel[0] * vel[0] + (grid1->header->n_rows - row) * grid1->header->inc[GMT_Y] * vel[0] / prf[0]) /
				    dec);
				if (range >= 0 && range <= grid11->header->n_columns && azimuth >= 0 && azimuth <= grid11->header->n_rows) {
					node11 = GMT_Get_Index(API, grid11->header, azimuth, range);
					grid1->data[node1] = grid11->data[node11];
				}
			}
		}

		/* subswath 2 */
		for (col = 0; col < grid2->header->n_columns; col++) {
			for (row = 0; row < grid2->header->n_rows; row++) {
				node2 = GMT_Get_Index(API, grid2->header, row, col);
				range = round((near_range_rel[1] + col * grid2->header->inc[GMT_X] * dr) / dec);
				azimuth = round(
				    (clock_start_rel[1] * vel[1] + (grid2->header->n_rows - row) * grid2->header->inc[GMT_Y] * vel[1] / prf[1]) /
				    dec);
				if (range >= 0 && range <= grid11->header->n_columns && azimuth >= 0 && azimuth <= grid11->header->n_rows) {
					node11 = GMT_Get_Index(API, grid11->header, azimuth, range);
					grid2->data[node2] = grid11->data[node11];
				}
			}
		}

		/* subswath 3 */
		for (col = 0; col < grid3->header->n_columns; col++) {
			for (row = 0; row < grid3->header->n_rows; row++) {
				node3 = GMT_Get_Index(API, grid3->header, row, col);
				range = round((near_range_rel[2] + col * grid3->header->inc[GMT_X] * dr) / dec);
				azimuth = round(
				    (clock_start_rel[2] * vel[2] + (grid3->header->n_rows - row) * grid3->header->inc[GMT_Y] * vel[2] / prf[2]) /
				    dec);
				if (range >= 0 && range <= grid11->header->n_columns && azimuth >= 0 && azimuth <= grid11->header->n_rows) {
					node11 = GMT_Get_Index(API, grid11->header, azimuth, range);
					grid3->data[node3] = grid11->data[node11];
				}
			}
		}

		/* subswath 4 */
		for (col = 0; col < grid4->header->n_columns; col++) {
			for (row = 0; row < grid4->header->n_rows; row++) {
				node4 = GMT_Get_Index(API, grid4->header, row, col);
				range = round((near_range_rel[3] + col * grid4->header->inc[GMT_X] * dr) / dec);
				azimuth = round(
				    (clock_start_rel[3] * vel[3] + (grid4->header->n_rows - row) * grid4->header->inc[GMT_Y] * vel[3] / prf[3]) /
				    dec);
				if (range >= 0 && range <= grid11->header->n_columns && azimuth >= 0 && azimuth <= grid11->header->n_rows) {
					node11 = GMT_Get_Index(API, grid11->header, azimuth, range);
					grid4->data[node4] = grid11->data[node11];
				}
			}
		}

		/* subswath 5 */
		for (col = 0; col < grid5->header->n_columns; col++) {
			for (row = 0; row < grid5->header->n_rows; row++) {
				node5 = GMT_Get_Index(API, grid5->header, row, col);
				range = round((near_range_rel[4] + col * grid5->header->inc[GMT_X] * dr) / dec);
				azimuth = round(
				    (clock_start_rel[4] * vel[4] + (grid5->header->n_rows - row) * grid5->header->inc[GMT_Y] * vel[4] / prf[4]) /
				    dec);
				if (range >= 0 && range <= grid11->header->n_columns && azimuth >= 0 && azimuth <= grid11->header->n_rows) {
					node11 = GMT_Get_Index(API, grid11->header, azimuth, range);
					grid5->data[node5] = grid11->data[node11];
				}
			}
		}

		/* subswath 6 */
		for (col = 0; col < grid6->header->n_columns; col++) {
			for (row = 0; row < grid6->header->n_rows; row++) {
				node6 = GMT_Get_Index(API, grid6->header, row, col);
				range = round((near_range_rel[5] + col * grid6->header->inc[GMT_X] * dr) / dec);
				azimuth = round(
				    (clock_start_rel[5] * vel[5] + (grid6->header->n_rows - row) * grid6->header->inc[GMT_Y] * vel[5] / prf[5]) /
				    dec);
				if (range >= 0 && range <= grid11->header->n_columns && azimuth >= 0 && azimuth <= grid11->header->n_rows) {
					node11 = GMT_Get_Index(API, grid11->header, azimuth, range);
					grid6->data[node6] = grid11->data[node11];
				}
			}
		}

		/* subswath 7 */
		for (col = 0; col < grid7->header->n_columns; col++) {
			for (row = 0; row < grid7->header->n_rows; row++) {
				node7 = GMT_Get_Index(API, grid7->header, row, col);
				range = round((near_range_rel[6] + col * grid7->header->inc[GMT_X] * dr) / dec);
				azimuth = round(
				    (clock_start_rel[6] * vel[6] + (grid7->header->n_rows - row) * grid7->header->inc[GMT_Y] * vel[6] / prf[6]) /
				    dec);
				if (range >= 0 && range <= grid11->header->n_columns && azimuth >= 0 && azimuth <= grid11->header->n_rows) {
					node11 = GMT_Get_Index(API, grid11->header, azimuth, range);
					grid7->data[node7] = grid11->data[node11];
				}
			}
		}

		/* subswath 8 */
		for (col = 0; col < grid8->header->n_columns; col++) {
			for (row = 0; row < grid8->header->n_rows; row++) {
				node8 = GMT_Get_Index(API, grid8->header, row, col);
				range = round((near_range_rel[7] + col * grid8->header->inc[GMT_X] * dr) / dec);
				azimuth = round(
				    (clock_start_rel[7] * vel[7] + (grid8->header->n_rows - row) * grid8->header->inc[GMT_Y] * vel[7] / prf[7]) /
				    dec);
				if (range >= 0 && range <= grid11->header->n_columns && azimuth >= 0 && azimuth <= grid11->header->n_rows) {
					node11 = GMT_Get_Index(API, grid11->header, azimuth, range);
					grid8->data[node8] = grid11->data[node11];
				}
			}
		}

		/* subswath 9 */
		for (col = 0; col < grid9->header->n_columns; col++) {
			for (row = 0; row < grid9->header->n_rows; row++) {
				node9 = GMT_Get_Index(API, grid9->header, row, col);
				range = round((near_range_rel[8] + col * grid9->header->inc[GMT_X] * dr) / dec);
				azimuth = round(
				    (clock_start_rel[8] * vel[8] + (grid9->header->n_rows - row) * grid9->header->inc[GMT_Y] * vel[8] / prf[8]) /
				    dec);
				if (range >= 0 && range <= grid11->header->n_columns && azimuth >= 0 && azimuth <= grid11->header->n_rows) {
					node11 = GMT_Get_Index(API, grid11->header, azimuth, range);
					grid9->data[node9] = grid11->data[node11];
				}
			}
		}

		/* subswath 10 */
		for (col = 0; col < grid10->header->n_columns; col++) {
			for (row = 0; row < grid10->header->n_rows; row++) {
				node10 = GMT_Get_Index(API, grid10->header, row, col);
				range = round((near_range_rel[9] + col * grid10->header->inc[GMT_X] * dr) / dec);
				azimuth = round(
				    (clock_start_rel[9] * vel[9] + (grid10->header->n_rows - row) * grid10->header->inc[GMT_Y] * vel[9] / prf[9]) /
				    dec);
				if (range >= 0 && range <= grid11->header->n_columns && azimuth >= 0 && azimuth <= grid11->header->n_rows) {
					node11 = GMT_Get_Index(API, grid11->header, azimuth, range);
					grid10->data[node10] = grid11->data[node11];
				}
			}
		}

		if (GMT_Write_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, argv[2], grid1))
			return EXIT_FAILURE;

		if (GMT_Write_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, argv[4], grid2))
			return EXIT_FAILURE;

		if (GMT_Write_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, argv[6], grid3))
			return EXIT_FAILURE;

		if (GMT_Write_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, argv[8], grid4))
			return EXIT_FAILURE;

		if (GMT_Write_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, argv[10], grid5))
			return EXIT_FAILURE;

		if (GMT_Write_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, argv[12], grid6))
			return EXIT_FAILURE;

		if (GMT_Write_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, argv[14], grid7))
			return EXIT_FAILURE;

		if (GMT_Write_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, argv[16], grid8))
			return EXIT_FAILURE;

		if (GMT_Write_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, argv[18], grid9))
			return EXIT_FAILURE;

		if (GMT_Write_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, argv[20], grid10))
			return EXIT_FAILURE;
	}

	if (GMT_Destroy_Session(API))
		return EXIT_FAILURE; /* Remove the GMT machinery */

	return (EXIT_SUCCESS);
}
