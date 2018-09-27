/**************************************************************************/
/*                                                                        */
/*   Mosaic 5 subswaths of ALOS ScanSAR in radar coordinates              */
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
/*     May 30 2018: divide mosaic into 5 subswath                         */
/*                                                                        */
/**************************************************************************/

char *USAGE = " ALOS_mosaic_ss prmfile1 grd1 prmfile2 grd2 prmfile3 grd3 prmfile4 grd4 "
              "prmfile5 grd5 mosaic dec dir \n\n"
              " USAGE:  \n"
              "  mosaic 5 subswath of ALOS ScanSAR in radar coordinates\n\n"
              "  prmfile[1-5]    the PRM files \n"
              "  grd[1-5]        the GMT grid files (corr, phase, phasefilt) \n"
              "  dec             decimation factor in radar coordinates (100) \n"
              "  dir             dir is 1 means mosaic 5 subswaths \n"
              "                  dir is 0 means divide the mosaic into 5 subswaths \n";

#include "gmt.h"
#include "gmtsar.h"

void set_ALOS_defaults(struct PRM *);
void get_sio_struct(FILE *, struct PRM *);

int main(int argc, char **argv) {

	/* define variables */
	FILE *prmfile;
	struct PRM prm;
	void *API = NULL; /* GMT control structure */
	struct GMT_GRID *grid1 = NULL, *grid2 = NULL, *grid3 = NULL, *grid4 = NULL, *grid5 = NULL, *grid6 = NULL;
	double inc[2], wesn[4];
	double clock_start[5], clock_start_rel[5], near_range[5], near_range_rel[5];
	double prf[5], rng_samp_rate[5], vel[5];
	int num_rng[5], num_azi[5];
	double clock_start_min;
	int ra[5][4];
	int range_ub, azimuth_ub;
	double sol = 3.0e8;
	double rsr, dr;
	double dec;
	int row, col;
	int i;
	int range, azimuth;
	int range_boundary1, range_boundary2, range_boundary3, range_boundary4, range_boundary5;
	int node1, node2, node3, node4, node5, node6;
	int dir;
	double num_x, num_y;

	/* parse arguments */
	if (argc != 14)
		die("\n", USAGE);

	/* Begin: Initializing new GMT session */
	if ((API = GMT_Create_Session(argv[0], 2, 0, NULL)) == NULL)
		return EXIT_FAILURE;

	/* read in dec */
	dec = atof(argv[12]);

	/* read in dir */
	dir = atoi(argv[13]);

	num_x = 1000;
	num_y = 5000;

	/* open and read PRM files */
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
		        grid1->header->nx, grid1->header->inc[GMT_X], grid1->header->ny, grid1->header->inc[GMT_Y]);

		if ((num_rng[0] != round(grid1->header->nx * grid1->header->inc[GMT_X])) ||
		    (num_azi[0] != round(grid1->header->ny * grid1->header->inc[GMT_Y]))) {
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
		        grid2->header->nx, grid2->header->inc[GMT_X], grid2->header->ny, grid2->header->inc[GMT_Y]);
		if ((num_rng[1] != round(grid2->header->nx * grid2->header->inc[GMT_X])) ||
		    (num_azi[1] != round(grid2->header->ny * grid2->header->inc[GMT_Y]))) {
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
		        grid3->header->nx, grid3->header->inc[GMT_X], grid3->header->ny, grid3->header->inc[GMT_Y]);
		if ((num_rng[2] != round(grid3->header->nx * grid3->header->inc[GMT_X])) ||
		    (num_azi[2] != round(grid3->header->ny * grid3->header->inc[GMT_Y]))) {
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
		        grid4->header->nx, grid4->header->inc[GMT_X], grid4->header->ny, grid4->header->inc[GMT_Y]);
		if ((num_rng[3] != round(grid4->header->nx * grid4->header->inc[GMT_X])) ||
		    (num_azi[3] != round(grid4->header->ny * grid4->header->inc[GMT_Y]))) {
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
		        grid5->header->nx, grid5->header->inc[GMT_X], grid5->header->ny, grid5->header->inc[GMT_Y]);
		if ((num_rng[4] != round(grid5->header->nx * grid5->header->inc[GMT_X])) ||
		    (num_azi[4] != round(grid5->header->ny * grid5->header->inc[GMT_Y]))) {
			die("\nPRM file and GMT grid don't match!\n", argv[10]);
		}
		if (GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, argv[10], grid5) == NULL)
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
		if (GMT_Set_Comment(API, GMT_IS_GRID, GMT_COMMENT_IS_REMARK, "subswath 1", grid1))
			return EXIT_FAILURE;
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
		if (GMT_Set_Comment(API, GMT_IS_GRID, GMT_COMMENT_IS_REMARK, "subswath 2", grid2))
			return EXIT_FAILURE;
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
		if (GMT_Set_Comment(API, GMT_IS_GRID, GMT_COMMENT_IS_REMARK, "subswath 3", grid3))
			return EXIT_FAILURE;
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
		if (GMT_Set_Comment(API, GMT_IS_GRID, GMT_COMMENT_IS_REMARK, "subswath 4", grid4))
			return EXIT_FAILURE;
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
		if (GMT_Set_Comment(API, GMT_IS_GRID, GMT_COMMENT_IS_REMARK, "subswath 5", grid5))
			return EXIT_FAILURE;
		/* fill the empty grid with nan */
		for (node5 = 0; node5 < grid5->header->size; node5++)
			grid5->data[node5] = NAN;
	}
	else {
		die("dir must be either 1 or 0", "");
	}

	/* find the maximum of the clock start */
	clock_start_min = 1e10;
	for (i = 0; i < 5; i++) {
		if (clock_start[i] < clock_start_min) {
			clock_start_min = clock_start[i];
		}
	}

	/* check the range sampling rate should be the same */
	rsr = rng_samp_rate[0];
	for (i = 1; i < 5; i++) {
		if (rsr != rng_samp_rate[i]) {
			die("\nrange sampling rate from PRM files are not the same.\n", " ");
		}
	}
	dr = sol / rsr / 2.0;

	fprintf(stderr, "dr = %f\n", dr);

	for (i = 0; i < 5; i++) {
		fprintf(stderr, "subswath %d, da = %f\n", i + 1, vel[i] / prf[i]);
	}

	/* calculate the relative clock in seconds */
	for (i = 0; i < 5; i++) {
		clock_start_rel[i] = clock_start[i] - clock_start_min;
		near_range_rel[i] = near_range[i] - near_range[0];
		fprintf(stderr, "subswath = %d, clock_start_rel = %f\n", i + 1, clock_start_rel[i]);
		fprintf(stderr, "subswath = %d, near_range_rel = %f\n", i + 1, near_range_rel[i]);
	}

	/* calculate the positions of the grids in dec meters */
	ra[0][0] = round(near_range_rel[0] / dec);
	ra[0][1] = ra[0][0] + round(grid1->header->nx * grid1->header->inc[GMT_X] * dr / dec);
	ra[0][2] = round(clock_start_rel[0] * vel[0] / dec);
	ra[0][3] = ra[0][2] + round(grid1->header->ny * grid1->header->inc[GMT_Y] * vel[0] / prf[0] / dec);

	ra[1][0] = round(near_range_rel[1] / dec);
	ra[1][1] = ra[1][0] + round(grid2->header->nx * grid2->header->inc[GMT_X] * dr / dec);
	ra[1][2] = round(clock_start_rel[1] * vel[1] / dec);
	ra[1][3] = ra[1][2] + round(grid2->header->ny * grid2->header->inc[GMT_Y] * vel[1] / prf[1] / dec);

	ra[2][0] = round(near_range_rel[2] / dec);
	ra[2][1] = ra[2][0] + round(grid3->header->nx * grid3->header->inc[GMT_X] * dr / dec);
	ra[2][2] = round(clock_start_rel[2] * vel[2] / dec);
	ra[2][3] = ra[2][2] + round(grid3->header->ny * grid3->header->inc[GMT_Y] * vel[2] / prf[2] / dec);

	ra[3][0] = round(near_range_rel[3] / dec);
	ra[3][1] = ra[3][0] + round(grid4->header->nx * grid4->header->inc[GMT_X] * dr / dec);
	ra[3][2] = round(clock_start_rel[3] * vel[3] / dec);
	ra[3][3] = ra[3][2] + round(grid4->header->ny * grid4->header->inc[GMT_Y] * vel[3] / prf[3] / dec);

	ra[4][0] = round(near_range_rel[4] / dec);
	ra[4][1] = ra[4][0] + round(grid5->header->nx * grid5->header->inc[GMT_X] * dr / dec);
	ra[4][2] = round(clock_start_rel[4] * vel[4] / dec);
	ra[4][3] = ra[4][2] + round(grid5->header->ny * grid5->header->inc[GMT_Y] * vel[4] / prf[4] / dec);

	/* debug */
	/*for (i=0;i<5;i++) {
	  for (j=0;j<4;j++) {
	    fprintf(stderr, "i=%d j=%d ra=%d\n", i,j,ra[i][j]);
	  }
	}*/

	range_boundary1 = round((ra[0][1] + ra[1][0]) / 2);
	range_boundary2 = round((ra[1][1] + ra[2][0]) / 2);
	range_boundary3 = round((ra[2][1] + ra[3][0]) / 2);
	range_boundary4 = round((ra[3][1] + ra[4][0]) / 2);
	range_boundary5 = ra[4][1];

	if (dir == 1) {
		fprintf(stderr, "mosaic the 5 subswath...\n");
		/* the lower boundary in the azimuth and range direction is 0 */
		/* find the upper boundary in the azimuth and range direction */
		azimuth_ub = -999999;
		for (i = 0; i < 5; i++) {
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

		if ((grid6 = GMT_Create_Data(API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, wesn, inc, GMT_GRID_PIXEL_REG, 0,
		                             NULL)) == NULL)
			die("could not allocate grid header", "");

		if (GMT_Set_Comment(API, GMT_IS_GRID, GMT_COMMENT_IS_REMARK, "mosaic", grid6))
			return EXIT_FAILURE;

		/* fill the empty grid with nan */
		for (node6 = 0; node6 < grid6->header->size; node6++)
			grid6->data[node6] = NAN;

		/* put in the subswath grids */
		/* debug */
		// fprintf(stderr,"range_boundary1=%d\n",range_boundary1);

		/* put in subswath 1 */
		for (range = ra[0][0]; range < range_boundary1; range++) {
			for (azimuth = ra[0][2]; azimuth < ra[0][3]; azimuth++) {
				node6 = GMT_Get_Index(API, grid6->header, azimuth, range);
				col = round((range * dec - near_range_rel[0]) / dr / grid1->header->inc[GMT_X]);
				row = grid1->header->ny -
				      round((azimuth * dec - clock_start_rel[0] * vel[0]) / (vel[0] / prf[0]) / grid1->header->inc[GMT_Y]);
				if (col >= 0 && col <= grid1->header->nx && row >= 0 && row <= grid1->header->ny) {
					node1 = GMT_Get_Index(API, grid1->header, row, col);
					grid6->data[node6] = grid1->data[node1];
				}
			}
		}

		/* put in subswath 2 */
		/* debug */
		// fprintf(stderr,"range_boundary2=%d\n",range_boundary2);
		for (range = range_boundary1; range < range_boundary2; range++) {
			for (azimuth = ra[1][2]; azimuth < ra[1][3]; azimuth++) {
				node6 = GMT_Get_Index(API, grid6->header, azimuth, range);
				col = round((range * dec - near_range_rel[1]) / dr / grid2->header->inc[GMT_X]);
				row = grid2->header->ny -
				      round((azimuth * dec - clock_start_rel[1] * vel[1]) / (vel[1] / prf[1]) / grid2->header->inc[GMT_Y]);
				if (col >= 0 && col <= grid2->header->nx && row >= 0 && row <= grid2->header->ny) {
					node2 = GMT_Get_Index(API, grid2->header, row, col);
					grid6->data[node6] = grid2->data[node2];
				}
			}
		}

		/* put in subswath 3 */
		/* debug */
		// fprintf(stderr,"range_boundary3=%d\n",range_boundary3);
		for (range = range_boundary2; range < range_boundary3; range++) {
			for (azimuth = ra[2][2]; azimuth < ra[2][3]; azimuth++) {
				node6 = GMT_Get_Index(API, grid6->header, azimuth, range);
				col = round((range * dec - near_range_rel[2]) / dr / grid3->header->inc[GMT_X]);
				row = grid3->header->ny -
				      round((azimuth * dec - clock_start_rel[2] * vel[2]) / (vel[2] / prf[2]) / grid3->header->inc[GMT_Y]);
				if (col >= 0 && col <= grid3->header->nx && row >= 0 && row <= grid3->header->ny) {
					node3 = GMT_Get_Index(API, grid3->header, row, col);
					grid6->data[node6] = grid3->data[node3];
				}
			}
		}

		/* put in subswath 4 */
		/* debug */
		// fprintf(stderr,"range_boundary4=%d\n",range_boundary4);
		for (range = range_boundary3; range < range_boundary4; range++) {
			for (azimuth = ra[3][2]; azimuth < ra[3][3]; azimuth++) {
				node6 = GMT_Get_Index(API, grid6->header, azimuth, range);
				col = round((range * dec - near_range_rel[3]) / dr / grid4->header->inc[GMT_X]);
				row = grid4->header->ny -
				      round((azimuth * dec - clock_start_rel[3] * vel[3]) / (vel[3] / prf[3]) / grid4->header->inc[GMT_Y]);
				if (col >= 0 && col <= grid4->header->nx && row >= 0 && row <= grid4->header->ny) {
					node4 = GMT_Get_Index(API, grid4->header, row, col);
					grid6->data[node6] = grid4->data[node4];
				}
			}
		}

		/* put in subswath 5 */
		/* debug */
		// fprintf(stderr,"range_boundary5=%d\n",range_boundary5);
		for (range = range_boundary4; range < range_boundary5; range++) {
			for (azimuth = ra[4][2]; azimuth < ra[4][3]; azimuth++) {
				node6 = GMT_Get_Index(API, grid6->header, azimuth, range);
				col = round((range * dec - near_range_rel[4]) / dr / grid5->header->inc[GMT_X]);
				row = grid5->header->ny -
				      round((azimuth * dec - clock_start_rel[4] * vel[4]) / (vel[4] / prf[4]) / grid5->header->inc[GMT_Y]);
				if (col >= 0 && col <= grid5->header->nx && row >= 0 && row <= grid5->header->ny) {
					node5 = GMT_Get_Index(API, grid5->header, row, col);
					grid6->data[node6] = grid5->data[node5];
				}
			}
		}

		/* output and clean up */
		if (GMT_Write_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, argv[11], grid6))
			return EXIT_FAILURE;
	}

	else if (dir == 0) {
		fprintf(stderr, "divide the mosaic into 5 subswath...\n");
		/* open the mosaic grid */
		if ((grid6 = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, argv[11], NULL)) ==
		    NULL)
			return EXIT_FAILURE;
		if (GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, argv[11], grid6) == NULL)
			return EXIT_FAILURE;

		/* subswath 1 */
		for (col = 0; col < grid1->header->nx; col++) {
			for (row = 0; row < grid1->header->ny; row++) {
				node1 = GMT_Get_Index(API, grid1->header, row, col);
				range = round((near_range_rel[0] + col * grid1->header->inc[GMT_X] * dr) / dec);
				azimuth = round(
				    (clock_start_rel[0] * vel[0] + (grid1->header->ny - row) * grid1->header->inc[GMT_Y] * vel[0] / prf[0]) /
				    dec);
				node6 = GMT_Get_Index(API, grid6->header, azimuth, range);
				grid1->data[node1] = grid6->data[node6];
			}
		}

		/* subswath 2 */
		for (col = 0; col < grid2->header->nx; col++) {
			for (row = 0; row < grid2->header->ny; row++) {
				node2 = GMT_Get_Index(API, grid2->header, row, col);
				range = round((near_range_rel[1] + col * grid2->header->inc[GMT_X] * dr) / dec);
				azimuth = round(
				    (clock_start_rel[1] * vel[1] + (grid2->header->ny - row) * grid2->header->inc[GMT_Y] * vel[1] / prf[1]) /
				    dec);
				node6 = GMT_Get_Index(API, grid6->header, azimuth, range);
				grid2->data[node2] = grid6->data[node6];
			}
		}

		/* subswath 3 */
		for (col = 0; col < grid3->header->nx; col++) {
			for (row = 0; row < grid3->header->ny; row++) {
				node3 = GMT_Get_Index(API, grid3->header, row, col);
				range = round((near_range_rel[2] + col * grid3->header->inc[GMT_X] * dr) / dec);
				azimuth = round(
				    (clock_start_rel[2] * vel[2] + (grid3->header->ny - row) * grid3->header->inc[GMT_Y] * vel[2] / prf[2]) /
				    dec);
				node6 = GMT_Get_Index(API, grid6->header, azimuth, range);
				grid3->data[node3] = grid6->data[node6];
			}
		}

		/* subswath 4 */
		for (col = 0; col < grid4->header->nx; col++) {
			for (row = 0; row < grid4->header->ny; row++) {
				node4 = GMT_Get_Index(API, grid4->header, row, col);
				range = round((near_range_rel[3] + col * grid4->header->inc[GMT_X] * dr) / dec);
				azimuth = round(
				    (clock_start_rel[3] * vel[3] + (grid4->header->ny - row) * grid4->header->inc[GMT_Y] * vel[3] / prf[3]) /
				    dec);
				node6 = GMT_Get_Index(API, grid6->header, azimuth, range);
				grid4->data[node4] = grid6->data[node6];
			}
		}

		/* subswath 5 */
		for (col = 0; col < grid5->header->nx; col++) {
			for (row = 0; row < grid5->header->ny; row++) {
				node5 = GMT_Get_Index(API, grid5->header, row, col);
				range = round((near_range_rel[4] + col * grid5->header->inc[GMT_X] * dr) / dec);
				azimuth = round(
				    (clock_start_rel[4] * vel[4] + (grid5->header->ny - row) * grid5->header->inc[GMT_Y] * vel[4] / prf[4]) /
				    dec);
				node6 = GMT_Get_Index(API, grid6->header, azimuth, range);
				grid5->data[node5] = grid6->data[node6];
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
	}

	if (GMT_Destroy_Session(API))
		return EXIT_FAILURE; /* Remove the GMT machinery */

	return (EXIT_SUCCESS);
}
