#include "gmtsar.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char *USAGE = "\nUsage: nearest_grid input.grd output.grd [search_radius]\n\n"
              "      NaNs will be interpolated to its nearest neighbour\n\n";

int create_grid(void *, char *, char *, int);

int main(int argc, char **argv) {

	void *API = NULL;
	int radius = 0;

	if (argc != 3 && argc != 4)
		die(USAGE, "");
	if (argc == 4) {
		radius = atoi(argv[3]);
		printf("Setting search radius to be %d ...\n", radius);
	}
	if ((API = GMT_Create_Session(argv[0], 0U, 0U, NULL)) == NULL)
		return EXIT_FAILURE;
	create_grid(API, argv[1], argv[2], radius);
	return (1);
}

int find_nearest(int i, int j, int *r2, int *is, int *js, int *xs, int *ys) {

	int ct = 0, nx, ny, nx1, ii, k = 0;
	int rr;

	rr = *r2 + 1e7;
	nx = (int)(sqrt((double)(*r2) / 2.0));
	for (nx1 = nx; nx1 <= (int)sqrt((double)(*r2)) + 1; nx1++) {
		if (nx1 * nx1 < *r2)
			ny = (int)(sqrt((*r2) - nx1 * nx1));
		else
			ny = 0;
		while (nx1 * nx1 + ny * ny <= (*r2) && ny <= nx1) {
			ny++;
		}
		if (ny <= nx1) {
			if (rr > nx1 * nx1 + ny * ny) {
				k = 0;
				rr = nx1 * nx1 + ny * ny;
				xs[k] = nx1;
				ys[k] = ny;
			}
			else if (rr == nx1 * nx1 + ny * ny) {
				k++;
				xs[k] = nx1;
				ys[k] = ny;
			}
		}
	}

	for (ii = 0; ii <= k; ii++) {
		nx = xs[ii];
		ny = ys[ii];

		if (ny == 0) {
			js[ct + 0] = 0;
			js[ct + 1] = 0;
			js[ct + 2] = nx;
			js[ct + 3] = -nx;
			is[ct + 0] = nx;
			is[ct + 1] = -nx;
			is[ct + 2] = 0;
			is[ct + 3] = 0;
			ct = ct + 4;
		}
		else if (nx != ny) {
			js[ct + 0] = ny;
			js[ct + 1] = ny;
			js[ct + 2] = -ny;
			js[ct + 3] = -ny;
			js[ct + 4] = nx;
			js[ct + 5] = -nx;
			js[ct + 6] = nx;
			js[ct + 7] = -nx;
			is[ct + 0] = nx;
			is[ct + 1] = -nx;
			is[ct + 2] = nx;
			is[ct + 3] = -nx;
			is[ct + 4] = ny;
			is[ct + 5] = ny;
			is[ct + 6] = -ny;
			is[ct + 7] = -ny;
			ct = ct + 8;
		}
		else {
			js[ct + 0] = nx;
			js[ct + 1] = nx;
			js[ct + 2] = -nx;
			js[ct + 3] = -nx;
			is[ct + 0] = nx;
			is[ct + 1] = -nx;
			is[ct + 2] = nx;
			is[ct + 3] = -nx;
			ct = ct + 4;
		}
	}
	/*
	    if (i==409 && j == 1422) {
	        fprintf(stderr,"Radius2: %d ",rr);
	        for (ii=0;ii<ct;ii++) fprintf(stderr,"(%d %d) ",is[ii],js[ii]);
	        fprintf(stderr,"\n");
	    }
	*/
	for (ii = 0; ii < ct; ii++) {
		is[ii] = is[ii] + i;
		js[ii] = js[ii] + j;
	}

	*r2 = rr;

	return (ct);
}

double nearest_interp(int nx, int ny, float *m, float *m_interp, int radius) {

	int i, j, flag, ct, k, kt, kk = 1, recx = 1, recy = 1;

	int *is, *js, *xs, *ys;

	int cs = 0;

	int rr;

	is = (int *)malloc(sizeof(int) * 2000);
	js = (int *)malloc(sizeof(int) * 2000);
	xs = (int *)malloc(sizeof(int) * 500);
	ys = (int *)malloc(sizeof(int) * 500);
	//    bufx = (int *)malloc(sizeof(int)*nx);
	//    bufy = (int *)malloc(sizeof(int)*nx);
	//    bufr = (int *)malloc(sizeof(int)*nx);
	//    idx = (int *)malloc(sizeof(int)*nx);
	//    idx2 = (int *)malloc(sizeof(int)*nx);
	/*
	    for (k=0;k<nx;k++) {
	        idx[k] = 0;
	        idx2[k] = 0;
	        bufx[k] = 0; bufy[k] = 0; bufr[k] = 0;
	    }
	*/
	printf("Interpolating to nearest neighbour...\n  ");
	fprintf(stderr, "Working on line  ");
	for (i = 0; i < ny; i++) {
		for (k = 0; k < kk; k++)
			fprintf(stderr, "\b");
		fprintf(stderr, "%d ...", i + 1);
		kk = floor(log10((double)(i + 1))) + 1 + 4;
		rr = 0;
		for (j = 0; j < nx; j++) {

			if (isnan(m[i * nx + j]) == 0) {
				m_interp[i * nx + j] = m[i * nx + j];
				rr = 0;
			}
			else {
				flag = 0;
				// if (rr >= 1) rr =
				// ((int)sqrt((double)rr)-1)*((int)sqrt((double)rr)-1);

				if (rr >= 4 && recy > 0 && recx > 0)
					rr = (recx - 1) * (recx - 1) + (recy - 1) * (recy - 1) - 1;
				// if (rr >= 4 && recy > 0) rr = ((int)sqrt(rr)-1)*((int)sqrt(rr)-1)-1;
				else if (rr >= 4 && recy == 0 && recx > 0)
					rr = (recx - 1) * (recx - 1) - 1;
				else if (rr >= 4 && recy > 0 && recx == 0)
					rr = (recy - 1) * (recy - 1) - 1;

				// else if (idx[j] != 0 && bufr[j] >= 4 && bufx[j] > 0) rr =
				// (bufx[j]-1)*(bufx[j]-1) + bufy[j]*bufy[j] -1; else if (idx[j] != 0 &&
				// bufr[j] >= 4 && bufx[j] == 0) rr = bufy[j]*bufy[j] -1;

				else
					rr = 0;

				// if (rr<0) rr = 0;
				// rr = 0;

				while (flag == 0 && rr <= (double)(radius * radius)) {
					ct = find_nearest(i, j, &rr, is, js, xs, ys);
					cs++;
					for (k = 0; k < ct; k++) {
						if (is[k] >= 0 && is[k] < ny && js[k] >= 0 && js[k] < nx) {
							if (isnan(m[is[k] * nx + js[k]]) == 0) {
								m_interp[i * nx + j] = m[is[k] * nx + js[k]];
								flag = 1;
								kt = k;
								recx = abs(is[k] - i);
								recy = abs(js[k] - j);
								//                                idx[j] = 1;
								//                                bufx[j] = abs(is[k]-i);
								//                                bufy[j] = abs(js[k]-j);
								//                                bufr[j] =
								//                                bufx[j]*bufx[j]+bufy[j]*bufy[j];
								break;
							}
						}
					}
				}
			}
			if (i == 0 && rr == -1)
				fprintf(stderr, "(%d %d %d %d)\n", j, rr, recx, recy);

			// if (i == 409 && j == 1422) fprintf(stderr,"(%.5f %d %d %.5f  %.6f %.6f
			// %.6f) ",m_interp[i*nx+j],is[kt]-i,js[kt]-j,m[is[kt]*nx+js[kt]],
			// m[0],m[1],m[2]);
		}
		//       for(k=0;k<nx;k++) {
		//           idx2[k] = idx[k];
		//           idx[k] = 0;
		//       }
	}
	fprintf(stderr, "\n");

	printf("%d number of searches used ...\n", cs);
	free(is);
	free(js);
	free(xs);
	free(ys);

	return (1.0);
}

int create_grid(void *API, char *file, char *output, int radius) {

	float *m, *m_interp;
	int i, j, nx, ny, rr;

	struct GMT_GRID *T = NULL, *OUT = NULL;

	rr = radius;

	printf("Reading in original grid...\n");
	if ((T = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, file, NULL)) == NULL)
		die("cannot open grdfile", file);
	if (GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, file, T) == NULL)
		return EXIT_FAILURE;

	nx = T->header->nx;
	ny = T->header->ny;
	// printf("nx: %d, ny: %d\n",nx,ny);
	m = T->data;
	m_interp = (float *)malloc(sizeof(float) * nx * ny);

	for (i = 0; i < ny; i++)
		for (j = 0; j < nx; j++)
			m_interp[i * nx + j] = m[i * nx + j];

	if (rr == 0) {
		rr = sqrt((double)(nx * nx + ny * ny)) + 1;
	}

	nearest_interp(nx, ny, m, m_interp, rr);

	printf("WRITING GRID IMAGE: Width x Heihgt = %d x %d...\n", nx, ny);
	if (OUT == NULL && (OUT = GMT_Duplicate_Data(API, GMT_IS_GRID, GMT_DUPLICATE_DATA, T)) == NULL)
		die("error creating output grid", "");
	OUT->data = m_interp;
	if (GMT_Write_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, output, OUT))
		die("Failed to write output grid", output);

	free(m_interp);

	return (1);
}
