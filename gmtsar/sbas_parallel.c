/*      $Id: sbas.h 39 2016-06-18 03/16/24 Xiaohua Xu $  */
/*****************************************************************************************
 *  Program to do InSAR time-series analysis. * Use Small Baseline Subset (SBAS)         *
 *  algorithm.                                                                           *
 *                                                                                       *
 *  Xiaohua Xu  Jan, 2021                                                                *
 *  Changing the code to parallel with openMP                                            *
 *                                                                                       *
 *  Anders Hogrelius Feb, 2022                                                           *
 *  Fixing memory allocation bugs and bugs in openMP code, optimizing code and           *
 *  adding code to place temporary files in /tmp directory (when running over NFS or     *
 *  iSCSI). Added disp variable to mmap.                                                 *
 *                                                                                       *
 *****************************************************************************************
 *  Program to do InSAR time-series analysis. * Use Small Baseline Subset (SBAS)         *
 *  algorithm.                                                                           *
 *                                                                                       *
 *  Xiaohua Xu and David T. Sandwell, Jul, 2016                                          *
 *                                                                                       *
 *  Taking the old sbas code to add atmospheric correction by means of common            *
 *  point stacking by Tymofeyeva & Fialko 2015.                                          *
 *                                                                                       *
 *  Xiaohua Xu and David Sandwell Jan 2021                                               *
 *  Adding option to mmap instead of using memory                                        *
 *                                                                                       *
 *****************************************************************************************
 **                                                                                     **
 **  Creator: Xiaopeng Tong and David Sandwell                                          **
 **  (Scripps Institution of Oceanography) * Date: 12/23/2012                           **
 **                                                                                     **
 *****************************************************************************************
 *  Modification history:                                                                *
 *  08/31/2013 debug the program                                                         *
 *  03/20/2014 add DEM error, mean velocity                                              *
 *  03/22/2014 add correlation, use weighted least-squares                               *
 *  04/01/2014 add seasonal term                                                         *
 *  08/19/2014 allocate memory with 1D array instead of multiple malloc                  *
 *  08/19/2014 do not require the velocity curve go through origin * 08/19/2014          *
 *  remove seasonal term                                                                 *
 *  08/19/2014 fix temporal smoothing                                                    *
 ****************************************************************************************/

/* Reference:
P. Berardino, G. Fornaro, R. Lanari, and E. Sansosti, “A new algorithm
for surface deformation monitoring based on small baseline differential
SAR interferograms,” IEEE Trans. Geosci. Remote Sensing, vol. 40, pp.
2375–2383, Nov. 2002.

Schmidt, D. A., and R. Bürgmann(2003),
Time-dependent land uplift and subsidence in the Santa Clara valley,
California, from a large interferometric synthetic aperture radar data set,
J. Geophys. Res., 108, 2416, doi:10.1029/2002JB002267, B9.
*/

/* Use DGELSY to solve the equations */
/* Calling DGELSY using column-major order */


#include "gmtsar.h"
#include<stdint.h>
#include<omp.h>
#include <sys/mman.h>
#include <fcntl.h>
//#define Malloc(type, n) (type *)malloc((n) * sizeof(type))

// It is better to use calloc as this allocates data in contigous memory blocks and zeroes them out on allocation
#define Malloc(type, n) (type *)calloc((n) , sizeof(type))

#define max(a, b) (((a) > (b)) ? (a) : (b))
#ifdef DEBUG
#define checkpoint() printf("Checkpoint64_t at line %lld in file %s\n", __LINE__, __FILE__)
#else
#define checkpoint()
#endif


char *USAGE = "USAGE: sbas_parallel intf.tab scene.tab N S xdim ydim [-atm ni] [-smooth sf] "
              "[-wavelength wl] [-incidence inc] [-range -rng] [-rms] [-dem]\n\n"
              " input: \n"
              "  intf.tab             --  list of unwrapped (filtered) interferograms:\n"
              "   format:   unwrap.grd  corr.grd  ref_id  rep_id  B_perp \n"
              "  scene.tab            --  list of the SAR scenes in chronological order\n"
              "   format:   scene_id   number_of_days \n"
              "   note:     the number_of_days is relative to a reference date \n"
              "  N                    --  number of the interferograms\n"
              "  S                    --  number of the SAR scenes \n"
              "  xdim and ydim        --  dimension of the interferograms\n"
              "  -smooth sf           --  smoothing factors, default=0 \n"
              "  -atm ni              --  number of iterations for atmospheric "
              "correction, default=0(skip atm correction) \n"
              "  -wavelength wl       --  wavelength of the radar wave (m) default=0.236 "
              "\n"
              "  -incidence theta     --  incidence angle of the radar wave (degree) "
              "default=37 \n"
              "  -range rng           --  range distance from the radar to the center of "
              "the interferogram (m) default=866000 \n"
              "  -rms                 --  output velocity uncertainty grids (mm/yr): "
              "rms.grd\n"
              "  -dem                 --  output DEM error (m): dem.grd \n"
              "  -mmap                --  use mmap to allocate disk space for less use of memory \n"
              "  -robust        --  only work with -atm turnned on, estimate velocity with records that has atm correction\n\n"
              " output: \n"
              "  disp_##.grd          --  cumulative displacement time series (mm) "
              "grids\n"
              "  vel.grd              --  linear regressed velocity (mm/yr) grids \n\n"
              " example:\n"
              "  sbas intf.tab scene.tab 88 28 700 1000 \n\n"
              "Before run: set the environment variable OMP_NUM_THREADS to a proper number"
              " that fits your machine\n\n"
              "REFERENCES: \n"
              "Berardino P., G. Fornaro, R. Lanari, and E. Sansosti, “A new algorithm "
              "for surface deformation monitoring based on small baseline differential "
              "SAR interferograms,” IEEE Trans. Geosci. Remote Sensing, vol. 40, pp. "
              "2375–2383, Nov. 2002. \n\n"
              "Schmidt, D. A., and R. Bürgmann 2003, Time-dependent land uplift and "
              "subsidence in the Santa Clara valley, California, from a large "
              "interferometric synthetic aperture radar data set, J. Geophys. Res., 108, "
              "2416, doi:10.1029/2002JB002267, B9. \n\n"
              "Tong, X. and Schmidt, D., 2016. Active movement of the Cascade landslide "
              "complex in Washington from a coherence-based InSAR time series method. "
              "Remote Sensing of Environment, 186, pp.405-415. \n\n"
              "Tymofyeyeva, E. and Fialko, Y., 2015. Mitigation of atmospheric phase "
              "delays in InSAR data, with application to the eastern California shear "
              "zone. Journal of Geophysical Research: Solid Earth, 120(8), "
              "pp.5952-5963.\n\n"
              "Xu, X., Sandwell, D. T., Tymofyeyeva, E., González-Ortega, A., & Tong, X."
              "(2017). Tectonic and anthropogenic deformation at the Cerro Prieto"
              "geothermal step-over revealed by Sentinel-1A InSAR. IEEE Transactions on"
              "Geoscience and Remote Sensing, 55(9), 5284-5292.\n\n";

void dgelsy_(const int64_t *m, const int64_t *n, const int64_t *nrhs, double *G, const int64_t *lda, double *b,
             const int64_t *ldb, int64_t *jpvt, const double *rcond, const int64_t *rank, double *work, const int64_t *lwork,
             const int64_t *info);

int parse_command_ts(int64_t agc, char **agv, float *sf, double *wl, double *theta, double *rng, int64_t *flag_rms,
                     int64_t *flag_dem, int64_t *atm, int64_t *flag_mmap, int64_t *flag_robust) {

	int64_t i;

	for (i = 7; i < agc; i++) {
		if (!strcmp(agv[i], "-smooth")) {
			i++;
			if (i == agc)
				die("no option after -smooth!\n", "");
			*sf = atof(agv[i]);
		}
		else if (!strcmp(agv[i], "-wavelength")) {
			i++;
			if (i == agc)
				die("no option after -wavelength!\n", "");
			*wl = atof(agv[i]);
		}
		else if (!strcmp(agv[i], "-incidence")) {
			i++;
			if (i == agc)
				die("no option after -incidence! \n", "");
			*theta = atof(agv[i]);
		}
		else if (!strcmp(agv[i], "-range")) {
			i++;
			if (i == agc)
				die("no option after -range \n", "");
			*rng = atof(agv[i]);
		}
		else if (!strcmp(agv[i], "-rms")) {
			*flag_rms = 1;
			fprintf(stderr, "compute RMS misfit\n");
		}
		else if (!strcmp(agv[i], "-dem")) {
			*flag_dem = 1;
			fprintf(stderr, "compute DEM error\n");
		}
        else if (!strcmp(agv[i], "-mmap")) {
            *flag_mmap = 1;
            fprintf(stderr, "mmap disk space for less use of memory\n");
        }  
        else if (!strcmp(agv[i], "-robust")) {
            *flag_robust = 1;
            fprintf(stderr, "using robust velocity estimate (require -atm)\n");
        } 
		else if (!strcmp(agv[i], "-atm")) {
			i++;
			if (i == agc)
				die("no option after -atm \n", "");
			*atm = atoi(agv[i]);
		}
		else {
			fprintf(stderr, " %s *** option not recognized ***\n\n", agv[i]);
			fprintf(stderr, "%s", USAGE);
			exit(1);
		}
	}

	fprintf(stderr, "\nsetting smoothing to %7.3lf \n", *sf);
	fprintf(stderr, "setting radar wavelength to %7.3lf m \n", *wl);
	fprintf(stderr, "setting radar incidence angle to %7.3lf degree \n", *theta);
	fprintf(stderr, "setting range to %9.3lf m \n", *rng);

	return (1);
}

int allocate_memory_ts(int64_t **jpvt, double **work, double **d, double **ds, float **bperp, char ***gfile, char ***cfile,
                       int64_t **L, double **time, int64_t **H, double **G, double **A, double **Gs, int64_t **flag, float **dem,
                       float **res, float **vel, float **phi, float **var, float **disp, int64_t n, int64_t m, int64_t lwork,
                       int64_t ldb, int64_t N, int64_t S, int64_t xdim, int64_t ydim, int64_t **hit, int64_t flag_mmap, int64_t n_atm) {

	int64_t i;
	char **p1, **p2;
    double Mem;

    Mem = xdim*ydim*(N*2 + S)*4 + n*4 + lwork*8 + N*8 + N*4 + N*8*2 + 256*N*2 + S*4 + S*8 + N*2*4 + m*n*8*2 + xdim*ydim*4*3 + S*S*4 + S*8;
    if (n_atm != 0) {
      Mem = Mem + xdim*ydim*(N*4+S*4) + N*4 + S*8 + S*4;
    }   
    if (flag_mmap != 0) {
      Mem = xdim*ydim*S*4 + n*4 + lwork*8 + N*8 + N*4 + N*8*2 + 256*N*2 + S*4 + S*8 + N*2*4 + m*n*8*2 + xdim*ydim*4*3 + S*S*4;
      if (n_atm != 0) {
        Mem = Mem + xdim*ydim*S*4 + N*4 + S*8 + S*4;
      }   
    }   
    Mem = Mem/1024/1024/1024;
    printf("Required Memory Usage is %.6f GB ...\n", Mem);

	if ((*jpvt = Malloc(int64_t, n)) == NULL)
		die("memory allocation!", "jpvt");
	if ((*work = Malloc(double, lwork)) == NULL)
		die("memory allocation!", "work");
	if ((*d = Malloc(double, ldb)) == NULL)
		die("memory allocation!", "d");
	if ((*ds = Malloc(double, N)) == NULL)
		die("memory allocation!", "ds");
	if ((*bperp = Malloc(float, N)) == NULL)
		die("memory allocation!", "bperp");
	if ((*gfile = Malloc(char *, N)) == NULL)
		die("memory allocation!", "gfile");
	if ((*cfile = Malloc(char *, N)) == NULL)
		die("memory allocation!", "cfile");
	p1 = *gfile;
	p2 = *cfile;
	for (i = 0; i < N; i++) {
		if ((p1[i] = Malloc(char, 256)) == NULL)
			die("memory allocation!", "gfile[i]");
		if ((p2[i] = Malloc(char, 256)) == NULL)
			die("memory allocation!", "cfile[i]");
	}
	if ((*L = Malloc(int64_t, S)) == NULL)
		die("memory allocation!", "L");
	if ((*time = Malloc(double, S)) == NULL)
		die("memory allocation!", "time");
	if ((*H = Malloc(int64_t, N * 2)) == NULL)
		die("memory allocation!", "H");
	if ((*G = Malloc(double, m *n)) == NULL)
		die("memory allocation!", "G");
	if ((*A = Malloc(double, m *n)) == NULL)
		die("memory allocation!", "A");
	if ((*Gs = Malloc(double, N *n)) == NULL)
		die("memory allocation!", "Gs");
	if ((*flag = Malloc(int64_t, xdim * ydim)) == NULL)
		die("memory allocation!", "flag");
	if ((*dem = Malloc(float, xdim *ydim)) == NULL)
		die("memory allocation!", "dem");
	if ((*res = Malloc(float, xdim *ydim)) == NULL)
		die("memory allocation!", "res");
	if ((*vel = Malloc(float, xdim *ydim)) == NULL)
		die("memory allocation!", "vel");
    if (flag_mmap == 0) {
	    if ((*phi = Malloc(float, N *xdim *ydim)) == NULL)
		    die("memory allocation!", "phi");
	    if ((*var = Malloc(float, N *xdim *ydim)) == NULL)
		    die("memory allocation!", "var");
            if ((*disp = Malloc(float, S *xdim *ydim)) == NULL)
                    die("memory allocation!", "disp");
    }

	if ((*hit = Malloc(int64_t, S * S)) == NULL)
		die("memory allocation!", "hit");
	printf("Memory Allocation Successful...\n");
	return (1);
}

int init_array_ts(double *G, double *Gs, float *res, float *dem, float *disp, int64_t n, int64_t m, int64_t xdim, int64_t ydim,
                  int64_t N, int64_t S, int64_t flag_mmap) {

        int64_t i, j, k, p;

        for (i = 0; i < m * n; i++)
                G[i] = 0;
        for (i = 0; i < n * N; i++)
                Gs[i] = 0;

	//fprintf(stderr,"S=%d N=%d n=%d m=%d\n",S,N,n,m);

	for (k = 0; k < (xdim*ydim); k++){
		//fprintf(stderr,"k=%d\n",k);
		res[k] = 0;
        	dem[k] = 0;
	}
	if (flag_mmap == 1) {
		for (p = 0; p < (xdim*ydim*S); p++){
			//fprintf(stderr,"p=%d\n",p);
			disp[p] = 0;
		}
	}
	return (1);
}

int read_table_data_ts(void *API, FILE *infile, FILE *datefile, char **gfile, char **cfile, int64_t *H, float *bperp,
                       int64_t *flag, float *var, float *phi, int64_t S, int64_t N, int64_t xdim, int64_t ydim,
                       struct GMT_GRID **Out, int64_t *L, double *time) {

	char tmp1[200], tmp2[200], tmp3[200];
	int64_t i, j, k, xin, yin, indx;
	float *corin, *grdin;
	struct GMT_GRID *CC = NULL, *GG = NULL;

	printf("read table file ...\n");
	/* read in scene.tab */
	i = 0;
	while (fscanf(datefile, "%s %s", &tmp1[0], &tmp2[0]) == 2) {
		if (feof(datefile))
			break;
		L[i] = atoi(tmp1);
		time[i] = atof(tmp2);
		i = i + 1;
	}
	if (i != S)
		die("S and number of the SAR scenes don't match!", "");
	fprintf(stderr, "number of SAR scenes is %lld \n", S);
	for (i = S - 1; i >= 0; i--)
		time[i] = time[i] - time[0];

	/* read in intf.tab */
	i = 0;
	while (fscanf(infile, "%s %s %s %s %s", gfile[i], cfile[i], &tmp1[0], &tmp2[0], &tmp3[0]) == 5) {
		if (feof(infile))
			break;
		H[i * 2 + 0] = atoi(tmp1);
		H[i * 2 + 1] = atoi(tmp2);
		bperp[i] = atof(tmp3);
		i = i + 1;
	}
	if (i != N)
		die("N and number of interferograms don't match!", "");
	fprintf(stderr, "number of interferograms is %lld \n", N);

	/* read in N 2-dimensional grd file into 3D array */
	printf("read phase and correlation grids ...\n");
	for (i = 0; i < N; i++) {
		if ((CC = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, cfile[i], NULL)) ==
		    NULL)
			die("Can't open ", cfile[i]);
		xin = CC->header->n_columns;
		yin = CC->header->n_rows;
		if (xin != xdim || yin != ydim)
			die("dimension don't match!", cfile[i]);
		if ((GG = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, gfile[i], NULL)) ==
		    NULL)
			die("Can't open ", gfile[i]);
		xin = GG->header->n_columns;
		yin = GG->header->n_rows;
		if (xin != xdim || yin != ydim)
			die("dimension don't match!", gfile[i]);
		if ((CC = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, cfile[i], CC)) == NULL)
			die("Can't read ", cfile[i]);
		if ((GG = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, gfile[i], GG)) == NULL)
			die("Can't read ", gfile[i]);
		corin = CC->data;
		grdin = GG->data;
		for (k = 0; k < ydim; k++) {
			for (j = 0; j < xdim; j++) {
                indx = i * xdim * ydim + ydim * j + k;
				phi[indx] = grdin[j + k * xdim];
				if (isnan(grdin[j + k * xdim]) != 0) {
					flag[j + xdim] = 1;
				}
				if (corin[j + k * xdim] >= 1e-2 && corin[j + k * xdim] <= 0.99) {
					/* Rosen et al., 2000 IEEE */
					var[indx] =
					    sqrt((1.0 - corin[j + k * xdim] * corin[j + k * xdim]) / (corin[j + k * xdim] * corin[j + k * xdim]));
				}
				else if (corin[j + k * xdim] < 1e-2) {
					var[indx] = 99.99;
				}
				else {
					var[indx] = 0.1;
				}
			}
		}
		if (*Out == NULL && (*Out = GMT_Duplicate_Data(API, GMT_IS_GRID, GMT_DUPLICATE_DATA, CC)) == NULL)
			die("error creating output grid", ""); /* Do this once */
		if (GMT_Destroy_Data(API, &CC))
			die("error freeing data", cfile[i]);
		if (GMT_Destroy_Data(API, &GG))
			die("error freeing data", gfile[i]);
	}
	return (1);
}

int init_G_ts(double *G, double *Gs, int64_t N, int64_t S, int64_t m, int64_t n, int64_t *L, int64_t *H, double *time, float sf,
              float *bperp, double scale) {

	int64_t i, j;

	printf("fill the G matrix ...\n");
	for (i = 0; i < N; i++) {
		for (j = 0; j < S - 1; j++) {
			if (L[j] >= H[i * 2 + 0] && L[j] < H[i * 2 + 1]) {
				G[i + m * j] = 1;
				Gs[i + N * j] = 1;
			}
		}
		G[i + m * (n - 1)] = bperp[i] * scale;
		Gs[i + N * (n - 1)] = bperp[i] * scale;
	}

	for (i = 0; i < S - 2; i++) {
		G[i + N + i * m] = sf / (time[i + 1] - time[i]);
		G[i + N + (i + 1) * m] = -sf / (time[i + 2] - time[i + 1]);
	}
	return (1);
}

int64_t lsqlin_sov_ts(int64_t xdim, int64_t ydim, float *disp, float *vel, int64_t *flag, double *d, double *ds, double *time,
                      double *G, double *Gs, double *A, float *var, float *phi, int64_t N, int64_t S, int64_t m, int64_t n,
                      double *work, int64_t lwork, int64_t flag_dem, float *dem, int64_t flag_rms, float *res, int64_t *jpvt,
                      double wl, double *atm_rms, int64_t flag_robust) {

    // float new,old;
    int64_t lda, ldb, zz;
    int64_t count, numt;
    // omp variables
    int64_t from, to, tid;
    int64_t i, j, k, p, info = 0;
    int64_t rank = 0, nrhs = 1;
    double rcond = 1e-3;
    double sumxx, sumxy, sumx, sumy, sumyy, aa;
    double *ddd,*GGG;

    lda = max(1, m);
    ldb = max(1, max(m, n));
    count = 0;
    for (zz = 0; zz < S; zz++) 
        if (atm_rms[zz] != 0.0 && zz != 0 && zz != 1 && zz != S - 1 && zz != S - 2)
            count++;

    printf("run least-squares problem over %lld by %lld pixel (%lld) ...\n", xdim, ydim, count);

    //Get max number of threads on this system
    int64_t numthreads_max;
    numthreads_max=sysconf(_SC_NPROCESSORS_ONLN);

    //Allocate separate doubles per thread for workwork variable
    double* workwork[numthreads_max];
        for (i = 0; i < numthreads_max; i++)
            workwork[i] = (double*)malloc(sizeof(double)*lwork);

        // the segment below needs some cleaning, deleting non-useful variables, etc.
    #pragma omp parallel private(from,to,tid,numt,i,j,k,p,rank,nrhs,rcond,sumxx,sumxy,sumx,sumy,sumyy,aa,ddd,GGG)
    {
    ddd = (double *) calloc(ldb,sizeof(double));
    GGG = (double *) calloc(m*n,sizeof(double));
    rank = 0; nrhs = 1; rcond = 1e-3;
    
    numt = omp_get_num_threads();
    tid = omp_get_thread_num();
    from = (ydim/numt)*tid;
    to = (ydim/numt)*(tid+1)-1;
    if (tid == numt-1)
        to = ydim-1;
    fprintf(stderr,"Initialing thread %d of %d, running rows from %d to %d\n", tid+1, numt, from, to);

    //This pragma necessary to keep data for all threads private within for clause
    #pragma omp parallel for collapse(2)
	for (k = from; k <= to; k++) {
		for (j = 0; j < xdim; j++) {
			/* check the dummy value of the grd file */
			if (flag[j + xdim * k] != 1) {
				for (i = 0; i < m; i++) {
					ddd[i] = 0;
					if (i < N) {
						ddd[i] = (double)phi[i * xdim * ydim + ydim * j + k] / var[i * xdim * ydim + ydim * j + k];
					}
				}

				for (i = 0; i < m; i++) {
					for (p = 0; p < n; p++) {
						if (i < N) {
							GGG[i + m * p] = A[i + m * p] / var[i * xdim * ydim + ydim * j + k];
						}
						else {
							GGG[i + m * p] = A[i + m * p];
						}
					}
				}

				dgelsy_(&m, &n, &nrhs, GGG, &lda, ddd, &ldb, jpvt, &rcond, &rank, workwork[tid], &lwork, &info);

				if (info != 0)
					printf("warning! input has an illegal value\n");
				if (j == 0 && k == 0) {
					if (rank == n) {
						printf("matrix is full rank: %lld\n\n", rank);
					}
					else if (rank < n) {
						printf("matrix is rank-deficient: %lld\n\n", rank);
						if (rank == 0) {
							printf("WARNING: rank is zero. Check scene.tab and intf.tab for "
							       "possible duplicates.\n\n");
							for (i = 0; i < N; i++) {
								for (p = 0; p < S - 1; p++) {
									printf("%.2f ", GGG[i + S * p]);
								}
								printf("\n");
							}
						}
					}
				}

				for (i = 0; i < S; i++) {
					for (p = 0; p < i; p++) {
						disp[i * xdim * ydim + j * ydim + k] = disp[i * xdim * ydim + j * ydim + k] + ddd[p];
					}
					// disp[i*xdim*ydim+j*ydim+k]=-79.58*wl*disp[i*xdim*ydim+j*ydim+k];
					// //1000/4/pi
				}

				if (flag_dem == 1)
					dem[j + xdim * k] = ddd[n - 1];
				/*
				                                if (flag_rms == 1) {
				                                        new=0;
				                                        old=0;
				                                // check the WRMS reduction
				                                        for (i=0;i<N;i++) {
				                                                pred=0;
				                                                for (p=0;p<n;p++) {
				                                                        pred=pred+d[p]*Gs[i+p*N]/var[i*xdim*ydim+ydim*j+k];
				                                                }
				                                                new=new+(ds[i]-pred)*(ds[i]-pred);
				                                                old=old+ds[i]*ds[i];
				                                        }
				                                        res[j+xdim*k]=(sqrt(old)-sqrt(new))/sqrt(old);

				                                        for (i=0;i<N;i++) {


				                                        }
				                                }
				*/
				// fitting a straight line by least squares
				sumxy = 0;
				sumxx = 0;
				sumx = 0;
				sumy = 0;
				sumyy = 0;
				if (count > 2 && flag_robust == 1) {
					for (i = 2; i < S - 2; i++) {
						if (atm_rms[i] != 0.0) {
							sumxy = sumxy + time[i] * disp[i * xdim * ydim + j * ydim + k];
							sumxx = sumxx + time[i] * time[i];
							sumy = sumy + disp[i * xdim * ydim + j * ydim + k];
							sumx = sumx + time[i];
						}
					}
					vel[j + xdim * k] = -79.58 * wl * (count * sumxy - sumx * sumy) / (count * sumxx - sumx * sumx) * 365.0;
					if (flag_rms == 1) {
						aa = sumy / count - (count * sumxy - sumx * sumy) / (count * sumxx - sumx * sumx) * sumx / count;
						for (i = 2; i < S - 2; i++) {
							if (atm_rms[i] != 0.0) {
								sumyy = sumyy + pow((disp[i * xdim * ydim + j * ydim + k] -
								                     time[i] * vel[j + xdim * k] / (-79.58 * wl * 365) - aa),
								                    2);
							}
						}
						res[j + xdim * k] =
						    sqrt(count * sumyy / ((count - 2) * (count * sumxx - sumx * sumx))) * (79.58 * wl * 365);
					}
				}
				else {
					for (i = 0; i < S; i++) {
						sumxy = sumxy + time[i] * disp[i * xdim * ydim + j * ydim + k];
						sumxx = sumxx + time[i] * time[i];
						sumy = sumy + disp[i * xdim * ydim + j * ydim + k];
						sumx = sumx + time[i];
					}
					vel[j + xdim * k] = -79.58 * wl * (S * sumxy - sumx * sumy) / (S * sumxx - sumx * sumx) * 365.0;
					if (flag_rms == 1) {
						aa = sumy / S - (S * sumxy - sumx * sumy) / (S * sumxx - sumx * sumx) * sumx / S;
						for (i = 2; i < S - 2; i++) {
							sumyy = sumyy + pow((disp[i * xdim * ydim + j * ydim + k] -
							                     time[i] * vel[j + xdim * k] / (-79.58 * wl * 365) - aa),
							                    2);
						}
						res[j + xdim * k] = sqrt(S * sumyy / ((S - 2) * (S * sumxx - sumx * sumx))) * (79.58 * wl * 365);
					}
				}
			}
			else {
				for (i = 0; i < S; i++)
					disp[i * xdim * ydim + j * ydim + k] = NAN;
				vel[j + xdim * k] = NAN;
				if (flag_rms == 1)
					res[j + xdim * k] = NAN;
				if (flag_dem == 1)
					dem[j + xdim * k] = NAN;
			}
		}
	}
    free(ddd); 
    free(GGG);
	}
	//Free up array of workwork doubles 
	for (i = 0; i < numthreads_max; i++)
        	free(workwork[i]);

	return (1);
}

int write_output_ts(void *API, struct GMT_GRID *Out, int64_t agc, char **agv, int64_t xdim, int64_t ydim, int64_t S,
                    int64_t flag_rms, int64_t flag_dem, float *disp, float *vel, float *res, float *dem, float *screen, double wl,
                    int64_t n_atm, int64_t *L) {

	int64_t i, j, k;
	float *grdin, *save_grid;
	char tmp1[200], outfile[200];

	printf("write output ...\n");

	strcpy(Out->header->title, "");
	strcpy(Out->header->remark, "");

	strcpy(tmp1, "");
	for (i = 0; i < agc; i++) {
		strcat(tmp1, agv[i]);
		strcat(tmp1, " ");
	}
	strcpy(Out->header->command, "");
	if (GMT_Set_Comment(API, GMT_IS_GRID, GMT_COMMENT_IS_COMMAND, tmp1, Out))
		die("could not set title", "");

	// strcpy(Out->header->title,"displacement time series (mm)");
	grdin = Out->data;
	for (i = 0; i < S; i++) {
		strcpy(outfile, "disp_");
		for (k = 0; k < ydim; k++) {
			for (j = 0; j < xdim; j++) {
				// disp[i*xdim*ydim+j*ydim+k]=-79.58*wl*disp[i*xdim*ydim+j*ydim+k];
				// //1000/4/pi
				grdin[j + k * xdim] = -79.58 * wl * disp[i * xdim * ydim + j * ydim + k];
			}
		}
		sprintf(tmp1, "%07lld", L[i]);
		strcat(outfile, tmp1);
		strcat(outfile, ".grd");
		sprintf(tmp1, "Displacement Time Series %03lld", i + 1);
		strcpy(Out->header->title, "");
		strcpy(Out->header->remark, "");
		if (GMT_Set_Comment(API, GMT_IS_GRID, GMT_COMMENT_IS_REMARK, "Displacement Time Series (mm)", Out))
			die("could not set title", "");
		if (GMT_Set_Comment(API, GMT_IS_GRID, GMT_COMMENT_IS_TITLE, tmp1, Out))
			die("could not set title", "");
		if (GMT_Write_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, outfile, Out)) {
			die("Failed to write output grid", outfile);
		}
	}
	save_grid = Out->data; /* Remember the original grid */
	if (flag_rms == 1) {
		Out->data = res;
		sprintf(outfile, "rms.grd");
		strcpy(Out->header->title, "");
		strcpy(Out->header->remark, "");
		if (GMT_Set_Comment(API, GMT_IS_GRID, GMT_COMMENT_IS_REMARK, "WRMS reduction from SBAS (mm)", Out))
			die("could not set title", "");
		if (GMT_Set_Comment(API, GMT_IS_GRID, GMT_COMMENT_IS_TITLE, "Weighed Root Mean Square of Fitting", Out))
			die("could not set title", "");
		if (GMT_Write_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, outfile, Out)) {
			die("Failed to write output grid", outfile);
		}
	}

	if (flag_dem == 1) {
		Out->data = dem;
		sprintf(outfile, "dem.grd");
		strcpy(Out->header->title, "");
		strcpy(Out->header->remark, "");
		if (GMT_Set_Comment(API, GMT_IS_GRID, GMT_COMMENT_IS_REMARK, "DEM error estimated from SBAS (m)", Out))
			die("could not set title", "");
		if (GMT_Set_Comment(API, GMT_IS_GRID, GMT_COMMENT_IS_TITLE, "Digital Elevation Model Error", Out))
			die("could not set title", "");
		if (GMT_Write_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, outfile, Out)) {
			die("Failed to write output grid", outfile);
		}
	}
	if (n_atm != 0) {
		grdin = Out->data;
		for (i = 0; i < S; i++) {
			strcpy(outfile, "aps_");
			for (k = 0; k < ydim; k++) {
				for (j = 0; j < xdim; j++) {
					grdin[j + k * xdim] = screen[i * xdim * ydim + j * ydim + k];
				}
			}
			sprintf(tmp1, "%07lld", L[i]);
			strcat(outfile, tmp1);
			strcat(outfile, ".grd");
			sprintf(tmp1, "Atmospheric Phase Screen %03lld", i + 1);
			strcpy(Out->header->title, "");
			strcpy(Out->header->remark, "");
			if (GMT_Set_Comment(API, GMT_IS_GRID, GMT_COMMENT_IS_REMARK, "Atmospheric Phase Screen", Out))
				die("could not set title", "");
			if (GMT_Set_Comment(API, GMT_IS_GRID, GMT_COMMENT_IS_TITLE, tmp1, Out))
				die("could not set title", "");
			if (GMT_Write_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, outfile, Out)) {
				die("Failed to write output grid", outfile);
			}
		}
	}

	Out->data = vel;
	sprintf(outfile, "vel.grd");
	strcpy(Out->header->title, "");
	strcpy(Out->header->remark, "");
	if (GMT_Set_Comment(API, GMT_IS_GRID, GMT_COMMENT_IS_REMARK, "Mean LOS velocity from SBAS (mm/yr)", Out))
		die("could not set title", "");
	if (GMT_Set_Comment(API, GMT_IS_GRID, GMT_COMMENT_IS_TITLE, "Mean Line-Of-Sight velocity from SBAS", Out))
		die("could not set title", "");
	if (GMT_Write_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, outfile, Out)) {
		die("Failed to write output grid", outfile);
	}

	Out->data = save_grid; /* Restore the pointer before destroying */

	return (1);
}

int free_memory_ts(int64_t N, float *phi, float *var, char **gfile, char **cfile, float *disp, double *G, double *A, double *Gs,
                   int64_t *H, double *d, double *ds, int64_t *L, float *res, float *vel, double *time, int64_t *flag,
                   float *bperp, float *dem, double *work, int64_t *jpvt, int64_t *hit, int64_t flag_mmap) {

	int64_t i;

	for (i = 0; i < N; i++) {
		free(gfile[i]);
		free(cfile[i]);
	}
	free(gfile);
	free(cfile);
	free(G);
	free(A);
	free(Gs);
	free(H);
	free(d);
	free(ds);
	free(L);
	free(res);
	free(vel);
	free(time);
	free(flag);
	free(bperp);
	free(dem);
	free(work);
	free(jpvt);
	free(hit);
    //this memory could be allocated with mmap
    if (flag_mmap == 0) {
	    free(phi);
	    free(var);
	    free(disp);
    }

	return (1);
}

int sum_intfs(float *phi, int64_t *mark, float *screen, int64_t xdim, int64_t ydim, int64_t N) {

	int64_t n, i, j, sum = 0;

	for (i = 0; i < N; i++)
		sum += llabs(mark[i]);
	// fprintf(stderr,"%lld sums       ",sum);

	for (i = 0; i < ydim; i++) {
		for (j = 0; j < xdim; j++) {
			screen[i * xdim + j] = 0.0;
		}
	}

	if (sum != 0) {
		for (n = 0; n < N; n++) {
			if (mark[n] == 0)
				continue;
			for (i = 0; i < ydim; i++) {
				for (j = 0; j < xdim; j++) {
					screen[i * xdim + j] += -phi[n * xdim * ydim + i * xdim + j] * (float)mark[n] / (float)sum;
				}
			}
		}
	}

	return (1);
}

int connect(int64_t *L, int64_t *H, double *time, int64_t *hit, int64_t *mark, int64_t N, int64_t S, int64_t n, int64_t mode) {

	// mode = 0 for all connections, mode = 1 for even connections

	int64_t i, j; //,sum = 0;
	for (i = 0; i < N; i++)
		mark[i] = 0;

	for (i = 0; i < S; i++) {
		if (hit[i * S + n] == 1) {
			for (j = 0; j < N; j++) {
				if (L[i] == H[2 * j] && L[n] == H[2 * j + 1])
					mark[j] = -1;
			}
		}
		if (hit[n * S + i] == 1) {
			for (j = 0; j < N; j++) {
				if (L[n] == H[2 * j] && L[i] == H[2 * j + 1])
					mark[j] = 1;
			}
		}
	}

	/* trim the result */
	if (mode == 1) {
		for (i = 0; i < S; i++) {
			if (hit[i * S + n] == 1) {
				for (j = n; j < S; j++) {
					if (hit[n * S + j] == 1 && fabs((time[n] - time[i]) - (time[j] - time[n])) < 5)
						break;
				}
				if (j == S) {
					for (j = 0; j < N; j++)
						if (L[i] == H[2 * j] && L[n] == H[2 * j + 1])
							mark[j] = 0;
				}
			}
			if (hit[n * S + i] == 1) {
				for (j = 0; j < n; j++) {
					if (hit[j * S + n] == 1 && fabs((time[n] - time[j]) - (time[i] - time[n])) < 5)
						break;
				}
				if (j == n) {
					for (j = 0; j < N; j++)
						if (L[n] == H[2 * j] && L[i] == H[2 * j + 1])
							mark[j] = 0;
				}
			}
		}
	}

	/*
	        for (i=0;i<N;i++) sum += mark[i];
	        while (sum != 0) {
	                if (sum > 0) {
	                        j=N-1;
	                        while (j>=0 && mark[j] == 0) j--;
	                        if (j<0) die("incorrect process of trimming connection
	   verctor",""); mark[j] = 0; sum -= 1;
	                }
	                if (sum < 0) {
	                        j=0;
	                        while (j<N && mark[j] == 0) j++;
	                        if (j>=N) die("incorrect process of trimming
	   connection verctor",""); mark[j] = 0; sum += 1;
	                }

	        }
	*/

	return (1);
}

double compute_noise(float *screen, int64_t xdim, int64_t ydim) {

	int64_t i, j, n;
	double sum = 0.0, rms = 0.0;

	n = 0;
	for (i = 0; i < ydim; i++) {
		for (j = 0; j < xdim; j++) {
			if (isnan(screen[i * xdim + j]) == 0) {
				sum = sum + (double)screen[i * xdim + j];
				n++;
			}
		}
	}

	if (sum != 0) {
		sum = sum / (double)n;

		for (i = 0; i < ydim; i++) {
			for (j = 0; j < xdim; j++) {
				if (isnan(screen[i * xdim + j]) == 0) {
					rms = rms + ((double)screen[i * xdim + j] - sum) * ((double)screen[i * xdim + j] - sum);
				}
			}
		}

		rms = sqrt(rms / (double)n);
	}
	// fprintf(stderr,"num points calculated: %lld\n",n);
	return (rms);
}

int apply_screen(float *screen, float *phi, int64_t xdim, int64_t ydim, int64_t N, int64_t *mark) {
	// also correct for the ones not used in estimation of aps.
	int64_t i, j, n;

	for (n = 0; n < N; n++) {
		if (mark[n] != 0) {
			// fprintf(stderr,"applying atm screen to intf %lld...\n",n);
			for (i = 0; i < ydim; i++) {
				for (j = 0; j < xdim; j++) {
					phi[n * xdim * ydim + i * xdim + j] = phi[n * xdim * ydim + i * xdim + j] + screen[i * xdim + j] * mark[n];
				}
			}
		}
	}

	return (1);
}

int remove_ts(float *phi, float *ts, int64_t xdim, int64_t ydim, int64_t N, int64_t S, int64_t *H, int64_t *L) {

	int64_t i, j, n, h1, h2;

	for (n = 0; n < N; n++) {
		for (i = 0; i < S; i++) {
			if (H[2 * n] == L[i])
				h1 = i;
			if (H[2 * n + 1] == L[i])
				h2 = i;
		}
		// fprintf(stderr,"removing deformation for pair
		// %lld_%lld...\n",L[h1],L[h2]);
		for (i = 0; i < ydim; i++) {
			for (j = 0; j < xdim; j++) {
				phi[n * xdim * ydim + i * xdim + j] = phi[n * xdim * ydim + i * xdim + j] - ts[h2 * xdim * ydim + i * xdim + j] +
				                                      ts[h1 * xdim * ydim + i * xdim + j];
			}
		}
	}

	return (1);
}

int rank_double(double *nums, int64_t *seq, int64_t n) {

	// rank the numbers, seq will record the sequence
	int64_t i, j, recn, *mark;
	double *nums2, rec;

	nums2 = (double *)malloc(n * sizeof(double));
	mark = (int64_t *)malloc(n * sizeof(int64_t));

	for (i = 0; i < n; i++)
		nums2[i] = nums[i];
	for (i = 0; i < n; i++) {
		if (nums2[i] == 0.0)
			mark[i] = 1;
		else
			mark[i] = 0;
		// fprintf(stderr,"%lld ",mark[i]);
	}
	// fprintf(stderr,"\n\n");

	for (i = 0; i < n; i++) {
		rec = 0.0;
		for (j = 0; j < n; j++) {
			if (rec < fabs(nums2[j])) {
				rec = fabs(nums2[j]);
				recn = j;
			}
		}
		if (rec != 0.0) {
			seq[i] = recn;
			nums2[recn] = 0.0;
		}
		else {
			break;
		}
	}

	recn = 0;
	for (i = 0; i < n; i++)
		if (mark[i] != 0)
			recn++;
	for (i = 0; i < recn; i++) {
		for (j = 0; j < n; j++) {
			if (mark[j] != 0) {
				seq[n - recn + i] = j;
				mark[j] = 0;
				break;
			}
		}
	}

	return (1);
}

int main(int argc, char **argv) {

	/* define variables */
	double EE = 2.718281828459046;

	char **gfile = NULL, **cfile = NULL;
	int64_t i, j, m, n, nrhs = 1, xdim, lwork, ydim, k1, k2;
	int64_t N, S;
	int64_t ldb, lda, *flag = NULL, *jpvt = NULL, *H = NULL, *L = NULL, *hit = NULL, *mark = NULL;
	int64_t flag_rms = 0, flag_dem = 0, flag_mmap = 0, flag_robust = 0;
	float *phi = NULL, *tmp_phi = NULL, sf, *disp = NULL, *res = NULL, *dem = NULL, *bperp = NULL, *vel = NULL, *screen = NULL,
	      *tmp_screen = NULL;
	float *var = NULL;
	double *G = NULL, *A = NULL, *Gs = NULL, *d = NULL, *ds = NULL;
	double *work = NULL, *time = NULL;
	double rng, wl, theta, scale, aa, bb;
	FILE *infile = NULL, *datefile = NULL;
	void *API = NULL;            /* GMT control structure */
	struct GMT_GRID *Out = NULL; /* For the output grid */

	double *atm_rms;
	int64_t *atm_rank, n_atm = 0, kk;
	float *sfs;

	size_t mm_size_N;
	size_t mm_size_S;
	int ftmp_phi = 0, fphi = 0, fvar = 0, fdisp = 0;

	char *sz_tmp_sbas_phi;
	char *sz_tmp_sbas_var;
	char *sz_tmp_sbas_disp;
	char *sz_tmp_sbas_tmp_phi;

	if (argc < 7)
		die("\n", USAGE);

        asprintf(&sz_tmp_sbas_phi,"./%d_tmp_sbas_phi.%d",geteuid(),getpid());
        asprintf(&sz_tmp_sbas_var,"./%d_tmp_sbas_var.%d",geteuid(),getpid());
        asprintf(&sz_tmp_sbas_disp,"./%d_tmp_sbas_disp.%d",geteuid(),getpid());
        asprintf(&sz_tmp_sbas_tmp_phi,"./%d_tmp_sbas_tmp_phi.%d",geteuid(),getpid());


        if (sz_tmp_sbas_phi == NULL || sz_tmp_sbas_var == NULL || sz_tmp_sbas_disp == NULL || sz_tmp_sbas_tmp_phi == NULL) {
                fprintf(stderr, "Error in asprintf when allocation memory.\n");
                return EXIT_FAILURE;
        }


	/* Begin: Initializing new GMT session */
	if ((API = GMT_Create_Session(argv[0], 0U, 0U, NULL)) == NULL)
		return EXIT_FAILURE;

	/* default parameters are for ALOS-1 */
	/* use range at the center of the image*/
	rng = 866000;
	/* wavelength of the radar wave*/
	wl = 0.236;
	/* incidence angle at the center of the image*/
	theta = 37;
	/* smoothing factor */
	sf = 0;

	/* reading in some parameters and open corresponding files */
	if ((infile = fopen(argv[1], "r")) == NULL)
		die("Can't open file", argv[1]);
	if ((datefile = fopen(argv[2], "r")) == NULL)
		die("Can't open file", argv[2]);
	N = (int64_t)atoi(argv[3]);
	S = (int64_t)atoi(argv[4]);
	xdim = (int64_t)atoi(argv[5]);
	ydim = (int64_t)atoi(argv[6]);

	fprintf(stderr, "\n");

	/* read in the parameters from command line */
	parse_command_ts(argc, argv, &sf, &wl, &theta, &rng, &flag_rms, &flag_dem, &n_atm, &flag_mmap, &flag_robust);

	/* setting up some parameters */
	scale = 4.0 * M_PI / wl / rng / sin(theta / 180.0 * M_PI);
	m = N + S - 2; // number of rows in the G matrix
	n = S;         // number of columns in the G matrix
	lwork = max(1, m * n + max(m * n, nrhs) * 16);
	lda = max(1, m);
	ldb = max(1, max(m, n));

	/* memory allocation */ // also malloc for atm(nx,ny,S), hit(N,S), sum_vec(N)
	                        // and atm_rms(S)
	mm_size_N = 4 * (size_t)N * (size_t)xdim * (size_t)ydim;
	mm_size_S = 4 * (size_t)S * (size_t)xdim * (size_t)ydim;
	allocate_memory_ts(&jpvt, &work, &d, &ds, &bperp, &gfile, &cfile, &L, &time, &H, &G, &A, &Gs, &flag, &dem, &res, &vel, &phi,
	                   &var, &disp, n, m, lwork, ldb, N, S, xdim, ydim, &hit, flag_mmap, n_atm);

        /* mmap the phi and var arrays. this must be done in the main program  */
    if (flag_mmap == 1) {
	    remove(sz_tmp_sbas_phi);
	    if ((fphi = open(sz_tmp_sbas_phi, O_RDWR | O_CREAT | O_EXCL, (mode_t)0755)) < 0)
		    die("can't open %s for reading", sz_tmp_sbas_phi);
	    lseek(fphi,mm_size_N-1, SEEK_SET);
	    write(fphi, "",1);
	    if ((phi = mmap( NULL , mm_size_N, PROT_READ | PROT_WRITE, MAP_SHARED, fphi, 0)) == MAP_FAILED)
		    die("mmap error for input", "phi");

	    remove(sz_tmp_sbas_var);
	    if ((fvar = open(sz_tmp_sbas_var, O_RDWR | O_CREAT | O_EXCL, (mode_t)0755)) < 0)
		    die("can't open %s for reading", sz_tmp_sbas_var);
	    lseek(fvar,mm_size_N-1, SEEK_SET);
	    write(fvar, "",1);
	    if ((var = mmap( NULL , mm_size_N, PROT_READ | PROT_WRITE, MAP_SHARED, fvar, 0)) == MAP_FAILED)
		    die("mmap error for input", "var");

            remove(sz_tmp_sbas_disp);
            if ((fdisp = open(sz_tmp_sbas_disp, O_RDWR | O_CREAT | O_EXCL, (mode_t)0755)) < 0)
                    die("can't open %s for reading", sz_tmp_sbas_disp);
            lseek(fdisp,mm_size_S-1, SEEK_SET);
            write(fdisp, "",1);
            if ((disp = mmap( NULL , mm_size_S, PROT_READ | PROT_WRITE, MAP_SHARED, fdisp, 0)) == MAP_FAILED)
                    die("mmap error for input", "disp");

    }

	// initialization
	//fprintf(stderr,"Before 1 call init_array_ts: xdim=%d ydim=%d N=%d S=%d\n",xdim, ydim, N, S);
	init_array_ts(G, Gs, res, dem, disp, n, m, xdim, ydim, N, S, flag_mmap);

	// reading in the table files
	read_table_data_ts(API, infile, datefile, gfile, cfile, H, bperp, flag, var, phi, S, N, xdim, ydim, &Out, L, time);

	printf("%.6f %.6f %.6f %.6f\n", sf, scale, time[0], bperp[0]);

	if (n_atm == 0) {
        flag_robust = 0;
		atm_rms = (double *)malloc(S * sizeof(double));
		for (i = 0; i < S; i++)
			atm_rms[i] = 0.0;

		init_G_ts(G, Gs, N, S, m, n, L, H, time, sf, bperp, scale);
		for (i = 0; i < m * n; i++)
			A[i] = G[i];
		for (i = 0; i < xdim * ydim * S; i++)
			disp[i] = 0.0;
		lsqlin_sov_ts(xdim, ydim, disp, vel, flag, d, ds, time, G, Gs, A, var, phi, N, S, m, n, work, lwork, flag_dem, dem,
		              flag_rms, res, jpvt, wl, atm_rms, flag_robust);
	}
	else {
		fprintf(stderr, "\n\nApplying atmospheric correction by common point stacking...\n\n");
		mark = (int64_t *)malloc(N * sizeof(int64_t));
		screen = (float *)malloc(xdim * ydim * sizeof(float) * S);
		atm_rms = (double *)malloc(S * sizeof(double));
		atm_rank = (int64_t *)malloc(S * sizeof(int64_t));
		sfs = (float *)malloc((n_atm + 2) * sizeof(float));
		tmp_screen = (float *)malloc(xdim * ydim * sizeof(float));

		/* allocate the memory of tmp_phi using mmap */
        if (flag_mmap == 0) {
		    tmp_phi = (float *)malloc(xdim * ydim * sizeof(float) * N);
        }
        else {
		    remove(sz_tmp_sbas_tmp_phi);
		    if ((ftmp_phi = open(sz_tmp_sbas_tmp_phi, O_RDWR | O_CREAT | O_EXCL, (mode_t)0755)) < 0)
			    die("can't open %s for reading", sz_tmp_sbas_tmp_phi);
		    lseek(ftmp_phi,mm_size_N-1, SEEK_SET);
		    write(ftmp_phi, "",1);
		    if ((tmp_phi = mmap( NULL , mm_size_N, PROT_READ | PROT_WRITE, MAP_SHARED, ftmp_phi, 0)) == MAP_FAILED)
			    die("mmap error for input", "tmp_phi");
        }

		sfs[0] = 1000.0;
		sfs[n_atm] = sf;

		if (n_atm >= 2) {
            if (sf > 0) {
			  bb = (log(1000) - log(sf)) / (double)(n_atm);
			  aa = pow(EE, log(sf) - bb);
            }
            else {
              bb = (log(1000) - log(0.01)) / (double)(n_atm);
              aa = pow(EE, log(0.01) - bb);
            }           
			for (i = 1; i <= n_atm - 1; i++)
				sfs[n_atm - i] = aa * pow(EE, bb * (double)(i + 1));
		}

		// get the hit matrix which records the pairs processed
		for (i = 0; i < S * S; i++)
			hit[i] = 0;
		fprintf(stderr, "\n\n\nHit Matrix:\n");
		for (i = 0; i < N; i++) {
			for (j = 0; j < S; j++) {
				if (H[i * 2] == L[j])
					k1 = j;
				if (H[i * 2 + 1] == L[j])
					k2 = j;
			}
			hit[k1 * S + k2] = 1;
		}

		for (i = 0; i < S; i++) {
			fprintf(stderr, "%lld ", L[i]);
			for (j = 0; j < S; j++) {
				fprintf(stderr, "%lld ", hit[i * S + j]);
			}
			fprintf(stderr, "\n");
		}
		fprintf(stderr, "\n");
		fprintf(stderr, "\n");

		// for debuging
		// example of connect
		//        connect(L,H,time,hit,mark,N,S,10,0);
		//        fprintf(stderr,"\n\n%lld ",L[10]);
		//        for (i=0;i<N;i++) {
		//                fprintf(stderr," %lld",mark[i]);
		//                if(mark[i] != 0)
		//                fprintf(stderr,"(%d)",(int)(floor((H[i*2+1]-H[i*2])/1000.0)*365+H[i*2+1]-floor(H[i*2+1]/1000.0)*1000-(H[i*2]-floor(H[i*2]/1000.0)*1000)));
		//        }
		//        fprintf(stderr,"\n");
		//        connect(L,H,time,hit,mark,N,S,10,1);
		//        fprintf(stderr,"\n\n%lld ",L[10]);
		//        for (i=0;i<N;i++) {
		//                fprintf(stderr," %lld",mark[i]);
		//                if(mark[i] != 0)
		//                fprintf(stderr,"(%d)",(int)(floor((H[i*2+1]-H[i*2])/1000.0)*365+H[i*2+1]-floor(H[i*2+1]/1000.0)*1000-(H[i*2]-floor(H[i*2]/1000.0)*1000)));
		//        }
		//        fprintf(stderr,"\n");

		// compute time series with tons of smoothing

		fprintf(stderr, "Applying exponential relaxation on smoothing parameters\n");
		for (i = 0; i < xdim * ydim * N; i++)
			tmp_phi[i] = phi[i];

		for (kk = 1; kk <= n_atm; kk++) {
			sf = sfs[kk - 1];
			fprintf(stderr, "\nSetting smoothing parameter to %f...\n", sf);
			init_G_ts(G, Gs, N, S, m, n, L, H, time, sf, bperp, scale);
			for (i = 0; i < m * n; i++)
				A[i] = G[i];
			for (i = 0; i < S; i++)
				atm_rms[i] = 0.0;

			fprintf(stderr, "Computing deformation time-series...\n");
			// progam below is paralleled
			lsqlin_sov_ts(xdim, ydim, disp, vel, flag, d, ds, time, G, Gs, A, var, tmp_phi, N, S, m, n, work, lwork, flag_dem,
			              dem, flag_rms, res, jpvt, wl, atm_rms, flag_robust);
			// remove the very smooth deformation signal from the data
			if (kk > 1)
				for (i = 0; i < xdim * ydim * N; i++)
					tmp_phi[i] = phi[i];
			fprintf(stderr, "Removing deformation time-series from original "
			                "unwrapped phase...\n");
			remove_ts(tmp_phi, disp, xdim, ydim, N, S, H, L);

			// compute atmospheric phase screen and the noise rms, 1st time, do not
			// update the phase during computation
			fprintf(stderr, "Computing atmospheric phase screen by common-point stacking...\n");
			if (kk == 1) {
				fprintf(stderr, "Initial estimate of APS...\n");

				for (i = 0; i < S; i++) {
					connect(L, H, time, hit, mark, N, S, i, 1);
					// compute atm with original interferograms
					sum_intfs(phi, mark, tmp_screen, xdim, ydim, N);
					atm_rms[i] = compute_noise(tmp_screen, xdim, ydim);
					for (j = 0; j < xdim * ydim; j++)
						screen[i * xdim * ydim + j] = tmp_screen[j];
					// apply_screen(tmp_screen,phi,xdim,ydim,N,mark);
				}

				rank_double(atm_rms, atm_rank, S);
				for (i = 0; i < S; i++)
					fprintf(stderr, "atm_noise(NO.%lld) = %lf\n ", atm_rank[i], atm_rms[atm_rank[i]]);
				fprintf(stderr, "\n\n");
			}

			// compute and apply aps and update as you go
			for (i = 0; i < S; i++) {
				connect(L, H, time, hit, mark, N, S, atm_rank[i], 1);
				sum_intfs(tmp_phi, mark, tmp_screen, xdim, ydim, N);
				atm_rms[atm_rank[i]] = compute_noise(tmp_screen, xdim, ydim);
				for (j = 0; j < xdim * ydim; j++)
					screen[atm_rank[i] * xdim * ydim + j] = tmp_screen[j];
				connect(L, H, time, hit, mark, N, S, atm_rank[i], 0);
				apply_screen(tmp_screen, tmp_phi, xdim, ydim, N, mark);
			}
			rank_double(atm_rms, atm_rank, S);
			for (i = 0; i < S; i++)
				fprintf(stderr, "atm_noise(NO.%lld) = %lf\n ", atm_rank[i], atm_rms[atm_rank[i]]);
			fprintf(stderr, "\n\n");

			// start agian with aps correction
			fprintf(stderr, "Applying atmospheric phase screen to original unwrapped phase...\n");
			for (i = 0; i < xdim * ydim * N; i++)
				tmp_phi[i] = phi[i];
			for (i = 0; i < S; i++) {
				connect(L, H, time, hit, mark, N, S, i, 0);
				for (j = 0; j < xdim * ydim; j++)
					tmp_screen[j] = screen[i * xdim * ydim + j];
				apply_screen(tmp_screen, tmp_phi, xdim, ydim, N, mark);
			}
		}

		// lastly compute time-series
		sf = sfs[n_atm];
		fprintf(stderr, "Setting smoothing parameter to %f...\n", sf);
		//fprintf(stderr,"Before 2 call init_array_ts: xdim=%d ydim=%d N=%d S=%d\n",xdim, ydim, N, S);
		init_array_ts(G, Gs, res, dem, disp, n, m, xdim, ydim, N, S, flag_mmap);
		init_G_ts(G, Gs, N, S, m, n, L, H, time, sf, bperp, scale);
		for (i = 0; i < m * n; i++)
			A[i] = G[i];
		for (i = 0; i < xdim * ydim * S; i++)
			disp[i] = 0.0;
		lsqlin_sov_ts(xdim, ydim, disp, vel, flag, d, ds, time, G, Gs, A, var, tmp_phi, N, S, m, n, work, lwork, flag_dem, dem,
		              flag_rms, res, jpvt, wl, atm_rms, flag_robust);
	}

	// write output
	write_output_ts(API, Out, argc, argv, xdim, ydim, S, flag_rms, flag_dem, disp, vel, res, dem, screen, wl, n_atm, L);

	/* free memory */

	free_memory_ts(N, phi, var, gfile, cfile, disp, G, A, Gs, H, d, ds, L, res, vel, time, flag, bperp, dem, work, jpvt, hit, flag_mmap);

	if (n_atm != 0) {
		free(mark);
		free(screen);
		free(tmp_screen);
		free(atm_rank);
        if (flag_mmap == 0) {
		    free(tmp_phi);
        }
        else {
		    munmap(tmp_phi, mm_size_N);
        }
		free(atm_rms);
	}
    if (flag_mmap == 1) {
		munmap(phi, mm_size_N);
		munmap(var, mm_size_N);
		munmap(disp, mm_size_S);
    }

	fclose(infile);
	fclose(datefile);
	close(fphi);
	close(fvar);
	close(fdisp);
	close(ftmp_phi);

	free(sz_tmp_sbas_phi);
	free(sz_tmp_sbas_var);
	free(sz_tmp_sbas_disp);
	free(sz_tmp_sbas_tmp_phi);

	if (GMT_Destroy_Session(API))
		return EXIT_FAILURE; /* Remove the GMT machinery */

	return (EXIT_SUCCESS);
}
