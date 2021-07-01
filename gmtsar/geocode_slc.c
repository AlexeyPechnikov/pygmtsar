/***************************************************************************
 * geocode_slc reads in one complex SAR_SLC image and sample it to the     *
 * given DEM file. Propogation delay is geometrically computed and removed *
 ***************************************************************************/
/***************************************************************************
 * Creator:  Xiaohua Xu                                                    *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  6/7/21                                                        *
 ***************************************************************************/
#include "gmtsar.h"
#include "llt2xyz.h"
#include "orbit.h"
#include <fcntl.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/types.h>

char *USAGE = "geocode_slc [GMTSAR] - Sample slc to DEM and remove propogation delay\n\n"
              "Usage: "
              "geocode_slc your_file.PRM dem.grd \n"
              "(Put your .LED .SLC .PRM file in PWD)\n \n";

void fbisinc(double *, short *, int, int, double *);
void fbisinc_tops(double *, short *, float *, int , int , double *);
void read_orb(FILE *, struct PRM *, struct SAT_ORB *); 
void set_prm_defaults(struct PRM *); 
void hermite_c(double *, double *, double *, int, int, double, double *, int *); 
void set_prm_defaults(struct PRM *); 
void interpolate_SAT_orbit_slow(struct SAT_ORB *orb, double time, double *, double *, double *, int *); 
void polyfit(double *, double *, double *, int *, int *); 
int calorb_alos(struct SAT_ORB *, double **, double, double, int );


#define R 0.61803399
#define C 0.382
#define SHFT2(a, b, c)                                                                                                           \
    (a) = (b);                                                                                                                   \
    (b) = (c);
#define SHFT3(a, b, c, d)                                                                                                        \
    (a) = (b);                                                                                                                   \
    (b) = (c);                                                                                                                   \
    (c) = (d);
#define TOL 2

int npad = 8000;

dcomplex Cmuld(dcomplex x, dcomplex y) {
    dcomplex z;
    z.r = x.r * y.r - x.i * y.i;
    z.i = x.i * y.r + x.r * y.i;
    return z;
}
dcomplex Cexpd(double theta) {
    dcomplex z;
    z.r = cos(theta);
    z.i = sin(theta);
    return z;
}

int main (int argc, char **argv) {
    
    FILE *SLCfile = NULL, *ldrfile = NULL, *RMPfile = NULL;
    int i,j,k,ii,jj;
    struct PRM p1;
    double **orb_pos = NULL;
    struct SAT_ORB *orb = NULL;
    double dr, t1, t11, t2, tm;
    double ts, rng0,ras[2],cnst,pha;
    dcomplex r_i, r_i2, pshif;
    double xp[3];
    double xt[3];
    double rp[3];
    double r0, rf, a0, af;
    double fll, rdd, daa, drr, dopc;
    double dt, dtt, xs, ys, zs;
    double time[20], rng[20], d[3]; /* arrays used for polynomial refinement of min range */
    int ir, ntt = 10, nc = 3;    /* size of arrays used for polynomial refinement */
    int nrec;
    int goldop();
    int stai, endi, midi;
    struct PRM prm;
    void *API = NULL;
    struct GMT_GRID *DEM = NULL, *OUT_R = NULL, *OUT_I = NULL;
    float *real, *imag;
    int fdin,pdin;
    size_t st_size;
    short *sinn = NULL;
    float *pinn = NULL;
    char tmp1[256];

    if (argc != 3) {
        fprintf(stderr, "%s\n", USAGE);
        exit(-1);
    } 
    

    if ((API = GMT_Create_Session(argv[0], 0U, 0U, NULL)) == NULL)
        return EXIT_FAILURE;

    get_prm(&prm,argv[1]);

    if ((fdin = open(prm.SLC_file, O_RDONLY)) < 0)
        die("can't open %s for reading", prm.SLC_file);
    st_size = (size_t)4 * (size_t)prm.num_rng_bins * (size_t)prm.num_patches*prm.num_valid_az;
    if ((sinn = mmap(0, st_size, PROT_READ, MAP_SHARED, fdin, 0)) == MAP_FAILED)
        die("mmap error for input", " ");
    // read in ramp file for sampling purpose
    if (prm.SC_identity == 10) {
        strcpy(tmp1,argv[1]);
        tmp1[strlen(tmp1)-4] = '\0';
        strcat(tmp1,".RMP");
        if ((pdin = open(tmp1, O_RDONLY)) < 0)
            die("can't open %s for reading", tmp1);
        st_size = (size_t)4 * (size_t)prm.num_rng_bins * (size_t)prm.num_patches*prm.num_valid_az;
        if ((pinn = mmap(0, st_size, PROT_READ, MAP_SHARED, pdin, 0)) == MAP_FAILED)
            die("mmap error for input", " ");
    }

    if ((ldrfile = fopen(prm.led_file, "r")) == NULL)
        die("Can't open LEDfile", p1.led_file);
    if ((DEM = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, argv[2], NULL)) == NULL)
        die("cannot open DEM", argv[2]);
    if (GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, argv[2], DEM) == NULL)
        return EXIT_FAILURE;

    real = (float *)malloc(DEM->header->n_columns*DEM->header->n_rows*sizeof(float));
    imag = (float *)malloc(DEM->header->n_columns*DEM->header->n_rows*sizeof(float));

    dr = 0.5 * SOL / prm.fs;
    r0 = -10.;
    rf = prm.num_rng_bins + 10.;
    a0 = -20.;
    af = prm.num_patches * prm.num_valid_az + 20.;

    /* compute the flattening */

    fll = (prm.ra - prm.rc) / prm.ra;

    /* compute the start time, stop time and increment */

    t1 = 86400. * prm.clock_start + (prm.nrows - prm.num_valid_az) / (2. * prm.prf);
    t2 = t1 + prm.num_patches * prm.num_valid_az / prm.prf;

    /* sample the orbit only every 2th point or about 8 m along track */
    /* if this is S1A which has a low PRF sample 2 times more often */

    ts = 2. / prm.prf;
    if (prm.prf < 600.) {
        ts = 2. / (2. * prm.prf);
        npad = 20000;
    }
    nrec = (int)((t2 - t1) / ts);

    /* allocate memory for an array of pointers*/
    orb = (struct SAT_ORB *)malloc(sizeof(struct SAT_ORB));
    /* for each pointer, allocate storage for an array of floats  */
    read_orb(ldrfile, &prm, orb);
    orb_pos = malloc(4 * sizeof(double *));
    for (j = 0; j < 4; j++) {
         orb_pos[j] = malloc((nrec + 2 * npad) * sizeof(double));
    }

    /* read in the postion of the orbit */
    (void)calorb_alos(orb, orb_pos, ts, t1, nrec);

    /* go through every point from the DEM */
    for (ii=0;ii<DEM->header->n_rows;ii++) {
        for (jj=0;jj<DEM->header->n_columns;jj++) {
            rp[0] = DEM->header->wesn[3]-ii*DEM->header->inc[1];    // latitude
            rp[1] = DEM->header->wesn[0]+jj*DEM->header->inc[0];    // longitude
            rp[2] = DEM->data[ii*DEM->header->n_columns+jj];              // elevation
            
            plh2xyz(rp, xp, prm.ra, fll);
            if (rp[1] > 180.)
                rp[1] = rp[1] - 360.;
            xt[0] = -1.0;

            /* compute the topography due to the difference between the local radius and center radius */
            rp[2] = sqrt(xp[0] * xp[0] + xp[1] * xp[1] + xp[2] * xp[2]) - prm.RE;
            /* minimum for each point */
            stai = 0;
            endi = nrec + npad * 2 - 1;
            midi = (stai + (endi - stai) * C); 
            (void)goldop(ts, t1, orb_pos, stai, endi, midi, xp[0], xp[1], xp[2], &rng0, &tm);
        
            /* refine this minimum range and azimuth with a polynomial fit */
            dt = 1. / ntt; /* make the polynomial 1 second long */
            for (k = 0; k < ntt; k++) {
                time[k] = dt * (k - ntt / 2 + .5);
                t11 = tm + time[k];
                interpolate_SAT_orbit_slow(orb, t11, &xs, &ys, &zs, &ir);
                rng[k] = sqrt((xp[0] - xs) * (xp[0] - xs) + (xp[1] - ys) * (xp[1] - ys) + (xp[2] - zs) * (xp[2] - zs)) - rng0;
            } 

            /* fit a second order polynomial to the range versus time function and update the tm and rng0 */
            polyfit(time, rng, d, &ntt, &nc);
            dtt = -d[1] / (2. * d[2]);
            tm = tm + dtt;
            interpolate_SAT_orbit_slow(orb, tm, &xs, &ys, &zs, &ir);
            rng0 = sqrt((xp[0] - xs) * (xp[0] - xs) + (xp[1] - ys) * (xp[1] - ys) + (xp[2] - zs) * (xp[2] - zs));

            /* compute the range and azimuth in pixel space */
            xt[0] = rng0;
            xt[1] = tm;
            xt[0] = (xt[0] - prm.near_range) / dr - (prm.rshift + prm.sub_int_r) + prm.chirp_ext;
            xt[1] = prm.prf * (xt[1] - t1) - (prm.ashift + prm.sub_int_a);

            /* For Envisat correct for biases based on Pinon reflector analysis */
            if (prm.SC_identity == 4) {
                xt[0] = xt[0] + 8.4;
                xt[1] = xt[1] + 4;
            }


            /* compute the azimuth and range correction if the Doppler is not zero */
            if (prm.fd1 != 0.) {
                dopc = prm.fd1 + prm.fdd1 * (prm.near_range + dr * prm.num_rng_bins / 2.);
                rdd = (prm.vel * prm.vel) / rng0;
                daa = -0.5 * (prm.lambda * dopc) / rdd;
                drr = 0.5 * rdd * daa * daa / dr;
                daa = prm.prf * daa;
                xt[0] = xt[0] + drr;
                xt[1] = xt[1] + daa;
            }
            
            /* write into the output grids. */ 
            ras[0] = xt[0];
            ras[1] = xt[1];

            r_i.r = 0.; r_i.i = 0.;
            if (prm.SC_identity == 10) {
                fbisinc_tops(ras,sinn,pinn,prm.num_patches*prm.num_valid_az,prm.num_rng_bins,&r_i.r);
                //fbisinc(ras,sinn,prm.num_patches*prm.num_valid_az,prm.num_rng_bins,&r_i.r);
            }
            else {
                fbisinc(ras,sinn,prm.num_patches*prm.num_valid_az,prm.num_rng_bins,&r_i.r);
            }
            
            //cnst = -4.0 * PI / prm.lambda;
            cnst = 4.0 * PI / prm.lambda;
            pha = cnst*xt[0]*dr;
            pshif = Cexpd(pha);
            r_i2 = Cmuld(r_i,pshif);
            real[ii*DEM->header->n_columns+jj] = r_i2.r;
            imag[ii*DEM->header->n_columns+jj] = r_i2.i;
        }
    }

    if (OUT_R == NULL && (OUT_R = GMT_Duplicate_Data(API, GMT_IS_GRID, GMT_DUPLICATE_DATA, DEM)) == NULL)
        die("error creating output grid", "real.grd");
    if (OUT_I == NULL && (OUT_I = GMT_Duplicate_Data(API, GMT_IS_GRID, GMT_DUPLICATE_DATA, DEM)) == NULL)
        die("error creating output grid", "imag.grd");
    
    if (GMT_Set_Comment(API, GMT_IS_GRID, GMT_COMMENT_IS_REMARK, "Geocoded SLC Real", OUT_R))
        die("could not set title", "");
    if (GMT_Set_Comment(API, GMT_IS_GRID, GMT_COMMENT_IS_REMARK, "Geocoded SLC Imag", OUT_I))
        die("could not set title", "");
    strcpy(tmp1, "");
    for (i = 0; i < argc; i++) {
        strcat(tmp1, argv[i]);
        strcat(tmp1, " ");
    }
    strcpy(OUT_R->header->command, "");
    strcpy(OUT_I->header->command, "");
    if (GMT_Set_Comment(API, GMT_IS_GRID, GMT_COMMENT_IS_COMMAND, tmp1, OUT_R))
        die("GMT[ERROR]: could not set title", "");
    if (GMT_Set_Comment(API, GMT_IS_GRID, GMT_COMMENT_IS_COMMAND, tmp1, OUT_I))
        die("GMT[ERROR]: could not set title", "");

    strcpy(OUT_R->header->title, "");
    strcpy(OUT_I->header->title, "");
    strcpy(tmp1, "Produced by geocode_slc");
    if (GMT_Set_Comment(API, GMT_IS_GRID, GMT_COMMENT_IS_TITLE, tmp1, OUT_R))
        die("GMT[ERROR]: could not set title", "");
    if (GMT_Set_Comment(API, GMT_IS_GRID, GMT_COMMENT_IS_TITLE, tmp1, OUT_I))
        die("GMT[ERROR]: could not set title", "");

    OUT_R->data = real;
    OUT_I->data = imag;


    if (GMT_Write_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, "real.grd", OUT_R))
        die("Failed to write output grid ", "real.grd");
    if (GMT_Write_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, "imag.grd", OUT_I))
        die("Failed to write output grid ", "imag.grd");
    
    for (j = 0; j < 4; j++) {
        free(orb_pos[j]);
    } 
    free(orb_pos);
    free(orb);
    fclose(SLCfile);
    fclose(ldrfile);
    if (prm.SC_identity == 10)
        fclose(RMPfile);

    if (GMT_Destroy_Data(API, &OUT_R))
        die("error freeing data ", "real.grd");
    if (GMT_Destroy_Data(API, &OUT_I))
        die("error freeing data ", "imag.grd");

    if (GMT_Destroy_Session(API))
        return EXIT_FAILURE; /* Remove the GMT machinery */

    return(1);
}

#define PI 3.1415926535897932

double sinc_kernel(double x) { 
    double arg, f;

    arg = fabs(PI * x);
    if (arg > 0.) {
        f = sin(arg) / arg;
    }
    else {
        f = 1.;
    }
    return (f);
}


void sinc_one(double *rdata, double *idata, double x, double y, double *cz) {
    int i, j, ij, ns2 = NS / 2 - 1; 
    double wx[NS], wy[NS];
    double arg, w, wsum, rsum, isum;

    for (i = 0; i < NS; i++) {
        arg = fabs(x + ns2 - i);
        wx[i] = sinc_kernel(arg);
        arg = fabs(y + ns2 - i);
        wy[i] = sinc_kernel(arg);
    }    

    rsum = isum = wsum = 0.0; 
    ij = 0; 
    for (j = 0; j < NS; j++) {
        for (i = 0; i < NS; i++) {
            w = wx[i] * wy[j];
            rsum += rdata[ij + i] * w; 
            isum += idata[ij + i] * w; 
            wsum += w;
        }    
        ij += NS;
    }    
    if (wsum <= 0.0) 
        printf(" error wsum is zero \n");
    cz[0] = rsum / wsum;
    cz[1] = isum / wsum;
}

void fbisinc(double *ras, short *s_in, int ydims, int xdims, double *sout) {
    double dr, da, ns2 = NS / 2 - 1;
    double rdata[NS * NS], idata[NS * NS], cz[2];
    int i, j, k, kk;
    int i0, j0;
    int nclip;

    /* compute the residual offsets */
    nclip = 0;
    j0 = (int)floor(ras[0]);
    i0 = (int)floor(ras[1]);
    dr = ras[0] - (double)j0;
    da = ras[1] - (double)i0;
    if (dr < 0. || dr > 1. || da < 0. || da > 1)
        fprintf(stderr, " dr or da out of bounds %f %f \n", dr, da);

    /* make sure all 4 corners are within the bounds of the aligned array */

    if ((i0 - ns2) < 0 || (i0 + ns2 + 1) >= ydims || (j0 - ns2) < 0 || (j0 + ns2 + 1) >= xdims) {
        sout[0] = NAN;
        sout[1] = NAN;
    }
    else {

        /* safe to do the interpolation */

        for (i = 0; i < NS; i++) {
            for (j = 0; j < NS; j++) {
                k = i * NS + j;
                kk = xdims * (i0 - ns2 + i) * 2 + (j0 - ns2 + j) * 2;
                rdata[k] = (double)s_in[kk];
                idata[k] = (double)s_in[kk+1];
            }
        }

        /* interpolate the real and imaginary data */

        sinc_one(rdata, idata, dr, da, cz);

        if ((int)fabs(cz[0]) > I2MAX)
            nclip = nclip + 1;
        sout[0] = cz[0];
        if ((int)fabs(cz[1]) > I2MAX)
            nclip = nclip + 1;
        sout[1] = cz[1];
    }    
    // if(nclip > 0) fprintf(stderr," %d integers were clipped \n",nclip);
}

void sample_one_p(double *pdata, double x, double y, double *p) {
    int i, j, ij, ns2 = NS / 2 - 1;
    double wx[NS], wy[NS];
    double arg, w, wsum, psum;
    
    for (i = 0; i < NS; i++) {
        arg = fabs(x + ns2 - i);
        wx[i] = 1/arg;
        arg = fabs(y + ns2 - i);
        wy[i] = 1/arg;
    }
    
    psum = wsum = 0.0;
    ij = 0;
    for (j = 0; j < NS; j++) {
        for (i = 0; i < NS; i++) {
            w = wx[i] * wy[j];
            psum += pdata[ij + i] * w;
            wsum += w;
        }
        ij += NS;
    }
    if (wsum <= 0.0)
        printf(" error wsum is zero \n");
    *p = psum / wsum;
}

void fbisinc_tops(double *ras, short *s_in, float *p_in, int ydims, int xdims, double *sout) {
    double dr, da, ns2 = NS / 2 - 1;
    double rdata[NS * NS], idata[NS * NS], pdata[NS * NS], cz[2], pp;
    dcomplex dd,dd2,dd3;
    int i, j, k, kk;
    int i0, j0;
    int nclip;
    
    /* compute the residual offsets */
    nclip = 0;
    j0 = (int)floor(ras[0]);
    i0 = (int)floor(ras[1]);
    dr = ras[0] - (double)j0;
    da = ras[1] - (double)i0;
    if (dr < 0. || dr > 1. || da < 0. || da > 1)
        fprintf(stderr, " dr or da out of bounds %f %f \n", dr, da);
    
    /* make sure all 4 corners are within the bounds of the aligned array */
    
    if ((i0 - ns2) < 0 || (i0 + ns2 + 1) >= ydims || (j0 - ns2) < 0 || (j0 + ns2 + 1) >= xdims) {
        sout[0] = NAN;
        sout[1] = NAN;
    }
    else {
        
        /* safe to do the interpolation */
        
        for (i = 0; i < NS; i++) {
            for (j = 0; j < NS; j++) {
                k = i * NS + j;
                kk = xdims * (i0 - ns2 + i)  + (j0 - ns2 + j);
                pdata[k] = (double)p_in[kk];
                dd = Cexpd((double)pdata[k]);
                dd2.r = (double)s_in[kk*2];
                dd2.i = (double)s_in[kk*2+1];
                dd3 = Cmuld(dd2,dd);
                rdata[k] = dd3.r;
                idata[k] = dd3.i;
                
                //rdata[k] = (double)s_in[kk*2];
                //idata[k] = (double)s_in[kk*2+1];
            }
        }
        
        /* interpolate the real and imaginary data */
        
        sinc_one(rdata, idata, dr, da, cz);

        sample_one_p(pdata, dr, da, &pp);
        dd = Cexpd(pp);
        dd.i = -dd.i;
        dd2.r = cz[0];
        dd2.i = cz[1];
        
        dd3 = Cmuld(dd2,dd);
        cz[0] = dd3.r;
        cz[1] = dd3.i;

        if ((int)fabs(cz[0]) > I2MAX)
            nclip = nclip + 1;
        sout[0] = cz[0];
        if ((int)fabs(cz[1]) > I2MAX)
            nclip = nclip + 1;
        sout[1] = cz[1];
    }
    // if(nclip > 0) fprintf(stderr," %d integers were clipped \n",nclip);
}

/*    subfunctions    */

int goldop(double ts, double t1, double **orb_pos, int ax, int bx, int cx, double xpx, double xpy, double xpz, double *rng,
           double *tm) {

        /* use golden section search to find the minimum range between the target and
         * the orbit */
        /* xpx, xpy, xpz is the position of the target in cartesian coordinate */
        /* ax is stai; bx is endi; cx is midi it's easy to tangle */

        double f1, f2;
        int x0, x1, x2, x3;
        int xmin;
        double dist();

        x0 = ax;
        x3 = bx;
        //      if (fabs(bx-cx) > fabs(cx-ax)) {
        if (abs(bx - cx) > abs(cx - ax)) {
                x1 = cx;
                x2 = cx + (int)fabs((C * (bx - cx)));
        }
        else {
                x2 = cx;
                x1 = cx - (int)fabs((C * (cx - ax))); /* make x0 to x1 the smaller segment */
        }

        f1 = dist(xpx, xpy, xpz, x1, orb_pos);
        f2 = dist(xpx, xpy, xpz, x2, orb_pos);

        while ((x3 - x0) > TOL && (x2 != x1)) {
                if (f2 < f1) {
                        SHFT3(x0, x1, x2, (int)(R * x3 + C * x1));
                        SHFT2(f1, f2, dist(xpx, xpy, xpz, x2, orb_pos));
                }
                else {
                        SHFT3(x3, x2, x1, (int)(R * x0 + C * x2));
                        SHFT2(f2, f1, dist(xpx, xpy, xpz, x1, orb_pos));
                }
        }

        if (f1 < f2) {
        if (x1 <= bx && x1 >= ax) {
                    xmin = x1;
        }
        else{
            xmin = abs(x1-bx) > abs(x1-ax) ? ax : bx;
        }
                *tm = orb_pos[0][x1];
                *rng = f1;
        }
        else {
        if (x2 <= bx && x2 >= ax) {
                    xmin = x2;
        }
        else {
            xmin = abs(x2-bx) > abs(x2-ax) ? ax : bx;
        }
                *tm = orb_pos[0][x2];
                *rng = f2;
        }

        return (xmin);
}

double dist(double x, double y, double z, int n, double **orb_pos) {

        double d, dx, dy, dz;

        dx = x - orb_pos[1][n];
        dy = y - orb_pos[2][n];
        dz = z - orb_pos[3][n];
        d = sqrt(dx * dx + dy * dy + dz * dz);

        return (d);
}

int calorb_alos(struct SAT_ORB *orb, double **orb_pos, double ts, double t1, int nrec)
/* function to calculate every position in the orbit   */

{
        int i, k, nval;
        // int     npad = 8000;   /* number of buffer points to add before and after
        // the acquisition */
        int ir;            /* return code: 0 = ok; 1 = interp not in center; 2 = time out of
                              range */
        double xs, ys, zs; /* position at time */
        double *pt, *px, *py, *pz, *pvx, *pvy, *pvz;
        double pt0;
        double time;

        px = (double *)malloc(orb->nd * sizeof(double));
        py = (double *)malloc(orb->nd * sizeof(double));
        pz = (double *)malloc(orb->nd * sizeof(double));
        pvx = (double *)malloc(orb->nd * sizeof(double));
        pvy = (double *)malloc(orb->nd * sizeof(double));
        pvz = (double *)malloc(orb->nd * sizeof(double));
        pt = (double *)malloc(orb->nd * sizeof(double));

        pt0 = 86400. * orb->id + orb->sec;
        for (k = 0; k < orb->nd; k++) {
                pt[k] = pt0 + k * orb->dsec;
                px[k] = orb->points[k].px;
                py[k] = orb->points[k].py;
                pz[k] = orb->points[k].pz;
                pvx[k] = orb->points[k].vx;
                pvy[k] = orb->points[k].vy;
                pvz[k] = orb->points[k].vz;
        }

        nval = 6;

        /* loop to get orbit position of every point and store them into orb_pos */
        for (i = 0; i < nrec + npad * 2; i++) {
                time = t1 - npad * ts + i * ts;
                orb_pos[0][i] = time;

                hermite_c(pt, px, pvx, orb->nd, nval, time, &xs, &ir);
                hermite_c(pt, py, pvy, orb->nd, nval, time, &ys, &ir);
                hermite_c(pt, pz, pvz, orb->nd, nval, time, &zs, &ir);

                orb_pos[1][i] = xs;
                orb_pos[2][i] = ys;
                orb_pos[3][i] = zs;
        }

        free((double *)px);
        free((double *)py);
        free((double *)pz);
        free((double *)pt);
        free((double *)pvx);
        free((double *)pvy);
        free((double *)pvz);

        return orb->nd;
}

