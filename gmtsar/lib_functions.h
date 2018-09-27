/*	$Id: lib_functions.h 39 2013-04-07 00:49:34Z pwessel $	*/
/* include files to define sarleader structure */
#ifndef LIB_FUNCTIONS_H
#define LIB_FUNCTIONS_H
#include "sarleader_ALOS.h"
#include "sarleader_fdr.h"
#include "sfd_complex.h"
#include "xcorr.h"
#include "PRM.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* function prototypes 				*/
void null_sio_struct(struct PRM *);
void get_sio_struct(FILE *, struct PRM *);
void put_sio_struct(struct PRM, FILE *);
void get_string(char *, char *, char *, char *);
void get_int(char *, char *, char *, int *);
void get_double(char *, char *, char *, double *);
void ALOS_ldr_prm(struct SAR_info, struct PRM *);
int is_big_endian_(void);
int is_big_endian__(void);
void die(char *, char *);
int get_prm(struct PRM *p, char *filename);
void cross3(double *, double *, double *);
void get_seconds(struct PRM, double *, double *);
void plh2xyz(double *, double *, double, double);
void xyz2plh(double *, double *, double, double);
void find_unit_vector(double *, double *);
double find_length(double *);
double find_distance(double *, double *);
double find_distance3(double, double, double, double, double, double);
int geo2latlon(double *, double *, struct PRM);
void geoxyz(double, double, double, double *, double *);
int spline_(int *istart, int *nn, double *x, double *u, double *s, double *a);
int evals_(int *istart, double *y, int *nn, double *x, double *u, double *s, double *eval);
int find_fft_length(int n);
void aastretch(fcomplex **fdata, int ipatch, int nrows, int num_valid_az, int num_rng_bins, float coef);
void acpatch(void *API, fcomplex **data, int nrows, double delr, double fd, double fdd, double fddd);
void conv2d(float *rdat, int *ni, int *nj, float *filt, int *nif, int *njf, float *fdat, int *ic, int *jc, float *rnorm);
void do_freq_corr(void *API, struct xcorr *xc, int iloc);
void do_time_corr(struct xcorr *xc, int iloc);
double calc_time_corr(struct xcorr *xc, int ioff, int joff);
int fft_bins(int num);
void fft_interpolate_1d(void *API, struct FCOMPLEX *in, int N, struct FCOMPLEX *out, int ifactor);
void fft_interpolate_2d(void *API, struct FCOMPLEX *in, int N1, int M1, struct FCOMPLEX *out, int N, int M, int ifactor);
void print_prm_params(struct PRM p1, struct PRM p2);
void fix_prm_params(struct PRM *p, char *s);
void get_locations(struct xcorr *xc);
void get_params(FILE *fh);
void do_highres_corr(void *API, struct xcorr *xc, int iloc);
void intp_coef(int nfilter, float *xintp);
void print_params(struct xcorr *xc);
void set_defaults(struct xcorr *xc);
void parse_command_line(int na, char **a, struct xcorr *xc, int *nfiles, int *input_flag, char *USAGE);
void handle_prm(void *, char **argv, struct xcorr *xc, int nfiles);
void print_results(struct xcorr *xc, int iloc);
void print_complex(struct FCOMPLEX *a, int ny, int nx, int real_flag);
void print_float(float *a, int ny, int nx);
void print_double(double *a, int ny, int nx);
void print_int(int *a, int ny, int nx);
void radopp(double *fd, double *fdd, double *fddd, double r, double del);
void read_complex_short(FILE *f, int *d, int iy, int jx, int npx, int npy, int nx);
void read_real_float(FILE *f, int *d, int iy, int jx, int npx, int npy, int nx);
void read_data(struct xcorr xc);
void read_complex_short2float(FILE *f, float *d, int iy, int jx, int npx, int npy, int nx);
void read_optional_args(void *API, int argc, char **argv, struct PRM *tp, int *topoflag, struct PRM *mp, int *modelflag);
void read_xcorr_data(struct xcorr *xc, int iloc);
void rmpatch(fcomplex **data, int nrows, double delr, double fd, double fdd, double fddd);
void rng_cmp(void *API, int ranfft, fcomplex *data, fcomplex *ref);
void rng_ref(void *API, int ranfft, float delr, fcomplex *ref1);
void rng_filter(void *API, fcomplex *cin, int nffti, fcomplex *cout);
void shift(void *API, int ranfft, fcomplex *data, double shift);
int trans_col(void *API, int xnum, int ynum, fcomplex **data);
int read_SLC_short2float(FILE *SLCfile, char *name, short *sdata, fcomplex *cdata, int xdim, int psize, double dfact);
int read_SLC_short2double(FILE *SLCfile, char *name, short *sdata, dcomplex *cdata, int xdim, int psize, double dfact);
void handle_input(char *, struct xcorr *);
void read_params(struct xcorr *, FILE *);
void make_mask(struct xcorr *);
void do_highres(struct xcorr *, int);
void allocate_arrays(struct xcorr *);
char *trimwhitespace(char *str);

#endif /* LIB_FUNCTIONS_H */
