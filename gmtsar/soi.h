/************************************************************************
 * soi.h is the include file for the esarp SAR processor.		*
 ************************************************************************/
/************************************************************************
 * Creator: Evelyn J. Price	(Scripps Institution of Oceanography)	*
 * Date   : 11/18/96							*
 ************************************************************************/
/************************************************************************
 * Modification History							*
 *									*
 * Date									*
 *									*
 *  4/23/97- 	added parameters for orbit calculations: x_target,      *
 *		y_target,z_target,baseline,alpha,sc_identity,		*
 *		ref_identity,SC_clock_start,SC_clock_stop,              *
 *		clock_start,clock_stop   				*
 *		-DTS							*
 *									*
 * 4/23/97-	added parameters: rec_start, rec_stop			*
 *		-EJP							*
 *									*
 * 8/28/97-	added parameters baseline_start baseline_end		*
 *		alpha_start alpha_end					*
 *									*
 * 9/12/97	added clipi2 function to clip to short int		*
 *									*
 * 4/26/06	added nrows, num_lines					*
 ************************************************************************/
#ifndef SOI_H
#define SOI_H
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define SOL 299792456.0
#define PI 3.1415926535897932
#define PI2 6.2831853071795864
#define I2MAX1 32767.
#define I2SCALE 4.e6
#define TRUE 1
#define FALSE 0
#define RW 0666
#define MULT_FACT 1000.0
#define sgn(A) ((A) >= 0.0 ? 1.0 : -1.0)
#define clipi22(A) (((A) > I2MAX1) ? I2MAX1 : (((A) < -I2MAX1) ? -I2MAX1 : A))
#include "sfd_complex.h"

char *input_file;
char *led_file;
char *out_amp_file;
char *out_data_file;
char *deskew;
char *iqflip;
char *off_vid;
char *srm;
char *ref_file;
char *orbdir;
char *lookdir;

int debug_flag;
int bytes_per_line;
int good_bytes;
int first_line;
int num_patches;
int first_sample;
int num_valid_az;
int st_rng_bin;
int num_rng_bins;
int nextend;
int nlooks;
int xshift;
int yshift;
int fdc_ystrt;
int fdc_strt;

/*New parameters 4/23/97 -EJP */
int rec_start;
int rec_stop;
/* End new parameters 4/23/97 -EJP */

/* New parameters 4/23/97 -DTS */
int SC_identity;       /* (1)-ERS1 (2)-ERS2 (3)-Radarsat (4)-Envisat (5)-ALOS
                          (6)-Envisat_SLC  (7)-TSX (8)-CSK (9)-RS2 (10)-S1A*/
int ref_identity;      /* (1)-ERS1 (2)-ERS2 (3)-Radarsat (4)-Envisat (5)-ALOS
                          (6)-Envisat_SLC  (7)-TSX (8)-CSK (9)-RS2 (10)-S1A*/
double SC_clock_start; /* YYDDD.DDDD */
double SC_clock_stop;  /* YYDDD.DDDD */
double icu_start;      /* onboard clock counter */
double clock_start;    /* DDD.DDDDDDDD  clock without year has more precision */
double clock_stop;     /* DDD.DDDDDDDD  clock without year has more precision */
/* End new parameters 4/23/97 -DTS */

double caltone;
double RE;   /* Local Earth radius */
double raa;  /* ellipsoid semi-major axis - added by RJM */
double rcc;  /* ellipsoid semi-minor axis - added by RJM */
double vel1; /* Equivalent SC velocity */
double ht1;  /* (SC_radius - RE) center of frame*/
double ht0;  /* (SC_radius - RE) start of frame */
double htf;  /* (SC_radius - RE) end of frame */
double near_range;
double far_range;
double prf1;
double xmi1;
double xmq1;
double az_res;
double fs;
double slope;
double pulsedur;
double lambda;
double rhww;
double pctbw;
double pctbwaz;
double fd1;
double fdd1;
double fddd1;
double sub_int_r;
double sub_int_a;
double stretch_r;
double stretch_a;
double a_stretch_r;
double a_stretch_a;

/* New parameters 8/28/97 -DTS */
double baseline_start;
double baseline_center;
double baseline_end;
double alpha_start;
double alpha_center;
double alpha_end;
/* New parameters 9/25/18 -EXU */
double B_offset_start;
double B_offset_center;
double B_offset_end;
/* End new parameters 8/28/97 -DTS */
double bparaa; /* parallel baseline - added by RJM */
double bperpp; /* perpendicular baseline - added by RJM */

/* New parameters 4/26/06 */
int nrows;
int num_lines;

/* New parameters 09/18/08 */
double TEC_start;
double TEC_end;
#endif /* SOI_H	*/
