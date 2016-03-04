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
*		clock_start,clock_stop	         			*
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
	    
#include <stdio.h>
#include <string.h> 
#include <math.h>
#define SOL 299792456.0
#define PI 3.1415926535897932
#define PI2 6.2831853071795864
#define I2MAX 32767.
#define I2SCALE 4.e6
#define TRUE 1
#define FALSE 0
#define RW 0666
#define MULT_FACT 1000.0 
#define sgn(A) ((A) >= 0.0 ? 1.0 : -1.0)
#define clipi2(A) ( ((A) > I2MAX) ? I2MAX : (((A) < -I2MAX) ? -I2MAX : A) )

typedef struct SCOMPLEX {short r,i;} scomplex;
typedef struct FCOMPLEX {float r,i;} fcomplex;
typedef struct DCOMPLEX {double r,i;} dcomplex;

char *input_file;
char *out_amp_file;
char *out_data_file;
char *deskew;
char *iqflip;
char *off_vid;
char *srm;
char *ref_file;
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
int SC_identity;	/* (1)-ERS1 (2)-ERS2 (3)-Radarsat (4)-Envisat (5)-ALOS */
int ref_identity;	/* (1)-ERS1 (2)-ERS2 (3)-Radarsat (4)-Envisat (5)-ALOS */
double SC_clock_start;	/* YYDDD.DDDD */
double SC_clock_stop;	/* YYDDD.DDDD */
double icu_start;       /* onboard clock counter */
double clock_start;     /* DDD.DDDDDDDDDD same as SC_clock_start but no date and more precise */
double clock_stop;      /* DDD.DDDDDDDDDD same as SC_clock_start but no date and more precise */
/* End new parameters 4/23/97 -DTS */

double caltone;
double RE;		/* Local Earth radius */
double vel1;		/* Equivalent SC velocity */
double ht1;		/* (SC_radius - RE) */
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
double baseline_end;
double alpha_start;
double alpha_end;
/* End new parameters 8/28/97 -DTS */

/* New parameters 4/26/06 */
int nrows;
int num_lines;
