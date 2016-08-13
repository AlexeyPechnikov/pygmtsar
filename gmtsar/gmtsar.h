/*	$Id: gmtsar.h 109 2015-01-19 23:01:24Z sandwell $	*/
/* taken from soi.h */
#ifndef GMTSAR_H
#define GMTSAR_H
#include "gmt.h"
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "sfd_complex.h"
#include "lib_functions.h"
#ifndef M_PI
#define M_PI          3.14159265358979323846
#endif
#ifndef SOI_H
#define SOL 299792456.0
#define PI 3.1415926535897932
#define PI2 6.2831853071795864
#define I2MAX 32767.0
#define I2SCALE 4.e6
#define DFACT 2.5e-07 
#define NS  4  /*number of samples in the sinc function of 2D interpolation     */
#define TRUE 1
#define FALSE 0
#define RW 0666
#define MULT_FACT 1000.0
#define sgn(A) ((A) >= 0.0 ? 1.0 : -1.0)
#define clipi2(A) ( ((A) > I2MAX) ? I2MAX : (((A) < -I2MAX) ? -I2MAX : A) )
#endif	/* SOI_H */

#define nint(x) (int)rint(x)
#define ERS1 1
#define ERS2 2
#define RSAT 3
#define ENVS 4
#define ALOS 5

#define EXIT_FLAG 1
#define paka(p) {perror((p)); exit(EXIT_FLAG);}
#define MALLOC(p,s) if (((p) = malloc(s)) == NULL) {paka("error: malloc()  ");}

#define NULL_DATA 15
#define NULL_INT -99999
#define NULL_DOUBLE -99999.9999
#define NULL_CHAR  "XXXXXXXX"

#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))	/* min and max value macros */
#endif
#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif

#include "siocomplex.h"
#include "llt2rat.h"
#include "PRM.h"

/* global variables */
int	verbose;	/* controls minimal level of output 	*/ 
int	debug; 		/* more output 				*/
int	swap; 		/* whether to swap bytes 		*/
int	quad_pol; 	/* quad polarization data 		*/
int	force_slope; 	/* whether to force the slope 		*/
int	dopp;		/* whether to calculate doppler 	*/
int	roi_flag;	/* whether to write roi.in 		*/
int	sio_flag;	/* whether to write PRM file 		*/
int	nodata;
int	quiet_flag;
double	forced_slope;	/* value to set chirp_slope to		*/
int     SAR_mode;       /* 0 => high-res                        */
			/* 1 => wide obs                        */
			/* 2 => polarimetry                     */
			/* from ALOS Product Format 3-2         */
#endif /* GMTSAR_H */
