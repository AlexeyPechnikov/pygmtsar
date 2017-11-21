/*	$Id: xcorr.h 39 2013-04-07 00:49:34Z pwessel $	*/
#ifndef XCORR_H
#define XCORR_H
#include <stdio.h>
struct locs {
	int	x;		/* x pixel location */
	int	y;		/* y pixel location */
	int	qflag;		/* quality flag 1 = good */
	float	corr;		/* correlation, normalized time domain */
	float	yoff;		/* estimated y offset */
	float	xoff;		/* estimated x offset */
	float	xfrac;		/* estimated x offset (fraction) */
	float	yfrac;		/* estimated y offset (fraction) */
	int	m1;		/* mean value */
	int	m2;		/* mean value */
	};

struct xcorr{
	int	format;		/* type of input data [0 short complex, 1 real float, 3 real grd] */
	int	corr_flag;	/* type of correlation flag 0 = standard; 1 = Gatelli; 2 = fft */
	int	offset_flag;	/* set offset to zero (ignore prm)  */
	int	interp_flag;	/* interpolation flag 1 = yes 0 = no */
	int	interp_factor;	/* interpolation factor (power of 2) */
	int	nx_corr;	/* size of correlation window */
	int	ny_corr;	/* size of correlation window */
	int	xsearch;	/* size of x search offset */
	int	ysearch;	/* size of y search offset */
	int	m_nx;		/* x size of master file */
	int	m_ny;		/* y size of master file */
	int	s_nx;		/* x size of slave file */
	int	s_ny;		/* y size of slave file */
	int	nxc;		/* x size of correlation */
	int	nyc;		/* y size of correlation */
	int	npx;		/* x size of patch */
	int	npy;		/* y size of patch */
	int	nxl;		/* x number of locations */
	int	nyl;		/* y number of locations */
	int	x_offset;	/* intial starting offset in x */
	int	y_offset;	/* intial starting offset in y */
	int	nlocs;		/* number of locations */
	int	x_inc;		/* x distance between locations */
	int	y_inc;		/* y distance between locations */
	int	ri;		/* range interpolation factor (must be power of two) */
	short	*mask;		/* mask file (short integer) */
	int	*i1;		/* data matrix 1 (integer) */
	int	*i2;		/* data matrix 2 (integer) */
	int	n2x;		/* size of interpolation */
	int	n2y;		/* size of interpolation */
	struct FCOMPLEX	*d1;	/* data 1 (amplitude in real, imag = 0)*/
	struct FCOMPLEX	*d2;	/* data 2 (amplitude in real, imag = 0)*/
	struct FCOMPLEX	*c1;	/* subset data patch 1 (complex float) */
	struct FCOMPLEX	*c2;	/* subset data patch 2 (complex float) */
	struct FCOMPLEX	*c3;	/* c1 * c2 (complex float) */
	struct FCOMPLEX	*ritmp;	/* tmp array for range interpoation */
	double	*corr;		/* correlation (real double) */
	double  astretcha;	/* azimuth stretch parameter estimated from (prf2-prf1)/prf1 */
	struct FCOMPLEX	*md;	/* interpolation file */
	struct FCOMPLEX	*md_exp;	/* interpolation file */
	struct FCOMPLEX	*sd;	/* interpolation file */
	struct FCOMPLEX	*sd_exp;	/* interpolation file */
	struct FCOMPLEX	*cd_exp;	/* interpolation file */
	double	*interp_corr;	/* interpolation file */
	FILE	*data1;		/* data file 1 */
	FILE	*data2;		/* data file 2 */
	FILE	*param;		/* input parameters file */
	FILE	*file;		/* output file (offsets) */
	char	param_name[128];
	char	data1_name[128];
	char	data2_name[128];
	char	filename[128];
	struct  locs *loc;
    struct GMT_GRID *D1;
    struct GMT_GRID *D2;
	};
#endif /* XCORR_H */
