/*	$Id$	*/
/***************************************************************************
 * resamp resamples a slave image to match the geometry of a master image. * 
 **************************************************************************/

/***************************************************************************
 * Creator:  David T. Sandwell                                             *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  03/21/13                                                      *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * DATE                                                                    *
 * 11/18/96     Code largely based on phasediff.                           *
 * 01/11/14     Code modified to use mmap() instead of fread()             *
 * 01/06/15     Code modified to use only integer rshift, ashift for       *
 *              nearest (imode = 1) interpolation.                         *
 * 04/28/16 EXU Modified 4 resampling subroutine to shift pointer before   *
 *              interpolation, so that resamp won't fail on files larger   *
 *              than 4GB.                                                  *
 ***************************************************************************/

#include "gmtsar.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <fcntl.h>

char    *USAGE = "\nUsage: "
"resamp master.PRM slave.PRM new_slave.PRM new_slave.SLC intrp \n"
"   master.PRM       - PRM for master imagea \n"
"   slave.PRM        - PRM for slave image \n"
"   new_slave.PRM    - PRM for aligned slave image \n"
"   new_slave.SLC    - SLC for aligned slave image \n"

"   intrp            - interpolation method: 1-nearest; 2-bilinear; 3-biquadratic; 4-bisinc \n \n";

void print_prm_params(struct PRM, struct PRM);
void fix_prm_params(struct PRM *, char *);
void ram2ras(struct PRM , double *, double *);
void nearest(double *, short *, int , int , short *);
void bilinear(double *, short *, int , int , short *);
void bicubic(double *, short *, int , int , short *);
void bicubic_one(double *, double *, double , double , double *);
void bisinc (double *, short *, int , int , short *);
void sinc_one(double *, double *, double , double , double *);

int main (int argc, char **argv)
{
int	ii, jj;
int	debug, intrp;
int	xdimm, ydimm;		/* size of master SLC file */
int	xdims, ydims;		/* size of slave SLC file */
short   *sinn = NULL, *sout = NULL;		/* pointer to input (whole array) and output (row) files.*/
double  ram[2], ras[2] ;	/* range and azimuth locations for master and slave images */
FILE	*SLC_file2 = NULL, *prmout = NULL;
int	fdin;
double  sv_pr[6];
size_t  st_size;

struct PRM pm, ps;	

	debug = 0;
        intrp = 2;

	if (argc < 6) die (USAGE,"");

	/* read prm file into two pointers */
	get_prm(&pm, argv[1]);
	get_prm(&ps, argv[2]);
	intrp = atoi(argv[5]);

	if (debug) print_prm_params(pm, ps);

	/* set width and length of the master and slave images */
	xdimm = pm.num_rng_bins;
	ydimm = pm.num_patches * pm.num_valid_az;
	xdims = ps.num_rng_bins;
	ydims = ps.num_patches * ps.num_valid_az;

        /* force integer interpolation if this is nearest neighbor, needed for TOPS */
        if (intrp == 1) {
          sv_pr[0] = ps.sub_int_r;
	  ps.sub_int_r = 0.;
          sv_pr[1] = ps.stretch_r;
          ps.stretch_r = 0.;
          sv_pr[2] = ps.a_stretch_r;
          ps.a_stretch_r = 0.;
          sv_pr[3] = ps.sub_int_a;
          ps.sub_int_a = 0.;
          sv_pr[4] = ps.stretch_a;
          ps.stretch_a = 0.;
          sv_pr[5] = ps.a_stretch_a;
          ps.a_stretch_a = 0.;
        }
	
	/* allocate memory for one row of the slave image */
        if((sout = (short *) malloc(2 * xdimm * sizeof(short))) == NULL){
          fprintf(stderr,"Sorry, couldn't allocate memory for output indata.\n");
          exit(-1);
	}

	/* open the input file, determine its length and mmap the input file */
	if ((fdin = open(ps.SLC_file, O_RDONLY)) < 0)
	  die ("can't open %s for reading", ps.SLC_file);

        st_size = (size_t)4*(size_t)xdims*(size_t)ydims;
        
        /* mmap the file  */
	if ((sinn = mmap (0, st_size, PROT_READ, MAP_SHARED,fdin, 0)) == MAP_FAILED) 
	  die ("mmap error for input"," ");

	/* open the slave slc file for writing and write one row at a time */
 	if ((SLC_file2 = fopen(argv[4],"w")) == NULL) die("Can't open SLCfile for output",argv[4]);
	for(ii=0; ii<ydimm; ii++) {
	   for(jj=0; jj<xdimm; jj++) {

	/* convert master ra to slave ra */
		ram[0] = jj;
		ram[1] = ii; 
		ram2ras(ps,ram,ras);
        
        /*  do nearest, bilinear, bicubic, or sinc interpolation */
	
		if(intrp == 1) {
		nearest  (ras, sinn, ydims, xdims, &sout[2*jj]);
		}
		else if(intrp == 2) {
		bilinear (ras, sinn, ydims, xdims, &sout[2*jj]);
		}
		else if(intrp == 3) {
		bicubic  (ras, sinn, ydims, xdims, &sout[2*jj]);
		}
		else if(intrp == 4) {
		bisinc  (ras, sinn, ydims, xdims, &sout[2*jj]);
		}
	   }
           fwrite(sout, 2*sizeof(short), xdimm, SLC_file2);
	}

        /* restore the affine parameters if this is nearest interpolation */
        if (intrp == 1) {
          ps.sub_int_r = sv_pr[0];
          ps.stretch_r = sv_pr[1];
          ps.a_stretch_r = sv_pr[2];
          ps.sub_int_a = sv_pr[3];
          ps.stretch_a = sv_pr[4];
          ps.a_stretch_a = sv_pr[5];
        }

	/* update and write the slave PRM file */
	ps.num_rng_bins = pm.num_rng_bins;
	ps.fs = pm.fs;
	ps.bytes_per_line = pm.bytes_per_line;
	ps.good_bytes = pm.good_bytes;
	ps.prf = pm.prf;
	ps.num_valid_az = pm.num_valid_az;
	ps.num_lines = pm.num_lines;
	ps.num_patches = pm.num_patches;
        ps.nrows = pm.nrows;
        if ((prmout = fopen(argv[3],"w")) == NULL) die("can't open prfile",argv[3]);
	put_sio_struct(ps, prmout);
        fclose(prmout);
	if (munmap(sinn, st_size) == -1) die ("mmap error unmapping file"," ");
	close(fdin);
	fclose(SLC_file2);

	return(EXIT_SUCCESS);
}

/************************************************************************
* bi-cubic interpolation algorithm modified from GMT                    *
************************************************************************/
/************************************************************************
* Creator: David Sandwell       (Scripps Institution of Oceanography)   *
* Date   : 03/22/13                                                     *
************************************************************************/
#include <math.h>
#include <stdio.h>

double cubic_kernel( double, double);

void bicubic_one(double *rdata, double *idata, double x, double y, double *cz)
{
	int i, j, ij;
	double wx[4], wy[4];
        double arg, w, wsum, rsum, isum;
	double a = -0.3;

	/* These weights are based on the cubic convolution kernel, see for example
	   http://undergraduate.csse.uwa.edu.au/units/CITS4241/Handouts/Lecture04.html
	   These weights include a free parameter (a).
	*/

	for (i = 0; i < 4; i++) {
		arg = fabs(x + 1 - i);
		wx[i] = cubic_kernel(arg,a);
		arg = fabs(y + 1 - i);
		wy[i] = cubic_kernel(arg,a);
	}

	rsum = isum = wsum = 0.0;
	ij = 0;
	for (j = 0; j < 4; j++) {
		for (i = 0; i < 4; i++) {
		w = wx[i] * wy[j];
		rsum += rdata[ij+i] * w;
		isum += idata[ij+i] * w;
		wsum += w;
		}
		ij += 4;
	}
	if(wsum <= 0.0) printf(" error wsum is zero \n"); 
	cz[0] = rsum/wsum;
	cz[1] = isum/wsum;
}
/************************************************************************
  kernel computes a bi-cubic spline kernel using the formula given at 
  the following web page
  http://undergraduate.csse.uwa.edu.au/units/CITS4241/Handouts/Lecture04.html
************************************************************************/
/************************************************************************
* Creator: David Sandwell       (Scripps Institution of Oceanography)   *
* Date   : 03/22/13                                                     *
************************************************************************/
#include <stdio.h>

double cubic_kernel(double arg, double a)

/* note arg must be positive and a must be between -3 and 0. */
{

double arg2, arg3, f;

	arg2=arg*arg;
	arg3=arg2*arg;
	if(arg <= 1.){
		f=(a+2)*arg3-(a+3)*arg2+1.;
	}
	else if(arg <=2.){
		f=a*arg3-5*a*arg2+8*a*arg-4*a;
	}
	else {
		f=0.;
	}
        return (f);
}

/************************************************************************
* nearest, bilinear, and bicubic interpolations                         *
************************************************************************/
/************************************************************************
* Creator: David Sandwell       (Scripps Institution of Oceanography)   *
* Date   : 03/22/13                                                     *
************************************************************************/
#include <math.h>
#include <stdio.h>
#include "gmtsar.h"

void nearest(double *ras, short *s_in, int ydims, int xdims, short *sout)
{
	int i, j, k;
        short *tmp_sin;

	/* compute the indices of the upper left corner */

	j = (int)(ras[0] + 0.5);
	i = (int)(ras[1] + 0.5);
	//k = 2*xdims*i + 2*j;
	k = 2*j;

        /* shift the pointer by i0-ns2 lines  */
        tmp_sin = s_in;
        tmp_sin++;
        tmp_sin = s_in+(size_t)(2*xdims)*(size_t)i*(tmp_sin-s_in);

	/* use the nearest point if it is within the bounds of the slave array */

	if(i < 0 || i >= ydims || j < 0 || j >= xdims) {
	  sout[0] = 0;
	  sout[1] = 0;
	}
	else {
	  //sout[0] = s_in[k];
	  //sout[1] = s_in[k+1];
	  sout[0] = tmp_sin[k];
	  sout[1] = tmp_sin[k+1];
	}
}

void bilinear(double *ras, short *s_in, int ydims, int xdims, short *sout)
{
	double dr, da, real, imag;
	int k00, k01, k10, k11;
	int i0, j0;
	int nclip;
        short *tmp_sin;

	/* compute the residual offsets */
	nclip = 0;
	j0 = (int)floor(ras[0]);
	i0 = (int)floor(ras[1]);
        dr = ras[0] - (double)j0;
        da = ras[1] - (double)i0;
	if(dr < 0. || dr > 1. || da < 0. || da > 1) fprintf(stderr," dr or da out of bounds %f %f \n",dr,da);

        /* shift the pointer by i0-ns2 lines  */
        tmp_sin = s_in;
        tmp_sin++;
        tmp_sin = s_in+(size_t)(2*xdims)*(size_t)i0*(tmp_sin-s_in);

	/* compute the indices of the 4 corners */
        /*
	k00 = 2*xdims*i0     + 2*j0;
	k01 = 2*xdims*i0     + 2*(j0+1);
	k10 = 2*xdims*(i0+1) + 2*j0;
	k11 = 2*xdims*(i0+1) + 2*(j0+1);
        */
	k00 = 2*j0;
	k01 = 2*(j0+1);
	k10 = 2*xdims + 2*j0;
	k11 = 2*xdims + 2*(j0+1);

	/* do the interpolation if all 4 corners are within the bounds of the slave array */

	if(i0 < 0 || i0 >= (ydims-1) || j0 < 0 || j0 >= (xdims-1)) {
	  sout[0] = 0;
	  sout[1] = 0;
	}
	else {
        /*real = s_in[k00] * (1.0 - da) * (1.0 - dr)
             + s_in[k10] * (da)       * (1.0 - dr)
             + s_in[k01] * (1.0 - da) * (dr) 
             + s_in[k11] * (da)       * (dr);*/
        real = tmp_sin[k00] * (1.0 - da) * (1.0 - dr)
             + tmp_sin[k10] * (da)       * (1.0 - dr)
             + tmp_sin[k01] * (1.0 - da) * (dr) 
             + tmp_sin[k11] * (da)       * (dr);

        if((int)fabs(real) > I2MAX) nclip = nclip + 1;
	sout[0] = (short)clipi2(real + 0.5);
        /*imag = s_in[k00+1] * (1.0 - da) * (1.0 - dr)
             + s_in[k10+1] * (da)       * (1.0 - dr)
             + s_in[k01+1] * (1.0 - da) * (dr) 
             + s_in[k11+1] * (da)       * (dr);*/
        imag = tmp_sin[k00+1] * (1.0 - da) * (1.0 - dr)
             + tmp_sin[k10+1] * (da)       * (1.0 - dr)
             + tmp_sin[k01+1] * (1.0 - da) * (dr) 
             + tmp_sin[k11+1] * (da)       * (dr);

        if((int)fabs(imag) > I2MAX) nclip = nclip + 1;
	sout[1] = (short)clipi2(imag + 0.5);
	}
	/*if(nclip > 0) fprintf(stderr," %d integers were clipped \n",nclip);*/
}


void bicubic(double *ras, short *s_in, int ydims, int xdims, short *sout)
{
	double dr, da;
	double rdata[16], idata[16], cz[2];
	int i, j, k, kk;
	int i0, j0;
	int nclip;
        short *tmp_sin;

	/* compute the residual offsets */
	nclip = 0;
	j0 = (int)floor(ras[0]);
	i0 = (int)floor(ras[1]);
        dr = ras[0] - (double)j0;
        da = ras[1] - (double)i0;
	if(dr < 0. || dr > 1. || da < 0. || da > 1) fprintf(stderr," dr or da out of bounds %f %f \n",dr,da);

	/* make sure all 4 corners are within the bounds of the slave array */

	if((i0-1) < 0 || (i0+2) >= ydims || (j0-1) < 0 || (j0+2) >= xdims) {
	  sout[0] = 0;
	  sout[1] = 0;
	}
	else {

	/* safe to do the interpolation */

        /* shift the pointer by i0-1 lines  */
        tmp_sin = s_in;
        tmp_sin++;
        tmp_sin = s_in+(size_t)(2*xdims)*(size_t)(i0-1)*(tmp_sin-s_in);

	for (i=0; i<4; i++){
		for (j=0; j<4; j++){
		k = i*4 +j;
		//kk = 2*xdims*(i0-1+i)  + 2*(j0-1+j);
		kk = 2*xdims*i  + 2*(j0-1+j);
		//rdata[k] = s_in[kk];
		//idata[k] = s_in[kk+1];
                rdata[k] = tmp_sin[kk];
                idata[k] = tmp_sin[kk+1];
		}
	}

	/* interpolate the real and imaginary data */

	bicubic_one(rdata, idata, dr, da, cz);

        if((int)fabs(cz[0]) > I2MAX) nclip = nclip + 1;
	sout[0] = (short)clipi2(cz[0] + 0.5);
        if((int)fabs(cz[1]) > I2MAX) nclip = nclip + 1;
	sout[1] = (short)clipi2(cz[1] + 0.5);
	}
	/*if(nclip > 0) fprintf(stderr," %d integers were clipped \n",nclip);*/
}


void bisinc (double *ras, short *s_in, int ydims, int xdims, short *sout)
{
	double dr, da, ns2=NS/2-1;
	double rdata[NS*NS], idata[NS*NS], cz[2];
	int i, j, k, kk;
	int i0, j0;
	int nclip;
        short *tmp_sin;

	/* compute the residual offsets */
	nclip = 0;
	j0 = (int)floor(ras[0]);
	i0 = (int)floor(ras[1]);
        dr = ras[0] - (double)j0;
        da = ras[1] - (double)i0;
	if(dr < 0. || dr > 1. || da < 0. || da > 1) fprintf(stderr," dr or da out of bounds %f %f \n",dr,da);

	/* make sure all 4 corners are within the bounds of the slave array */

	if((i0-ns2) < 0 || (i0+ns2+1) >= ydims || (j0-ns2) < 0 || (j0+ns2+1) >= xdims) {
	  sout[0] = 0;
	  sout[1] = 0;
	}
	else {

	/* safe to do the interpolation */

        /* shift the pointer by i0-ns2 lines  */
        tmp_sin = s_in;
        tmp_sin++;
        tmp_sin = s_in+(size_t)(2*xdims)*(size_t)(i0-ns2)*(tmp_sin-s_in);

	for (i=0; i<NS; i++){
		for (j=0; j<NS; j++){
		k = i*NS +j; 
		//kk = 2*xdims*(i0-ns2+i)  + 2*(j0-ns2+j);
		kk = 2*xdims*i  + 2*(j0+j-ns2);
		rdata[k] = tmp_sin[kk];
		idata[k] = tmp_sin[kk+1];
		}
	}

	/* interpolate the real and imaginary data */

	sinc_one(rdata, idata, dr, da, cz);

        if((int)fabs(cz[0]) > I2MAX) nclip = nclip + 1;
	sout[0] = (short)clipi2(cz[0] + 0.5);
        if((int)fabs(cz[1]) > I2MAX) nclip = nclip + 1;
	sout[1] = (short)clipi2(cz[1] + 0.5);
	}
	/*if(nclip > 0) fprintf(stderr," %d integers were clipped \n",nclip);*/
}
#include "gmtsar.h"
#include "lib_functions.h"

/*--------------------------------------------------------------*/
void print_prm_params(struct PRM p1, struct PRM p2)
{
	fprintf(stderr," SLC 1: num_rng_bins %d num_lines %d \n",p1.num_rng_bins, p1.num_lines);
	fprintf(stderr," SLC 2: num_rng_bins %d num_lines %d \n",p2.num_rng_bins, p2.num_lines);
	fprintf(stderr," lambda %f \n",p2.lambda);
	fprintf(stderr," near_range %f \n",p2.near_range);
	fprintf(stderr," rng_samp_rate %.7f \n",p2.fs);
	fprintf(stderr," sc_clock_start %f \n",p2.SC_clock_start);
	fprintf(stderr," sc_clock_stop %f \n",p2.SC_clock_stop);
	fprintf(stderr," clock_start %f \n",p2.clock_start);
	fprintf(stderr," clock_stop %f \n",p2.clock_stop);
	fprintf(stderr," prfm %f \n",p1.prf);
	fprintf(stderr," prfs %f \n",p2.prf);
	fprintf(stderr," rshift %f \n",p2.rshift+p2.sub_int_r);
	fprintf(stderr," ashift %f \n",p2.ashift+p2.sub_int_a);
}

/*--------------------------------------------------------------*/
void fix_prm_params(struct PRM *p, char *s)
{
double delr;

	delr = SOL/p->fs/2.0;

	/* these are from prm2gips */
	p->near_range = p->near_range + (p->st_rng_bin - p->chirp_ext + p->rshift-1)*delr;
	p->SC_clock_start = p->SC_clock_start + p->ashift/(p->prf*86400.0) + (p->nrows-p->num_valid_az)/(2.0*p->prf*86400);
	p->SC_clock_stop  = p->SC_clock_start + (p->num_valid_az*p->num_patches)/(p->prf*86400.0);

}
/************************************************************************
* ram2ras maps range and azimuth from a master image location into the  *
* corresponding range and azimuth location of the slave image.          *
************************************************************************/
/************************************************************************
* Creator: David Sandwell       (Scripps Institution of Oceanography)   *
* Date   : 03/22/13                                                     *
************************************************************************/
/************************************************************************
* Modification history:                                                 *
*                                                                       *
* DATE                                                                  *
************************************************************************/
#include "gmtsar.h"
#include <stdio.h>
#include <stdlib.h>

void ram2ras(struct PRM ps, double *ram, double *ras)
{
	/* this is the range coordinate */
	ras[0] = ram[0] + ((ps.rshift+ps.sub_int_r)+ram[0]*ps.stretch_r+ram[1]*ps.a_stretch_r);

	/* this is the azimuth coordinate */
	ras[1] = ram[1] + ((ps.ashift+ps.sub_int_a)+ram[0]*ps.stretch_a+ram[1]*ps.a_stretch_a);
}

/************************************************************************
  computes sinc function kernel for interpolation
************************************************************************/
/************************************************************************
* Creator: David Sandwell       (Scripps Institution of Oceanography)   *
* Date   : 03/28/13                                                     *
************************************************************************/
#include <math.h>
#define PI 3.1415926535897932

double sinc_kernel (double x)
{
double arg, f;

	arg = fabs(PI*x);
        if(arg > 0.){
		f=sin(arg)/arg;
	}
	else {
		f=1.;
	}
        return (f);
}
/************************************************************************
* sinc function interpolation                                           *
************************************************************************/
/************************************************************************
* Creator: David Sandwell       (Scripps Institution of Oceanography)   *
* Date   : 03/28/13                                                     *
************************************************************************/
#include <math.h>
#include <stdio.h>
#include "gmtsar.h"

double sinc_kernel( double );

void sinc_one(double *rdata, double *idata, double x, double y, double *cz)
{
	int i, j, ij, ns2 = NS/2-1 ;
	double wx[NS], wy[NS];
        double arg, w, wsum, rsum, isum;

	for(i = 0; i < NS ; i++){
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
		rsum += rdata[ij+i] * w;
		isum += idata[ij+i] * w;
		wsum += w;
		}
		ij += NS;
	}
	if(wsum <= 0.0) printf(" error wsum is zero \n"); 
	cz[0] = rsum/wsum;
	cz[1] = isum/wsum;
}
