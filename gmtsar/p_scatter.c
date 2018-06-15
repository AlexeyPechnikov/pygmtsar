/*	$Id$	*/
/***************************************************************************
 * p_scatter computes the average amplitude or persistent scatter function *
 * from aligned SLCs.                                                      * 
 **************************************************************************/

/***************************************************************************
 * Creator:  David T. Sandwell                                             *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  06/06/18                                                      *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * DATE                                                                    *
 * 06/06/17     Code largely based on resamp.c                             *
 ***************************************************************************/

/* References:
Ferretti, A., Prati, C., & Rocca, F. (2001). Permanent scatterers in SAR interferometry. 
IEEE Transactions on geoscience and remote sensing, 39(1), 8-20.

Lyons, S., & Sandwell, D. (2003). Fault creep along the southern San Andreas 
from interferometric synthetic aperture radar, permanent scatterers, and stacking. 
Journal of Geophysical Research: Solid Earth, 108(B1).
*/

#include "gmtsar.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <fcntl.h>

char    *USAGE = "\nUsage: "
"p_scatter PRM_filelist fileout.grd mode \n"
"   PRM_filelist     - list of aligned SLCs \n"
"   fileout.grd      - output file either average amplitude or persistent_scatter \n"
"   mode             - (0) amplitude; (1) persistent_scatter \n \n"
"   Computes the average amplitude or persistent scattering function which is mu/(2*sig) \n"
"   The factor of 2 is used do the display_amplitude is not significantly changed. \n \n";

void read_input_file(char *, int, char **);
#define BUFSIZE 1024

int main (int argc, char **argv)
{
char    line[BUFSIZE];
char    **PRMname;
char    **SLCname;
int	i, ii, jj, kk, mm;
int	nfiles, imode = 1;
int	debug = 0 ;
int	xdimm, ydimm;		/* size of master SLC file */
short   *slc_rows = NULL;	/* pointer to a composite row of all the SLC files*/
float   *sum, *sum2;           
float   amp2, ave, sig, sig2;
FILE    *fin = NULL;
FILE    *prmfile = NULL;
FILE    **slcin = NULL;
struct PRM *r;	
double  inc[2], wesn[4];
void    *API = NULL; /* GMT control structure */
struct  GMT_GRID *scatter = NULL;        /* For the scatter grid */

	if (argc < 4) die (USAGE,"");

   /* Begin: Initializing new GMT session */
	if ((API = GMT_Create_Session (argv[0], 0U, 0U, NULL)) == NULL) return EXIT_FAILURE;

	imode = atoi(argv[3]);

   /* open the PRM_filelist count the number of lines */
        if ((fin = fopen(argv[1],"r")) == NULL) die("Can't open file",argv[1]);
        nfiles = 0;
	while(fgets(line, BUFSIZE, fin) != NULL) {
       	   nfiles++;
    	}
	fclose(fin);
        if(debug) fprintf(stderr,"number of PRM files %d \n",nfiles);

   /* allocate memory for the PRM, SLC, and the PRM structures */
        PRMname = malloc(nfiles*sizeof(char *));
        SLCname = malloc(nfiles*sizeof(char *));
        slcin   = malloc(nfiles*sizeof(FILE *));
        for (i=0; i<nfiles; i++) {
	    PRMname[i] = (char *)malloc(512*sizeof(char));
	    SLCname[i] = (char *)malloc(512*sizeof(char));
            slcin[i] = (FILE *)malloc(sizeof(FILE));
        }
        r = malloc(nfiles*sizeof(struct PRM));

   /* get the PRM filenames */
	read_input_file(argv[1], nfiles, PRMname);
        if(debug) for (i=0; i<nfiles; i++) fprintf(stderr,"%s \n",PRMname[i]);  

   /* read the first PRM file and get the dimensions of the SLC files */
	i = 0;
        if ((prmfile = fopen(PRMname[i],"r")) == NULL) die("Can't open prmfile ",PRMname[i]);
	get_sio_struct(prmfile, &r[i]);
        xdimm = r[i].num_rng_bins;
	ydimm = r[i].num_patches * r[i].num_valid_az;
	if (debug) fprintf(stderr," SLC dimensions %d %d \n",xdimm,ydimm);
	fclose(prmfile);

   /* prepare the gmt grd file */
        inc[GMT_X] = inc[GMT_Y] = 1.0;  /* Pixels */
        wesn[GMT_XLO] = 0.0;    wesn[GMT_XHI] = inc[GMT_X] * xdimm;
        wesn[GMT_YLO] = 0.0;    wesn[GMT_YHI] = inc[GMT_Y] * ydimm;
        if ((scatter = GMT_Create_Data (API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, wesn, inc, \
                GMT_GRID_PIXEL_REG, 0, NULL)) == NULL) die("could not allocate grid header","");
        if (GMT_Set_Comment (API, GMT_IS_GRID, GMT_COMMENT_IS_REMARK, "scatter", scatter)) return EXIT_FAILURE;

   /* allocate the memory for the input SLC files and output scattering (or amplitude) file */
        slc_rows = (short *)malloc(2 * xdimm * nfiles * sizeof(short));
        sum     = (float *)malloc(xdimm * sizeof(float));
        sum2    = (float *)malloc(xdimm * sizeof(float));

   /* open all the PRM files and get the names of the SLC files */
	for(kk=0; kk<nfiles; kk++) {
	   if ((prmfile = fopen(PRMname[kk],"r")) == NULL) die("Can't open prmfile ",PRMname[kk]);
           get_sio_struct(prmfile, &r[kk]);
           SLCname[kk] = r[kk].SLC_file;
           if (debug) fprintf(stderr," SLC name %s \n",SLCname[kk]);
           fclose(prmfile);
	}

   /* open all the SLC files for reading */
	for(kk=0; kk<nfiles; kk++) {
	   if ((slcin[kk] = fopen(SLCname[kk],"r")) == NULL) die("Can't open prmfile ",SLCname[kk]);
	}

   /* loopover the rows, read a row from each SLC, compute scatter, and save the output */
	for(ii=0; ii<ydimm; ii++) {

   /* read all the SLCs with the same row */
	   for(kk=0; kk<nfiles; kk++) {
	      fread(slc_rows+kk*2*xdimm, 2*sizeof(short), xdimm, slcin[kk]);
	   }
   
   /* zero all the sum and output arrays */
   /*  loop over all the columns and then the SLCs */
	   for(jj=0; jj<xdimm; jj++) {
		scatter->data[ii*xdimm+jj] = 0.;
		sum[jj] = 0.;
		sum2[jj] = 0.;
	   	for(kk=0; kk<nfiles; kk++) {
                        mm = kk*2*xdimm + jj;
           		amp2 =(float)slc_rows[mm]* (float)slc_rows[mm] + (float)slc_rows[mm+1]*(float)slc_rows[mm+1];
                        sum[jj] = sum[jj] + sqrt(amp2);
                        sum2[jj] = sum2[jj] + amp2;
		}
                ave = sum[jj]/nfiles;
                sig2 = sum2[jj]/nfiles - ave*ave;
                sig  = sqrt(sig2);

    /* output either the scatter function or the average amplitude */
                if(imode == 1) {
                   scatter->data[ii*xdimm+jj] = .1;
	           if(sig > .1) scatter->data[ii*xdimm+jj] = ave/(2.*sig2);
		}
		else {
	           scatter->data[ii*xdimm+jj]=ave;
		}
	   }
	}
	if (GMT_Write_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, argv[2], scatter)) {
                die("Failed to write output grd file","");}

   /* close all the files and remove the GMT machinery*/
	for(kk=0; kk<nfiles; kk++) {
	   fclose(slcin[kk]);
	}
	if (GMT_Destroy_Session (API)) return EXIT_FAILURE;   
	return(EXIT_SUCCESS);
}
/*---------------------------------------------------------------------------*/
void read_input_file(char *inputfilename, int nfiles, char **filename)
{
int     i;
FILE    *inputfile;
        if ((inputfile=fopen(inputfilename,"r")) == NULL) die("Can't open ",inputfilename);
        for (i=0; i<nfiles; i++) fscanf(inputfile," %s ",filename[i]);
        fclose(inputfile);
}
