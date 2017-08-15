/*	$Id: offset_topo.c 79 2013-06-10 23:43:27Z pwessel $	*/
/***************************************************************************/
/* offset_topo reads  an amplitude image and a topo_ra grid as well as     */
/* an initial guess on how to shift the topo_ra to match the master        */
/* The program uses cross correlation to estimate the refined shift to     */
/* make the topo_ra match the amplitude image more exactly.                */
/* There is an option to output a new shifted topo_ra.                     */
/***************************************************************************/

/***************************************************************************
 * Creator:  David T. Sandwell and Xiaopeng Tong                           *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  7/22/08                                                       *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * DATE                                                                    *
 ***************************************************************************/

#include "gmtsar.h"

int main (int argc, char **argv)
{
	int	i,j,k,i1,j1,k1;
	int	is,js,ns;
	int	ni,nj,ntot;
	int	xshft, yshft, ib = 200;  /* ib is the width of the edge of the images not used for corr. must be > 2 */
	int	imax = 0,jmax = 0;
	double	ra,rt,avea,suma,sumt,sumc,corr,denom,maxcorr=-9999.;
	void	*API = NULL; /* GMT control structure */
	struct	GMT_GRID *A = NULL, *T = NULL, *TS = NULL;	/* Grid structure containing ->header and ->data */

	/* get the information from the command line */
	if (argc < 6){
		printf("\noffset_topo [GMTSAR] - Determine topography offset\n \n");
		printf("\nUsage: offset_topo amp_master.grd topo_ra.grd rshift ashift ns [topo_shift.grd] \n \n");
		printf("   amp_master.grd - amplitude image of master \n");
		printf("   topo_ra.grd    - topo in range/azimuth coordinates of master \n");
		printf("   rshift         - guess at integer range shift \n");
		printf("   ashift         - guess at integer azimuth shift \n");
		printf("   ns             - integer search radius \n");
		printf("   topo_shift.grd - shifted topo_ra - optional, will be shifted by rshift, ashift \n \n");
		exit(-1);
	}

	/* Begin: Initializing new GMT session */
	if ((API = GMT_Create_Session (argv[0], 0U, 0U, NULL)) == NULL) return EXIT_FAILURE;

	xshft = atoi(argv[3]);
	yshft = atoi(argv[4]);
	ns = atoi(argv[5]);

	/* Get header from amplitude and topo grids */ 
	if ((A = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, argv[1], NULL)) == NULL) return EXIT_FAILURE;
	if ((T = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, argv[2], NULL)) == NULL) return EXIT_FAILURE;

/* make sure the dimensions match */
 	if (A->header->nx != T->header->nx) { 
                fprintf (stderr, "file dimensions do not match (must have same width)\n");
                exit (EXIT_FAILURE);
        }

 	if (argc >= 7) {
		if ((TS = GMT_Create_Data (API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, A->header->wesn, A->header->inc, \
			A->header->registration, GMT_NOTSET, NULL)) == NULL) return EXIT_FAILURE;
	}

	/* Read the two grids into A->data and T->data which automatically are allocated */
	if (GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, argv[1], A) == NULL) return EXIT_FAILURE;
	if (GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, argv[2], T) == NULL) return EXIT_FAILURE;
	if (A->header->ny < T->header->ny)
		ni = A->header->ny;
	else
		ni = T->header->ny;
	fprintf(stderr," %d %d %d \n",ni,A->header->ny,T->header->ny);
        nj=T->header->nx;

/* compute average */
	ntot=0;
	suma=sumt=0.0;
	for (i=0;i<ni;i++){
		for(j=0;j<nj;j++){
			k=i*nj+j;
			ntot++;
			suma=suma+A->data[k];
			sumt=sumt+T->data[k];
		}
	}
	avea=suma/ntot;

/*  compute the normalized cross correlation function */
	for(is=-ns+yshft; is<ns+1+yshft; is++){
	for(js=-ns+xshft; js<ns+1+xshft; js++){
		ntot=0;
		sumc=suma=sumt=0.0;
		for (i=0+ib;i<ni-ib;i++){
			i1=i-is;
			for(j=0+ib;j<nj-ib;j++){
				j1=j-js;
				k=i*nj+j;
				k1=i1*nj+j1;
				if(i1 >=0 && i1 < ni && j1 >=0 && j1 < nj) {
				ntot++;
				ra=A->data[k]-avea;
                                rt=T->data[k1+1]-T->data[k1-1];
				sumc=sumc+ra*rt;
				suma=suma+ra*ra;
				sumt=sumt+rt*rt;
				}
			}
		}
		corr=0;
		denom=suma*sumt;
		if(denom > 0.) corr=sumc/sqrt(denom);
		/*printf(" rshift = %d  ashift = %d  correlation = %f\n",js,is,corr);*/
                if(corr > maxcorr) {
			maxcorr=corr;
			imax=is;
			jmax=js;
		}
	}
	}
	printf(" optimal: rshift = %d  ashift = %d  max_correlation = %f\n",jmax,imax,maxcorr);


	if (argc >= 7) {	/* write the shifted topo phase file */
		for (i=0;i<ni;i++){
			i1=i-imax; 
			for(j=0;j<nj;j++){
				j1=j-jmax;  
				k=i*nj+j;
				k1=i1*nj+j1;
				TS->data[k]=0.0f;
				if(i1 >=0 && i1 < ni && j1 >=0 && j1 < nj) TS->data[k]=T->data[k1];
			}
		}

		/*   write the shifted grd-file */
		if (GMT_Write_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, argv[6], TS)) return EXIT_FAILURE;
	}
	
	if (GMT_Destroy_Session (API)) return EXIT_FAILURE;	/* Remove the GMT machinery */

	return (EXIT_SUCCESS);
}
