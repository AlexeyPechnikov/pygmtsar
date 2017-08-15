/*	$Id: phase2topo.c 74 2013-04-21 02:57:20Z pwessel $	*/
/***************************************************************************/
/* phase2topo reads residual phase and computes residual topography.       */
/***************************************************************************/

/***************************************************************************
 * Creator:  David T. Sandwell                                             *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  10 July, 1999                                                 *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * DATE: Nov 23 2012 - modified to Use GMT i/o functions [P.Wessel]       *
 * DATE: Nov 9 2010  - modified to read and write grd files                *
 *                   - use a linear relationship between the stacked       *
 *                     residual phase and DEM corrections                  *
 *                   - this linear relationship is valid because the       *
 *                     deviations from the SRTM is small (a few meters)    *
 *                                                         by Xiaopeng     *
 ***************************************************************************/
#include "gmtsar.h"

char *USAGE = "phase2topo [GMTSAR] - Compute residual topography\n\n"
"Usage: phase2topo master.PRM topo_in.grd res_phase.grd topo_out.grd \n \n"
"    master.PRM    - master PRM files used for mapping \n"
"    topo_in.grd   - name of input topography in the radar co-ordinates of the master. \n"
"    res_phase.grd - name of input phase per unit baseline\n"
"    topo_out.grd  - name of output corrected topography in the radar co-ordinates of the master. \n\n"
"    Note the residual phase should be scaled by the perpendicular baseline (see bperp). \n";

void calc_phase2topo(int xdim, double range0, double drange, double re, double height, float *scale, float *res)
{
	int     k;
	double rho,etat;
	double sint,c,c2,ret,ret2;

        c = re + height;
	c2 = c * c;

        for (k=0; k<xdim; k++){
                ret = re+scale[k];
		ret2 = ret * ret;
                rho = range0 + k*drange;
                etat = ((rho*rho + c2 - ret2)/(2.*rho*c));
                if ( etat >=1. ) die("eta >= 0","");

                sint = sqrt(1. - etat*etat);
                scale[k] = (float)(res[k] * rho * c * sint / ret);

/*              printf("ret = %f res[k] = %f rho = %f c = %f sint = %f scale = %f \n", ret, res[k], rho, c, sint, scale[k]); */
                }
}

int main (int argc, char **argv)
{
	unsigned int	col, row;
	uint64_t node;
	int ixdec; 
	float	*scale = NULL, *res = NULL; 
	double	drange, cnst; 
	struct	PRM prm; 
	void	*API = NULL; /* GMT control structure */
	struct	GMT_GRID *P = NULL, *T = NULL;	/* Grid structure containing ->header and ->data */
	
	debug = 0;

	if (argc < 5) die("\n", USAGE);

	/* Begin: Initializing new GMT session */
	if ((API = GMT_Create_Session (argv[0], 0U, 0U, NULL)) == NULL) return EXIT_FAILURE;

	/* Get header from grd rat file */ 
	if ((T = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, argv[2], NULL)) == NULL) return EXIT_FAILURE;

	/* Get header from residual phase file */ 
	if ((P = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, argv[3], NULL)) == NULL) return EXIT_FAILURE;

	/* Make sure the topo and residual phase have the same dimensions */
 	if (!(T->header->nx == P->header->nx && T->header->ny == P->header->ny)) { 
		die ("\n", "topo and residue phase should have the same dimensions. \n"); 
        }

	get_prm (&prm, argv[1]);	/* Open and read PRM file */ 

	/* allocate the memory for work arrays */
	scale = malloc (T->header->nx * sizeof(float));
	res = malloc (T->header->nx * sizeof(float));

	ixdec = (int) prm.num_rng_bins/T->header->nx; 
	drange = ixdec*SOL/(2.0*prm.fs);
	cnst = -prm.lambda/4.0/M_PI;

	/* Read the two grids into P->data and T->data which automatically are allocated */
	if (GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, argv[2], T) == NULL) return EXIT_FAILURE;
	if (GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, argv[3], P) == NULL) return EXIT_FAILURE;

	for (row = 0; row < T->header->ny; row++) {	/* For each row */
		for (col = 0; col < T->header->nx; col++) { 	/* For each column, get node number */
			node = GMT_Get_Index (API, T->header, row, col);
			scale[col] = T->data[node];
			res[col] = P->data[node];
		}
		calc_phase2topo (T->header->nx, prm.near_range, drange, prm.RE, prm.ht, scale, res); 
		for (col = 0; col < T->header->nx; col++) { 	/* For each column, get node number */
			node = GMT_Get_Index (API, T->header, row, col);
			T->data[node] += (float)(cnst * scale[col] + T->data[node]);
		}
	}

	/* Write the output grd file */ 
	if (GMT_Write_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, argv[4], T)) return EXIT_FAILURE;

	if (GMT_Destroy_Session (API)) return EXIT_FAILURE;	/* Remove the GMT machinery */

	return (EXIT_SUCCESS);
}
