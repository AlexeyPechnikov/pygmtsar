/*	$Id: bperp.c 79 2013-06-10 23:43:27Z pwessel $	*/
/***************************************************************************
 * bperp reads the PRM file of a repeat pass file and extracts the         *
 * baseline information.  Then it computes a complete array of             *
 * perpendicular baseline that depends on both range and azimuth.          *
 * This baseline array is used for taking sums and differences of          *
 * interferograms                                                          *
 ***************************************************************************/

/***************************************************************************
 * Creator:  David T. Sandwell                                             *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  7/7/95                                                        *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *     Modified from ihBperp to read and write grd files.                  *
 *                              Xiaopeng Nov 24 2010                       *
 * DATE                                                                    *
 ***************************************************************************/

#include "gmtsar.h"

char *USAGE = "bperp [GMT5SAR] - Perpendicular baseline processing\n\n"
"   Usage: bperp rep.PRM phase.grd bperp.grd \n\n"
"   rep.PRM    -  input repeat PRM file that provides the interferometric baseline \n"   
"   phase.grd  -  input phase grd file that provides the dimensions of the grd file \n"   
"   bperp.grd  -  output grd file of the perpendicular baseline\n\n"; 

void print_prm_params (struct PRM p1, struct PRM p2)
{
	fprintf(stderr," SLC 1: num_rng_bins %d num_lines %d \n",p1.num_rng_bins, p1.num_lines);
	fprintf(stderr," SLC 2: num_rng_bins %d num_lines %d \n",p2.num_rng_bins, p2.num_lines);
	fprintf(stderr," lambda %f \n",p2.lambda);
	fprintf(stderr," baseline_start %f \n",p2.baseline_start);
	fprintf(stderr," B_end %f \n",p2.baseline_end);
	fprintf(stderr," alpha_start %f \n",p2.alpha_start);
	fprintf(stderr," alpha_end %f \n",p2.alpha_end);
	fprintf(stderr," near_range %f \n",p2.near_range);
	fprintf(stderr," rng_samp_rate %f \n",p2.fs);
	fprintf(stderr," sc_clock_start %f \n",p2.SC_clock_start);
	fprintf(stderr," sc_clock_stop %f \n",p2.SC_clock_stop);
	fprintf(stderr," clock_start %f \n",p2.clock_start);
	fprintf(stderr," clock_stop %f \n",p2.clock_stop);
	fprintf(stderr," prf %f \n",p2.prf);
}

void fix_prm_params (struct PRM *p, char *s)
{
	double delr;

	delr = SOL/p->fs/2.0;

	/* these are from prm2gips */
	p->near_range = p->near_range + (p->st_rng_bin - p->chirp_ext + p->rshift-1)*delr;
	p->SC_clock_start = p->SC_clock_start + p->ashift/(p->prf*86400.0) + (p->nrows-p->num_valid_az)/(2.0*p->prf*86400);
	p->SC_clock_stop  = p->SC_clock_start + (p->num_valid_az*p->num_patches)/(p->prf*86400.0);

}

void calc_theta (int xdim, double range0, double drange, double re, double height, double alpha, float *theta)
{
        int k;
        double rho,eta;
        double c,c2,re2;

        c = re + height;
        c2 = c*c;
        re2 = re*re;

        for(k=0;k<xdim;k++){
                rho = range0 + k*drange;
                eta = ((rho*rho + c2 - re2)/(2.*rho*c));
		if ( eta >=1. ) die("eta >= 0","");
                theta[k]=(float)acos(eta);
        }
}
	
int main (int argc, char **argv)
{
	struct	PRM prm;
	float	*theta = NULL;
	double	drange, dBh, dBv, Bh0, Bhf, Bv0, Bvf, B, alpha, Bh, Bv, dht, height;
	unsigned int	row, col;
	uint64_t node;  
	void	*API = NULL; /* GMT API control structure */
	struct	GMT_GRID *P = NULL, *BP = NULL;	/* Grid structure containing ->header and ->data */

	if (argc < 4) die("\n", USAGE); 

	/* Begin: Initializing new GMT5 session */
	if ((API = GMT_Create_Session (argv[0], 0U, 0U, NULL)) == NULL) return EXIT_FAILURE;

	/* open PRM file */ 
	get_prm (&prm, argv[1]);

	/* near_range, SC_clock_start, and SC_clock_stop need to be changed */
        fix_prm_params (&prm, argv[1]);

	/* Get header from residual phase file */ 
	if ((P = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, argv[2], NULL)) == NULL) return EXIT_FAILURE;

	/* allocate memory */
	if ((BP = GMT_Create_Data (API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, P->header->wesn, P->header->inc, \
		P->header->registration, GMT_NOTSET, NULL)) == NULL) return EXIT_FAILURE;
	
	theta = (float *) malloc (P->header->nx * sizeof(float));

	/* calculate drange */

        drange = fabs (P->header->inc[GMT_X])*SOL/(2.0*prm.fs);

	/* compute starting and ending baselines in rectangular co-ordinates */

        dBh = dBv = 0.;
        Bh0 = prm.baseline_start*cos(prm.alpha_start*PI/180.);
        Bhf = prm.baseline_end*cos(prm.alpha_end*PI/180.);
        Bv0 = prm.baseline_start*sin(prm.alpha_start*PI/180.);
        Bvf = prm.baseline_end*sin(prm.alpha_end*PI/180.);
        dBh = (Bhf - Bh0)/P->header->ny;  
        dBv = (Bvf - Bv0)/P->header->ny;    
        dBh = (Bhf - Bh0)/P->header->ny;  
        dBv = (Bvf - Bv0)/P->header->ny;    

	/* calculate height increment if available */
        dht = 0.0;        

        if (prm.ht_start > 0.0 && prm.ht_end > 0.0) {
		dht = (prm.ht_end - prm.ht_start)/P->header->ny;
	}

	for (row = 0; row < BP->header->ny; row++) {	/* For each row */

		/* change the baseline and alpha along the satellite track */

        	Bh = Bh0 + dBh*row;
        	Bv = Bv0 + dBv*row;
        	B  = sqrt(Bh*Bh + Bv*Bv);
        	alpha = atan2(Bv,Bh);
		height = prm.ht_start + dht*row;

        	calc_theta (P->header->nx,prm.near_range,drange,prm.RE,height,alpha,theta);

		/* compute the perpendicular baseline */
		for (col = 0; col < BP->header->nx; col++){ 	/* For each column, get node number */
			node = GMT_Get_Index (API, BP->header, row, col);
                	BP->data[node] = (float)(B*cos((double)theta[col]-alpha));
		}
	}
	
	/* Write the output grd file */ 
	if (GMT_Set_Comment (API, GMT_IS_GRID, GMT_COMMENT_IS_TITLE, "perpendicular baseline", BP)) return EXIT_FAILURE;
	
	if (GMT_Write_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, argv[3], BP)) return EXIT_FAILURE;

	if (GMT_Destroy_Session (API)) return EXIT_FAILURE;	/* Remove the GMT machinery */

	return (EXIT_SUCCESS);
} 
