/*	$Id: file_stuff.c 109 2015-01-19 23:01:24Z sandwell $	*/
#include "gmtsar.h"

/*--------------------------------------------------------------*/
void print_prm_params(struct PRM p1, struct PRM p2)
{
	fprintf(stderr," SLC 1: num_rng_bins %d num_lines %d \n",p1.num_rng_bins, p1.num_lines);
	fprintf(stderr," SLC 2: num_rng_bins %d num_lines %d \n",p2.num_rng_bins, p2.num_lines);
	fprintf(stderr," lambda %f \n",p2.lambda);
	fprintf(stderr," baseline_start %f \n",p2.baseline_start);
	fprintf(stderr," baseline_end %f \n",p2.baseline_end);
	fprintf(stderr," alpha_start %f \n",p2.alpha_start);
	fprintf(stderr," alpha_end %f \n",p2.alpha_end);
	fprintf(stderr," near_range %f \n",p2.near_range);
	fprintf(stderr," rng_samp_rate %f \n",p2.fs);
	fprintf(stderr," SC_clock_start %f \n",p2.SC_clock_start);
	fprintf(stderr," SC_clock_stop %f \n",p2.SC_clock_stop);
	fprintf(stderr," clock_start %f \n",p2.clock_start);
	fprintf(stderr," clock_stop %f \n",p2.clock_stop);
	fprintf(stderr," prf %f \n",p2.prf);
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
	p->clock_start = p->clock_start + p->ashift/(p->prf*86400.0) + (p->nrows-p->num_valid_az)/(2.0*p->prf*86400);
	p->clock_stop  = p->clock_start + (p->num_valid_az*p->num_patches)/(p->prf*86400.0);

}

