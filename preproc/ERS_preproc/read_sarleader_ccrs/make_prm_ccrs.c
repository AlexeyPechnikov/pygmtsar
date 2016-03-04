/* write a PRM file */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/SARtape.h"
#include "../include/soi.h"
#define FACTOR 1000000

make_prm(sar)
struct	SAR_info	sar;
{
double	radar_wavelength,pulse_dur,rng_samp_rate;
double	I_mean,Q_mean,sc_clock_start,chirp_slope;

int     SC_identity; /* ERS-1 = 1, ERS2 = 2 */

/* translate from structure holding SAR info */
/* see SARtape.h and associated include files for more details */
/* KW - removed dpaf from sar.dpaf_dss 12/17/98 */
	sscanf((sar.dss->mission_identifier+4),"%d",&SC_identity);

	prf1 = atof(sar.dss->nominal_prf);
	radar_wavelength = atof(sar.dss->radar_wavelength);

/* convert into seconds from MHz */

	pulse_dur = (atof(sar.dss->range_pulse_length)/FACTOR);
	rng_samp_rate = FACTOR*(atof(sar.dss->sampling));

/* correct by factor of two */

	chirp_slope = 2.0*atof(sar.dss->range_pulse_phase_quad);

	I_mean = atof(sar.dss->dc_bias_i);
	Q_mean = atof(sar.dss->dc_bias_q);

/* pj: rsat leader has no values for dc_*; set to 15.5 (Sandwell) */
	
	if(I_mean == 0.)I_mean = 15.5;
	if(Q_mean == 0.)Q_mean = 15.5;

/* we don't want the first two character; (of year)) */
	sc_clock_start = atof((sar.dss->satelite_clock_time+2));
	

fprintf(stdout,"I_mean			= %lf\n",I_mean);
fprintf(stdout,"Q_mean			= %lf\n",Q_mean);
fprintf(stdout,"rng_samp_rate		= %lg\n",rng_samp_rate);
fprintf(stdout,"chirp_slope		= %lg\n",chirp_slope);
fprintf(stdout,"pulse_dur		= %lg\n",pulse_dur);
fprintf(stdout,"radar_wavelength		= %lg\n",radar_wavelength);
fprintf(stdout, "SC_identity          = %d\n",SC_identity);
}
