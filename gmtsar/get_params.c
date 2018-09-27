/************************************************************************
 * get_params retrieves parameters from the esar parameter file. *
 ************************************************************************/
/************************************************************************
 * Creator: Evelyn J. Price 	(Scripps Institution of Oceanography)	*
 * Date   : 11/18/96							*
 ************************************************************************/
/************************************************************************
 * Modification History: * 04/21/97 Added sc_clock_start,  sc_stop_clock and
 *sc_identity		* 06/12/00 Added icu_start * 04/26/06 Added nrows,
 *num_lines					* 06/09/07 Added
 *orbdir,rcc,raa,bpara,bperp                             * 08/01/07 Addes rshift
 *and ashift to be compatible with Mellors code   * 07/10/08 Added height0 and
 *heightf to make phase continuous           * 09/19/08 Added TEC_start and
 *TEC_end to account for ionosphere        * Date
 **
 ************************************************************************/

#include "soi.h"
#include <stdlib.h>
#include <string.h>
#define BUFSIZE 80

char *p2strcpy(char strng[]) {
	char *charptr;

	if ((charptr = (void *)malloc(strlen(strng) + 1)) != NULL)
		strcpy(charptr, strng);
	return (charptr);
}

void get_params(FILE *fh) {

	char variable[BUFSIZE], value[BUFSIZE];

	while (fscanf(fh, "%s = %s \n", variable, value) != EOF) {

		if (strcmp(variable, "input_file") == 0)
			input_file = p2strcpy(value);
		if (strcmp(variable, "led_file") == 0)
			led_file = p2strcpy(value);
		if (strcmp(variable, "out_amp_file") == 0)
			out_amp_file = p2strcpy(value);
		if (strcmp(variable, "out_data_file") == 0)
			out_data_file = p2strcpy(value);
		if (strcmp(variable, "bytes_per_line") == 0)
			bytes_per_line = atoi(value);
		if (strcmp(variable, "good_bytes_per_line") == 0)
			good_bytes = atoi(value);
		if (strcmp(variable, "first_line") == 0)
			first_line = atoi(value);
		if (strcmp(variable, "num_patches") == 0)
			num_patches = atoi(value);
		if (strcmp(variable, "first_sample") == 0)
			first_sample = atoi(value);
		if (strcmp(variable, "num_valid_az") == 0)
			num_valid_az = atoi(value);
		if (strcmp(variable, "icu_start") == 0)
			icu_start = atof(value);
		if (strcmp(variable, "SC_clock_start") == 0)
			SC_clock_start = atof(value);
		if (strcmp(variable, "SC_clock_stop") == 0)
			SC_clock_stop = atof(value);
		if (strcmp(variable, "SC_identity") == 0)
			SC_identity = atoi(value);
		if (strcmp(variable, "clock_start") == 0)
			clock_start = atof(value);
		if (strcmp(variable, "clock_stop") == 0)
			clock_stop = atof(value);
		if (strcmp(variable, "ref_identity") == 0)
			ref_identity = atoi(value);
		if (strcmp(variable, "deskew") == 0)
			deskew = p2strcpy(value);
		if (strcmp(variable, "caltone") == 0)
			caltone = atof(value);
		if (strcmp(variable, "st_rng_bin") == 0)
			st_rng_bin = atoi(value);
		if (strcmp(variable, "num_rng_bins") == 0)
			num_rng_bins = atoi(value);
		if (strcmp(variable, "earth_radius") == 0)
			RE = atof(value);
		if (strcmp(variable, "SC_vel") == 0)
			vel1 = atof(value);
		if (strcmp(variable, "SC_height") == 0)
			ht1 = atof(value);
		if (strcmp(variable, "near_range") == 0)
			near_range = atof(value);
		if (strcmp(variable, "PRF") == 0)
			prf1 = atof(value);
		if (strcmp(variable, "I_mean") == 0)
			xmi1 = atof(value);
		if (strcmp(variable, "Q_mean") == 0)
			xmq1 = atof(value);
		if (strcmp(variable, "Flip_iq") == 0)
			iqflip = p2strcpy(value);
		if (strcmp(variable, "offset_video") == 0)
			off_vid = p2strcpy(value);
		if (strcmp(variable, "az_res") == 0)
			az_res = atof(value);
		if (strcmp(variable, "nlooks") == 0)
			nlooks = atoi(value);
		if (strcmp(variable, "rng_samp_rate") == 0)
			fs = atof(value);
		if (strcmp(variable, "chirp_slope") == 0)
			slope = atof(value);
		if (strcmp(variable, "pulse_dur") == 0)
			pulsedur = atof(value);
		if (strcmp(variable, "chirp_ext") == 0)
			nextend = atoi(value);
		if (strcmp(variable, "scnd_rng_mig") == 0)
			srm = p2strcpy(value);
		if (strcmp(variable, "radar_wavelength") == 0)
			lambda = atof(value);
		if (strcmp(variable, "rng_spec_wgt") == 0)
			rhww = atof(value);
		if (strcmp(variable, "rm_rng_band") == 0)
			pctbw = atof(value);
		if (strcmp(variable, "rm_az_band") == 0)
			pctbwaz = atof(value);
		if (strcmp(variable, "fd1") == 0)
			fd1 = atof(value);
		if (strcmp(variable, "fdd1") == 0)
			fdd1 = atof(value);
		if (strcmp(variable, "fddd1") == 0)
			fddd1 = atof(value);
		if (strcmp(variable, "rshift") == 0)
			xshift = atoi(value);
		if (strcmp(variable, "ashift") == 0)
			yshift = atoi(value);
		if (strcmp(variable, "xshift") == 0)
			xshift = atoi(value);
		if (strcmp(variable, "yshift") == 0)
			yshift = atoi(value);
		if (strcmp(variable, "sub_int_r") == 0)
			sub_int_r = atof(value);
		if (strcmp(variable, "sub_int_a") == 0)
			sub_int_a = atof(value);
		if (strcmp(variable, "stretch_r") == 0)
			stretch_r = atof(value);
		if (strcmp(variable, "stretch_a") == 0)
			stretch_a = atof(value);
		if (strcmp(variable, "a_stretch_r") == 0)
			a_stretch_r = atof(value);
		if (strcmp(variable, "a_stretch_a") == 0)
			a_stretch_a = atof(value);
		if (strcmp(variable, "baseline_start") == 0)
			baseline_start = atof(value);
		if (strcmp(variable, "alpha_start") == 0)
			alpha_start = atof(value);
		if (strcmp(variable, "baseline_end") == 0)
			baseline_end = atof(value);
		if (strcmp(variable, "alpha_end") == 0)
			alpha_end = atof(value);
		if (strcmp(variable, "reference_file") == 0)
			ref_file = p2strcpy(value);
		if (strcmp(variable, "nrows") == 0)
			nrows = atoi(value);
		if (strcmp(variable, "num_lines") == 0)
			num_lines = atoi(value);
		/* new parameters added June 9, 2007 */
		if (strcmp(variable, "orbdir") == 0)
			orbdir = p2strcpy(value);
		if (strcmp(variable, "lookdir") == 0)
			lookdir = p2strcpy(value);
		if (strcmp(variable, "equatorial_radius") == 0)
			raa = atof(value);
		if (strcmp(variable, "polar_radius") == 0)
			rcc = atof(value);
		if (strcmp(variable, "B_parallel") == 0)
			bparaa = atof(value);
		if (strcmp(variable, "B_perpendicular") == 0)
			bperpp = atof(value);
		/* new parameters added July 11, 2008 */
		if (strcmp(variable, "SC_height_start") == 0)
			ht0 = atof(value);
		if (strcmp(variable, "SC_height_end") == 0)
			htf = atof(value);
		/* new parameters added September 18, 2008 */
		if (strcmp(variable, "TEC_start") == 0)
			TEC_start = atof(value);
		if (strcmp(variable, "TEC_end") == 0)
			TEC_end = atof(value);
	}
}
