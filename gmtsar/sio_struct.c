/************************************************************************
 * routines for handling sio struct in PRM files                        *
 ************************************************************************/
/************************************************************************
 * Creator: Rob J. Mellors, San Deigo State University                  *
 * Date   : December 18, 2007                                           *
 *                                                                      *
 * Modification history:                                                *
 *   some slight format changes - Dec 18, 2007 - RJM 		            *
 ************************************************************************/

/* null_sio_struct(struct PRM *p) sets all values to a NULL value       */
/* get_sio_struct(FILE *,struct PRM *p) reads values from FILE          */
/* put_sio_struct(struct PRM *p,FILE *,) writes all non-NULL values to FILE */

#include "gmtsar.h"
#include "lib_functions.h"

/*---------------------------------------------------------------*/
int get_prm(struct PRM *p, char *filename) {
	FILE *PRMfile;

	if ((PRMfile = fopen(filename, "r")) == NULL)
		die("can't open file", filename);

	/* set all prm parameters in structure to NULL values */
	null_sio_struct(p);

	/* read in prm parameters */
	get_sio_struct(PRMfile, p);

	return (EXIT_SUCCESS);
}
/*---------------------------------------------------------------*/
void null_sio_struct(struct PRM *p) {

	/* characters */
	strncpy(p->input_file, NULL_CHAR, 8);
	strncpy(p->SLC_file, NULL_CHAR, 8);
	strncpy(p->out_amp_file, NULL_CHAR, 8);
	strncpy(p->out_data_file, NULL_CHAR, 8);
	strncpy(p->deskew, NULL_CHAR, 8);
	strncpy(p->iqflip, NULL_CHAR, 8);
	strncpy(p->offset_video, NULL_CHAR, 8);
	strncpy(p->srm, NULL_CHAR, 8);
	strncpy(p->ref_file, NULL_CHAR, 8);
	strncpy(p->led_file, NULL_CHAR, 8);
	strncpy(p->orbdir, NULL_CHAR, 8);
	strncpy(p->lookdir, NULL_CHAR, 8);
	strncpy(p->dtype, NULL_CHAR, 8);

	/* ints	*/
	p->debug_flag = NULL_INT;
	p->bytes_per_line = NULL_INT;
	p->good_bytes = NULL_INT;
	p->first_line = NULL_INT;
	p->num_patches = NULL_INT;
	p->first_sample = NULL_INT;
	p->num_valid_az = NULL_INT;
	p->st_rng_bin = NULL_INT;
	p->num_rng_bins = NULL_INT;
	p->chirp_ext = NULL_INT;
	p->nlooks = NULL_INT;
	p->rshift = NULL_INT;
	p->ashift = NULL_INT;
	p->fdc_ystrt = NULL_INT;
	p->fdc_strt = NULL_INT;
	p->rec_start = NULL_INT;
	p->rec_stop = NULL_INT;
	p->SC_identity = NULL_INT;
	p->ref_identity = NULL_INT;
	p->nrows = NULL_INT;
	p->num_lines = NULL_INT;
	p->SLC_format = NULL_INT;

	/* doubles	*/
	p->SC_clock_start = NULL_DOUBLE;
	p->SC_clock_stop = NULL_DOUBLE;
	p->icu_start = NULL_DOUBLE;
	p->clock_start = NULL_DOUBLE;
	p->clock_stop = NULL_DOUBLE;
	p->caltone = NULL_DOUBLE;
	p->RE = NULL_DOUBLE;
	p->ra = NULL_DOUBLE;
	p->rc = NULL_DOUBLE;
	p->vel = NULL_DOUBLE;
	p->ht = NULL_DOUBLE;
	p->ht_start = NULL_DOUBLE;
	p->ht_end = NULL_DOUBLE;
	p->near_range = NULL_DOUBLE;
	p->far_range = NULL_DOUBLE;
	p->prf = NULL_DOUBLE;
	p->xmi = NULL_DOUBLE;
	p->xmq = NULL_DOUBLE;
	p->az_res = NULL_DOUBLE;
	p->fs = NULL_DOUBLE;
	p->chirp_slope = NULL_DOUBLE;
	p->pulsedur = NULL_DOUBLE;
	p->lambda = NULL_DOUBLE;
	p->rhww = NULL_DOUBLE;
	p->pctbw = NULL_DOUBLE;
	p->pctbwaz = NULL_DOUBLE;
	p->fd1 = NULL_DOUBLE;
	p->fdd1 = NULL_DOUBLE;
	p->fddd1 = NULL_DOUBLE;

	p->sub_int_r = NULL_DOUBLE;
	p->sub_int_a = NULL_DOUBLE;
	p->stretch_r = NULL_DOUBLE;
	p->stretch_a = NULL_DOUBLE;
	p->a_stretch_r = NULL_DOUBLE;
	p->a_stretch_a = NULL_DOUBLE;
	p->baseline_start = NULL_DOUBLE;
	p->baseline_center = NULL_DOUBLE;
	p->baseline_end = NULL_DOUBLE;
	p->alpha_start = NULL_DOUBLE;
	p->alpha_center = NULL_DOUBLE;
	p->alpha_end = NULL_DOUBLE;
	p->bpara = NULL_DOUBLE;
	p->bperp = NULL_DOUBLE;
	p->SLC_scale = NULL_DOUBLE;
	/* New parameters 9/25/18 -EXU */
	p->B_offset_start = NULL_DOUBLE;
	p->B_offset_center = NULL_DOUBLE;
	p->B_offset_end = NULL_DOUBLE;
};
/*--------------------------------------------------------------------*/
/*
        Read parameters into PRM structure from PRM file
        Based on get_params by Evelyn J. Price
        Modified by RJM
*/
/*--------------------------------------------------------------------*/

void get_sio_struct(FILE *fh, struct PRM *s) {
	char name[128], value[128], equal[128];
	char str[1024];

	if (debug) {
		fprintf(stderr, "get_sio_struct:\n");
		fprintf(stderr, "PRMname   (PRM value)     interpreted value\n");
	}

	while (fgets(str, 1024, fh) != NULL) {
		value[0] = '\0';
		sscanf(str, "%s %s %s", name, equal, value);

		/* strings */
		if (strcmp(name, "input_file") == 0)
			get_string(name, "input_file", value, s->input_file);
		if (strcmp(name, "led_file") == 0)
			get_string(name, "led_file", value, s->led_file);
		if (strcmp(name, "out_amp_file") == 0)
			get_string(name, "out_amp_file", value, s->out_amp_file);
		if (strcmp(name, "out_data_file") == 0)
			get_string(name, "out_data_file", value, s->out_data_file);
		if (strcmp(name, "scnd_rng_mig") == 0)
			get_string(name, "scnd_rng_mig", value, s->srm);
		if (strcmp(name, "deskew") == 0)
			get_string(name, "deskew", value, s->deskew);
		if (strcmp(name, "Flip_iq") == 0)
			get_string(name, "Flip_iq", value, s->iqflip);
		if (strcmp(name, "offset_video") == 0)
			get_string(name, "offset_video", value, s->offset_video);
		if (strcmp(name, "ref_file") == 0)
			get_string(name, "ref_file", value, s->ref_file);
		if (strcmp(name, "SLC_file") == 0)
			get_string(name, "SLC_file", value, s->SLC_file);
		if (strcmp(name, "orbdir") == 0)
			get_string(name, "orbdir", value, s->orbdir);
		if (strcmp(name, "lookdir") == 0)
			get_string(name, "lookdir", value, s->lookdir);
		if (strcmp(name, "dtype") == 0)
			get_string(name, "dtype", value, s->dtype);

		/* integers */
		if (strcmp(name, "nrows") == 0)
			get_int(name, "nrows", value, &s->nrows);
		if (strcmp(name, "num_lines") == 0)
			get_int(name, "num_lines", value, &s->num_lines);
		if (strcmp(name, "bytes_per_line") == 0)
			get_int(name, "bytes_per_line", value, &s->bytes_per_line);
		if (strcmp(name, "good_bytes_per_line") == 0)
			get_int(name, "good_bytes_per_line", value, &s->good_bytes);
		if (strcmp(name, "first_line") == 0)
			get_int(name, "first_line", value, &s->first_line);
		if (strcmp(name, "num_patches") == 0)
			get_int(name, "num_patches", value, &s->num_patches);
		if (strcmp(name, "first_sample") == 0)
			get_int(name, "first_sample", value, &s->first_sample);
		if (strcmp(name, "num_valid_az") == 0)
			get_int(name, "num_valid_az", value, &s->num_valid_az);
		if (strcmp(name, "SC_identity") == 0)
			get_int(name, "SC_identity", value, &s->SC_identity);
		if (strcmp(name, "chirp_ext") == 0)
			get_int(name, "chirp_ext", value, &s->chirp_ext);
		if (strcmp(name, "st_rng_bin") == 0)
			get_int(name, "st_rng_bin", value, &s->st_rng_bin);
		if (strcmp(name, "num_rng_bins") == 0)
			get_int(name, "num_rng_bins", value, &s->num_rng_bins);
		if (strcmp(name, "ref_identity") == 0)
			get_int(name, "ref_identity", value, &s->ref_identity);
		if (strcmp(name, "nlooks") == 0)
			get_int(name, "nlooks", value, &s->nlooks);
		if (strcmp(name, "rshift") == 0)
			get_int(name, "rshift", value, &s->rshift);
		if (strcmp(name, "ashift") == 0)
			get_int(name, "ashift", value, &s->ashift);
		/* backwards compatibility for xshift/rshift yshift/ashift */
		if (strcmp(name, "xshift") == 0)
			get_int(name, "rshift", value, &s->rshift);
		if (strcmp(name, "yshift") == 0)
			get_int(name, "ashift", value, &s->ashift);
		if (strcmp(name, "SLC_format") == 0)
			get_int(name, "SLC_format", value, &s->SLC_format);

		/* doubles */
		if (strcmp(name, "SC_clock_start") == 0)
			get_double(name, "SC_clock_start", value, &s->SC_clock_start);
		if (strcmp(name, "SC_clock_stop") == 0)
			get_double(name, "SC_clock_stop", value, &s->SC_clock_stop);
		if (strcmp(name, "icu_start") == 0)
			get_double(name, "icu_start", value, &s->icu_start);
		if (strcmp(name, "clock_start") == 0)
			get_double(name, "clock_start", value, &s->clock_start);
		if (strcmp(name, "clock_stop") == 0)
			get_double(name, "clock_stop", value, &s->clock_stop);
		if (strcmp(name, "caltone") == 0)
			get_double(name, "caltone", value, &s->caltone);
		if (strcmp(name, "earth_radius") == 0)
			get_double(name, "earth_radius", value, &s->RE);
		if (strcmp(name, "equatorial_radius") == 0)
			get_double(name, "equatorial_radius", value, &s->ra);
		if (strcmp(name, "polar_radius") == 0)
			get_double(name, "polar_radius", value, &s->rc);
		if (strcmp(name, "SC_vel") == 0)
			get_double(name, "SC_vel", value, &s->vel);
		if (strcmp(name, "SC_height") == 0)
			get_double(name, "SC_height", value, &s->ht);
		if (strcmp(name, "SC_height_start") == 0)
			get_double(name, "SC_height_start", value, &s->ht_start);
		if (strcmp(name, "SC_height_end") == 0)
			get_double(name, "SC_height_end", value, &s->ht_end);
		if (strcmp(name, "near_range") == 0)
			get_double(name, "near_range", value, &s->near_range);
		if (strcmp(name, "PRF") == 0)
			get_double(name, "PRF", value, &s->prf);
		if (strcmp(name, "I_mean") == 0)
			get_double(name, "I_mean", value, &s->xmi);
		if (strcmp(name, "Q_mean") == 0)
			get_double(name, "Q_mean", value, &s->xmq);
		if (strcmp(name, "az_res") == 0)
			get_double(name, "az_res", value, &s->az_res);
		if (strcmp(name, "rng_samp_rate") == 0)
			get_double(name, "rng_samp_rate", value, &s->fs);
		if (strcmp(name, "chirp_slope") == 0)
			get_double(name, "chirp_slope", value, &s->chirp_slope);
		if (strcmp(name, "pulse_dur") == 0)
			get_double(name, "pulse_dur", value, &s->pulsedur);
		if (strcmp(name, "radar_wavelength") == 0)
			get_double(name, "radar_wavelength", value, &s->lambda);
		if (strcmp(name, "rng_spec_wgt") == 0)
			get_double(name, "rng_spec_wgt", value, &s->rhww);
		if (strcmp(name, "rm_rng_band") == 0)
			get_double(name, "rm_rng_band", value, &s->pctbw);
		if (strcmp(name, "rm_az_band") == 0)
			get_double(name, "rm_az_band", value, &s->pctbwaz);
		if (strcmp(name, "fd1") == 0)
			get_double(name, "fd1", value, &s->fd1);
		if (strcmp(name, "fdd1") == 0)
			get_double(name, "fdd1", value, &s->fdd1);
		if (strcmp(name, "fddd1") == 0)
			get_double(name, "fddd1", value, &s->fddd1);
		if (strcmp(name, "sub_int_r") == 0)
			get_double(name, "sub_int_r", value, &s->sub_int_r);
		if (strcmp(name, "sub_int_a") == 0)
			get_double(name, "sub_int_a", value, &s->sub_int_a);
		if (strcmp(name, "stretch_r") == 0)
			get_double(name, "stretch_r", value, &s->stretch_r);
		if (strcmp(name, "stretch_a") == 0)
			get_double(name, "stretch_a", value, &s->stretch_a);
		if (strcmp(name, "a_stretch_r") == 0)
			get_double(name, "a_stretch_r", value, &s->a_stretch_r);
		if (strcmp(name, "a_stretch_a") == 0)
			get_double(name, "a_stretch_a", value, &s->a_stretch_a);
		if (strcmp(name, "baseline_start") == 0)
			get_double(name, "baseline_start", value, &s->baseline_start);
		if (strcmp(name, "baseline_center") == 0)
			get_double(name, "baseline_center", value, &s->baseline_center);
		if (strcmp(name, "baseline_end") == 0)
			get_double(name, "baseline_end", value, &s->baseline_end);
		if (strcmp(name, "alpha_start") == 0)
			get_double(name, "alpha_start", value, &s->alpha_start);
		if (strcmp(name, "alpha_center") == 0)
			get_double(name, "alpha_center", value, &s->alpha_center);
		if (strcmp(name, "alpha_end") == 0)
			get_double(name, "alpha_end", value, &s->alpha_end);
		if (strcmp(name, "SLC_scale") == 0)
			get_double(name, "SLC_scale", value, &s->SLC_scale);
		/* New parameters 9/25/18 -EXU */
		if (strcmp(name, "B_offset_start") == 0)
			get_double(name, "B_offset_start", value, &s->B_offset_start);
		if (strcmp(name, "B_offset_center") == 0)
			get_double(name, "B_offset_center", value, &s->B_offset_center);
		if (strcmp(name, "B_offset_end") == 0)
			get_double(name, "B_offset_end", value, &s->B_offset_end);
	}
}
/*--------------------------------------------------------------------------------*/
void get_string(char *s1, char *name, char *value, char *s2) {
	strcpy(s2, value);
	if (debug == 1)
		fprintf(stderr, " %s (%s) = %s\n", s1, name, value);
}
/*--------------------------------------------------------------------------------*/
void get_int(char *s1, char *name, char *value, int *iparam) {
	*iparam = atoi(value);
	if (debug == 1)
		fprintf(stderr, " %s (%s) = %s (%d)\n", s1, name, value, *iparam);
}
/*--------------------------------------------------------------------------------*/
void get_double(char *s1, char *name, char *value, double *param) {
	*param = atof(value);
	if (debug == 1)
		fprintf(stderr, " %s (%s) = %s (%lf)\n", s1, name, value, *param);
}
/*--------------------------------------------------------------------------------*/
/***************************************************************************/
void put_sio_struct(struct PRM prm, FILE *OUTFILE) {

	/* set by set_ALOS_defaults */
	if (prm.num_valid_az != NULL_INT)
		fprintf(OUTFILE, "num_valid_az   	= %d \n", prm.num_valid_az);
	if (prm.nrows != NULL_INT)
		fprintf(OUTFILE, "nrows   		= %d \n", prm.nrows);
	if (prm.first_line != NULL_INT)
		fprintf(OUTFILE, "first_line   		= %d \n", prm.first_line);
	if (strncmp(prm.deskew, NULL_CHAR, 8) != 0)
		fprintf(OUTFILE, "deskew   		= %s \n", prm.deskew);
	if (prm.caltone != NULL_DOUBLE)
		fprintf(OUTFILE, "caltone   		= %lf \n", prm.caltone);
	if (prm.st_rng_bin != NULL_INT)
		fprintf(OUTFILE, "st_rng_bin   		= %d \n", prm.st_rng_bin);
	if (strncmp(prm.iqflip, NULL_CHAR, 8) != 0)
		fprintf(OUTFILE, "Flip_iq   		= %s \n", prm.iqflip);
	if (strncmp(prm.offset_video, NULL_CHAR, 8) != 0)
		fprintf(OUTFILE, "offset_video   	= %s \n", prm.offset_video);
	if (prm.az_res != NULL_DOUBLE)
		fprintf(OUTFILE, "az_res   		= %lf \n", prm.az_res);
	if (prm.nlooks != NULL_INT)
		fprintf(OUTFILE, "nlooks   		= %d \n", prm.nlooks);
	if (prm.chirp_ext != NULL_INT)
		fprintf(OUTFILE, "chirp_ext   		= %d \n", prm.chirp_ext);
	if (strncmp(prm.srm, NULL_CHAR, 8) != 0)
		fprintf(OUTFILE, "scnd_rng_mig   	= %s \n", prm.srm);
	if (prm.rhww != NULL_DOUBLE)
		fprintf(OUTFILE, "rng_spec_wgt   	= %lf \n", prm.rhww);
	if (prm.pctbw != NULL_DOUBLE)
		fprintf(OUTFILE, "rm_rng_band   		= %lf \n", prm.pctbw);
	if (prm.pctbwaz != NULL_DOUBLE)
		fprintf(OUTFILE, "rm_az_band   		= %lf \n", prm.pctbwaz);
	if (prm.rshift != NULL_INT)
		fprintf(OUTFILE, "rshift  		= %d \n", prm.rshift);
	if (prm.ashift != NULL_INT)
		fprintf(OUTFILE, "ashift  	 	= %d \n", prm.ashift);
	if (prm.stretch_a != NULL_DOUBLE)
		fprintf(OUTFILE, "stretch_r   		= %g \n", prm.stretch_r);
	if (prm.stretch_a != NULL_DOUBLE)
		fprintf(OUTFILE, "stretch_a   		= %g \n", prm.stretch_a);
	if (prm.a_stretch_r != NULL_DOUBLE)
		fprintf(OUTFILE, "a_stretch_r   	= %g \n", prm.a_stretch_r);
	if (prm.a_stretch_a != NULL_DOUBLE)
		fprintf(OUTFILE, "a_stretch_a   	= %g \n", prm.a_stretch_a);
	if (prm.first_sample != NULL_INT)
		fprintf(OUTFILE, "first_sample   	= %d \n", prm.first_sample);
	if (prm.SC_identity != NULL_INT)
		fprintf(OUTFILE, "SC_identity   		= %d \n", prm.SC_identity);
	if (prm.fs != NULL_DOUBLE)
		fprintf(OUTFILE, "rng_samp_rate   	= %.6f \n", prm.fs);

	/* from read_ALOS_data */
	if (strncmp(prm.input_file, NULL_CHAR, 8) != 0)
		fprintf(OUTFILE, "input_file		= %s \n", prm.input_file);
	if (prm.num_rng_bins != NULL_INT)
		fprintf(OUTFILE, "num_rng_bins		= %d \n", prm.num_rng_bins);
	if (prm.bytes_per_line != NULL_INT)
		fprintf(OUTFILE, "bytes_per_line		= %d \n", prm.bytes_per_line);
	if (prm.good_bytes != NULL_INT)
		fprintf(OUTFILE, "good_bytes_per_line	= %d \n", prm.good_bytes);
	if (prm.prf != NULL_DOUBLE)
		fprintf(OUTFILE, "PRF			= %lf \n", prm.prf);
	if (prm.pulsedur != NULL_DOUBLE)
		fprintf(OUTFILE, "pulse_dur		= %e \n", prm.pulsedur);
	if (prm.near_range != NULL_DOUBLE)
		fprintf(OUTFILE, "near_range		= %lf \n", prm.near_range);
	if (prm.num_lines != NULL_INT)
		fprintf(OUTFILE, "num_lines		= %d \n", prm.num_lines);
	if (prm.num_patches != NULL_INT)
		fprintf(OUTFILE, "num_patches		= %d \n", prm.num_patches);
	if (prm.SC_clock_start != NULL_DOUBLE)
		fprintf(OUTFILE, "SC_clock_start		= %16.10lf \n", prm.SC_clock_start);
	if (prm.SC_clock_stop != NULL_DOUBLE)
		fprintf(OUTFILE, "SC_clock_stop		= %16.10lf \n", prm.SC_clock_stop);
	if (prm.clock_start != NULL_DOUBLE)
		fprintf(OUTFILE, "clock_start		= %16.12lf \n", prm.clock_start);
	if (prm.clock_stop != NULL_DOUBLE)
		fprintf(OUTFILE, "clock_stop			= %16.12lf \n", prm.clock_stop);
	if (strncmp(prm.led_file, NULL_CHAR, 8) != 0)
		fprintf(OUTFILE, "led_file		= %s \n", prm.led_file);

	/* from read_ALOS_ldrfile */
	if (strncmp(prm.orbdir, NULL_CHAR, 8) != 0)
		fprintf(OUTFILE, "orbdir	= %s \n", prm.orbdir);
	if (strncmp(prm.lookdir, NULL_CHAR, 8) != 0)
		fprintf(OUTFILE, "lookdir	= %s \n", prm.lookdir);
	if (prm.lambda != NULL_DOUBLE)
		fprintf(OUTFILE, "radar_wavelength	= %lg \n", prm.lambda);
	if (prm.chirp_slope != NULL_DOUBLE)
		fprintf(OUTFILE, "chirp_slope	= %lg \n", prm.chirp_slope);
	if (prm.fs != NULL_DOUBLE)
		fprintf(OUTFILE, "rng_samp_rate		= %.6f \n", prm.fs);
	if (prm.xmi != NULL_DOUBLE)
		fprintf(OUTFILE, "I_mean			= %lg \n", prm.xmi);
	if (prm.xmq != NULL_DOUBLE)
		fprintf(OUTFILE, "Q_mean			= %lg \n", prm.xmq);
	if (prm.vel != NULL_DOUBLE)
		fprintf(OUTFILE, "SC_vel			= %lf \n", prm.vel);
	if (prm.RE != NULL_DOUBLE)
		fprintf(OUTFILE, "earth_radius		= %lf \n", prm.RE);
	if (prm.ra != NULL_DOUBLE)
		fprintf(OUTFILE, "equatorial_radius	= %lf \n", prm.ra);
	if (prm.rc != NULL_DOUBLE)
		fprintf(OUTFILE, "polar_radius		= %lf \n", prm.rc);
	if (prm.ht != NULL_DOUBLE)
		fprintf(OUTFILE, "SC_height		= %lf \n", prm.ht);
	if (prm.ht_start != NULL_DOUBLE)
		fprintf(OUTFILE, "SC_height_start	= %lf \n", prm.ht_start);
	if (prm.ht_end != NULL_DOUBLE)
		fprintf(OUTFILE, "SC_height_end	= %lf \n", prm.ht_end);
	if (prm.fd1 != NULL_DOUBLE)
		fprintf(OUTFILE, "fd1			= %lf \n", prm.fd1);
	if (prm.fdd1 != NULL_DOUBLE)
		fprintf(OUTFILE, "fdd1			= %lf \n", prm.fdd1);
	if (prm.fddd1 != NULL_DOUBLE)
		fprintf(OUTFILE, "fddd1			= %lf \n", prm.fddd1);

	/* from calc_baseline */
	/*
	        if (prm.rshift != NULL_INT) fprintf(OUTFILE, "rshift = %d
	   \n",prm.rshift); if (prm.ashift != NULL_INT) fprintf(OUTFILE, "ashift =
	   %d\n",prm.ashift);
	*/
	if (prm.sub_int_r != NULL_DOUBLE)
		fprintf(OUTFILE, "sub_int_r               = %lf \n", prm.sub_int_r);
	if (prm.sub_int_a != NULL_DOUBLE)
		fprintf(OUTFILE, "sub_int_a               = %lf \n", prm.sub_int_a);
	if (prm.bpara != NULL_DOUBLE)
		fprintf(OUTFILE, "B_parallel              = %lf \n", prm.bpara);
	if (prm.bperp != NULL_DOUBLE)
		fprintf(OUTFILE, "B_perpendicular         = %lf \n", prm.bperp);
	if (prm.baseline_start != NULL_DOUBLE)
		fprintf(OUTFILE, "baseline_start          = %lf \n", prm.baseline_start);
	if (prm.baseline_center != NULL_DOUBLE)
		fprintf(OUTFILE, "baseline_center          = %lf \n", prm.baseline_center);
	if (prm.baseline_end != NULL_DOUBLE)
		fprintf(OUTFILE, "baseline_end            = %lf \n", prm.baseline_end);
	if (prm.alpha_start != NULL_DOUBLE)
		fprintf(OUTFILE, "alpha_start             = %lf \n", prm.alpha_start);
	if (prm.alpha_center != NULL_DOUBLE)
		fprintf(OUTFILE, "alpha_center             = %lf \n", prm.alpha_center);
	if (prm.alpha_end != NULL_DOUBLE)
		fprintf(OUTFILE, "alpha_end               = %lf \n", prm.alpha_end);
	if (prm.B_offset_start != NULL_DOUBLE)
		fprintf(OUTFILE, "B_offset_start          = %lf \n", prm.B_offset_start);
	if (prm.B_offset_center != NULL_DOUBLE)
		fprintf(OUTFILE, "B_offset_center         = %lf \n", prm.B_offset_center);
	if (prm.B_offset_end != NULL_DOUBLE)
		fprintf(OUTFILE, "B_offset_end            = %lf \n", prm.B_offset_end);

	/* from sarp */
	if (strncmp(prm.SLC_file, NULL_CHAR, 8) != 0)
		fprintf(OUTFILE, "SLC_file               = %s \n", prm.SLC_file);
	if (strncmp(prm.dtype, NULL_CHAR, 8) != 0)
		fprintf(OUTFILE, "dtype			= %.1s \n", prm.dtype);
	if (prm.SLC_scale != NULL_DOUBLE)
		fprintf(OUTFILE, "SLC_scale               = %lf \n", prm.SLC_scale);
}
/***************************************************************************/
