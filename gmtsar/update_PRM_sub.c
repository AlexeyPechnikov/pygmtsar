#include "gmtsar.h"
#include "orbit.h"
#include "update_PRM.h"
#include <assert.h>
#include <ctype.h>
#include <getopt.h>
#include <limits.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct list {
	char *string;
	char *value;
	struct list *next;
};

typedef struct list LIST;
typedef enum type { CHAR, INT, DOUBLE } type;

typedef struct PRMPRIMITIVE {
	char *name;
	char *alias;
	type etype;
	int offset;
	int size;
} prmdata;

/* Lookup table for types of PRM values */
/* Make sure the number in the define below matches the number of variables that
 * get initialized in init_lookup_table() */

#define LOOKUPTABLE_SIZE 86
static struct PRMPRIMITIVE lookuptable[LOOKUPTABLE_SIZE];

char *get_PRM_sub(char *szFile, char *szValueName);
int findlongest(char *szFile, int func);

int m_strcmp(char *szStringA, char *szStringB);

/* Hack to allow use of strcmp on NULL strings without segfault */
int m_strcmp(char *szStringA, char *szStringB) {

	if (szStringA != NULL && szStringB != NULL)
		return strcmp(szStringA, szStringB);
	else
		return -9999;
}

char *get_PRM_string(char *prmfile, char *valuename) { return get_PRM(prmfile, valuename); }

int get_PRM_int(char *prmfile, char *valuename) {

	int retval = 0;
	char *value;

	value = get_PRM(prmfile, valuename);
	retval = atoi(value);
	free(value);
	return retval;
}

double get_PRM_double(char *prmfile, char *valuename) {
	double retval = 0.0;
	char *value;

	value = get_PRM(prmfile, valuename);
	retval = atof(value);
	free(value);
	return retval;
}

int update_PRM_sub(char *prmfile, char *entry, char *value) {

	FILE *prm_file;
	struct PRM prm;

	if ((prm_file = fopen(prmfile, "r")) == NULL) {
		fprintf(stderr, "\n error -  could not open file %s \n\n", prmfile);
		return 1;
	}
	else {
		init_lookup_table();
		null_sio_struct(&prm);
		get_sio_struct(prm_file, &prm);
	}

	setvalue(&prm, entry, value);

	if ((prm_file = fopen(prmfile, "w+")) == NULL) {
		return 2;
	}
	else {
		put_sio_struct(prm, prm_file);
		fclose(prm_file);
	}

	return 0;
}

void init_lookup_table() {

	/* All char variables */
	lookuptable[0].name = "input_file";
	lookuptable[0].alias = "input_file";
	lookuptable[0].etype = CHAR;
	lookuptable[0].offset = offsetof(struct PRM, input_file);
	lookuptable[0].size = 256;

	lookuptable[1].name = "SLC_file";
	lookuptable[1].alias = "SLC_file";
	lookuptable[1].etype = CHAR;
	lookuptable[1].offset = offsetof(struct PRM, SLC_file);
	lookuptable[1].size = 256;

	lookuptable[2].name = "out_amp_file";
	lookuptable[2].alias = "out_amp_file";
	lookuptable[2].etype = CHAR;
	lookuptable[2].offset = offsetof(struct PRM, out_amp_file);
	lookuptable[2].size = 256;

	lookuptable[3].name = "out_data_file";
	lookuptable[3].alias = "out_data_file";
	lookuptable[3].etype = CHAR;
	lookuptable[3].offset = offsetof(struct PRM, out_data_file);
	lookuptable[3].size = 256;

	lookuptable[4].name = "deskew";
	lookuptable[4].alias = "deskew";
	lookuptable[4].etype = CHAR;
	lookuptable[4].offset = offsetof(struct PRM, deskew);
	lookuptable[4].size = 8;

	lookuptable[5].name = "iqflip";
	lookuptable[5].alias = "Flip_iq";
	lookuptable[5].etype = CHAR;
	lookuptable[5].offset = offsetof(struct PRM, iqflip);
	lookuptable[5].size = 8;

	lookuptable[6].name = "offset_video";
	lookuptable[6].alias = "offset_video";
	lookuptable[6].etype = CHAR;
	lookuptable[6].offset = offsetof(struct PRM, offset_video);
	lookuptable[6].size = 8;

	lookuptable[7].name = "srm";
	lookuptable[7].alias = "scnd_rng_mig";
	lookuptable[7].etype = CHAR;
	lookuptable[7].offset = offsetof(struct PRM, srm);
	lookuptable[7].size = 8;

	lookuptable[8].name = "ref_file";
	lookuptable[8].alias = "ref_file";
	lookuptable[8].etype = CHAR;
	lookuptable[8].offset = offsetof(struct PRM, ref_file);
	lookuptable[8].size = 128;

	lookuptable[9].name = "led_file";
	lookuptable[9].alias = "led_file";
	lookuptable[9].etype = CHAR;
	lookuptable[9].offset = offsetof(struct PRM, led_file);
	lookuptable[9].size = 128;

	lookuptable[10].name = "orbdir";
	lookuptable[10].alias = "orbdir";
	lookuptable[10].etype = CHAR;
	lookuptable[10].offset = offsetof(struct PRM, orbdir);
	lookuptable[10].size = 8;

	lookuptable[11].name = "lookdir";
	lookuptable[11].alias = "lookdir";
	lookuptable[11].etype = CHAR;
	lookuptable[11].offset = offsetof(struct PRM, lookdir);
	lookuptable[11].size = 8;

	lookuptable[12].name = "dtype";
	lookuptable[12].alias = "dtype";
	lookuptable[12].etype = CHAR;
	lookuptable[12].offset = offsetof(struct PRM, dtype);
	lookuptable[12].size = 8;

	lookuptable[13].name = "date";
	lookuptable[13].alias = "date";
	lookuptable[13].etype = CHAR;
	lookuptable[13].offset = offsetof(struct PRM, date);
	lookuptable[13].size = 16;

	/* All int variables */

	for (int i = 14; i <= 35; i++)
		lookuptable[i].size = sizeof(int);

	lookuptable[14].name = "debug_flag";
	lookuptable[14].alias = "debug_flag";
	lookuptable[14].etype = INT;
	lookuptable[14].offset = offsetof(struct PRM, debug_flag);

	lookuptable[15].name = "bytes_per_line";
	lookuptable[15].alias = "bytes_per_line";
	lookuptable[15].etype = INT;
	lookuptable[15].offset = offsetof(struct PRM, bytes_per_line);

	lookuptable[16].name = "good_bytes";
	lookuptable[16].alias = "good_bytes_per_line";
	lookuptable[16].etype = INT;
	lookuptable[16].offset = offsetof(struct PRM, good_bytes);

	lookuptable[17].name = "first_line";
	lookuptable[17].alias = "first_line";
	lookuptable[17].etype = INT;
	lookuptable[17].offset = offsetof(struct PRM, first_line);

	lookuptable[18].name = "num_patches";
	lookuptable[18].alias = "num_patches";
	lookuptable[18].etype = INT;
	lookuptable[18].offset = offsetof(struct PRM, num_patches);

	lookuptable[19].name = "first_sample";
	lookuptable[19].alias = "first_sample";
	lookuptable[19].etype = INT;
	lookuptable[19].offset = offsetof(struct PRM, first_sample);

	lookuptable[20].name = "num_valid_az";
	lookuptable[20].alias = "num_valid_az";
	lookuptable[20].etype = INT;
	lookuptable[20].offset = offsetof(struct PRM, num_valid_az);

	lookuptable[21].name = "st_rng_bin";
	lookuptable[21].alias = "st_rng_bin";
	lookuptable[21].etype = INT;
	lookuptable[21].offset = offsetof(struct PRM, st_rng_bin);

	lookuptable[22].name = "num_rng_bins";
	lookuptable[22].alias = "num_rng_bins";
	lookuptable[22].etype = INT;
	lookuptable[22].offset = offsetof(struct PRM, num_rng_bins);

	lookuptable[23].name = "chirp_ext";
	lookuptable[23].alias = "chirp_ext";
	lookuptable[23].etype = INT;
	lookuptable[23].offset = offsetof(struct PRM, chirp_ext);

	lookuptable[24].name = "nlooks";
	lookuptable[24].alias = "nlooks";
	lookuptable[24].etype = INT;
	lookuptable[24].offset = offsetof(struct PRM, nlooks);

	lookuptable[25].name = "rshift";
	lookuptable[25].alias = "rshift";
	lookuptable[25].etype = INT;
	lookuptable[25].offset = offsetof(struct PRM, rshift);

	lookuptable[26].name = "ashift";
	lookuptable[26].alias = "ashift";
	lookuptable[26].etype = INT;
	lookuptable[26].offset = offsetof(struct PRM, ashift);

	lookuptable[27].name = "fdc_ystrt";
	lookuptable[27].alias = "fdc_ystrt";
	lookuptable[27].etype = INT;
	lookuptable[27].offset = offsetof(struct PRM, fdc_ystrt);

	lookuptable[28].name = "fdc_strt";
	lookuptable[28].alias = "fdc_strt";
	lookuptable[28].etype = INT;
	lookuptable[28].offset = offsetof(struct PRM, fdc_strt);

	lookuptable[29].name = "rec_start";
	lookuptable[29].alias = "rec_start";
	lookuptable[29].etype = INT;
	lookuptable[29].offset = offsetof(struct PRM, rec_start);

	lookuptable[30].name = "rec_stop";
	lookuptable[30].alias = "rec_stop";
	lookuptable[30].etype = INT;
	lookuptable[30].offset = offsetof(struct PRM, rec_stop);

	lookuptable[31].name = "SC_identity";
	lookuptable[31].alias = "SC_identity";
	lookuptable[31].etype = INT;
	lookuptable[31].offset = offsetof(struct PRM, SC_identity);

	lookuptable[32].name = "ref_identity";
	lookuptable[32].alias = "ref_identity";
	lookuptable[32].etype = INT;
	lookuptable[32].offset = offsetof(struct PRM, ref_identity);

	lookuptable[33].name = "nrows";
	lookuptable[33].alias = "nrows";
	lookuptable[33].etype = INT;
	lookuptable[33].offset = offsetof(struct PRM, nrows);

	lookuptable[34].name = "num_lines";
	lookuptable[34].alias = "num_lines";
	lookuptable[34].etype = INT;
	lookuptable[34].offset = offsetof(struct PRM, num_lines);

	lookuptable[35].name = "SLC_format";
	lookuptable[35].alias = "SLC_format";
	lookuptable[35].etype = INT;
	lookuptable[35].offset = offsetof(struct PRM, SLC_format);

	/* All double variables */
	for (int i = 36; i <= 82; i++)
		lookuptable[i].size = sizeof(double);

	lookuptable[36].name = "SC_clock_start";
	lookuptable[36].alias = "SC_clock_start";
	lookuptable[36].etype = DOUBLE;
	lookuptable[36].offset = offsetof(struct PRM, SC_clock_start);

	lookuptable[37].name = "SC_clock_stop";
	lookuptable[37].alias = "SC_clock_stop";
	lookuptable[37].etype = DOUBLE;
	lookuptable[37].offset = offsetof(struct PRM, SC_clock_stop);

	lookuptable[38].name = "icu_start";
	lookuptable[38].alias = "icu_start";
	lookuptable[38].etype = DOUBLE;
	lookuptable[38].offset = offsetof(struct PRM, icu_start);

	lookuptable[39].name = "clock_start";
	lookuptable[39].alias = "clock_start";
	lookuptable[39].etype = DOUBLE;
	lookuptable[39].offset = offsetof(struct PRM, clock_start);

	lookuptable[40].name = "clock_stop";
	lookuptable[40].alias = "clock_stop";
	lookuptable[40].etype = DOUBLE;
	lookuptable[40].offset = offsetof(struct PRM, clock_stop);

	lookuptable[41].name = "caltone";
	lookuptable[41].alias = "caltone";
	lookuptable[41].etype = DOUBLE;
	lookuptable[41].offset = offsetof(struct PRM, caltone);

	lookuptable[42].name = "RE";            /*local earth eadius */
	lookuptable[42].alias = "earth_radius"; /*local earth eadius */
	lookuptable[42].etype = DOUBLE;
	lookuptable[42].offset = offsetof(struct PRM, RE);

	lookuptable[43].name = "rc";            /* polar radius */
	lookuptable[43].alias = "polar_radius"; /* polar radius */
	lookuptable[43].etype = DOUBLE;
	lookuptable[43].offset = offsetof(struct PRM, rc);

	lookuptable[44].name = "ra";                 /* equatorial radius */
	lookuptable[44].alias = "equatorial_radius"; /* equatorial radius */
	lookuptable[44].etype = DOUBLE;
	lookuptable[44].offset = offsetof(struct PRM, ra);

	lookuptable[45].name = "vel";     /* Equivalent SC velocity */
	lookuptable[45].alias = "SC_vel"; /* Equivalent SC velocity */
	lookuptable[45].etype = DOUBLE;
	lookuptable[45].offset = offsetof(struct PRM, vel);

	lookuptable[46].name = "ht";         /* (SC_radius - RE) center */
	lookuptable[46].alias = "SC_height"; /* (SC_radius - RE) center */
	lookuptable[46].etype = DOUBLE;
	lookuptable[46].offset = offsetof(struct PRM, ht);

	lookuptable[47].name = "ht_start";         /* (SC_radius - RE) start */
	lookuptable[47].alias = "SC_height_start"; /* (SC_radius - RE) start */
	lookuptable[47].etype = DOUBLE;
	lookuptable[47].offset = offsetof(struct PRM, ht_start);

	lookuptable[48].name = "ht_end";         /* (SC_radius - RE) end */
	lookuptable[48].alias = "SC_height_end"; /* (SC_radius - RE) end */
	lookuptable[48].etype = DOUBLE;
	lookuptable[48].offset = offsetof(struct PRM, ht_end);

	lookuptable[49].name = "near_range";
	lookuptable[49].alias = "near_range";
	lookuptable[49].etype = DOUBLE;
	lookuptable[49].offset = offsetof(struct PRM, near_range);

	lookuptable[50].name = "far_range";
	lookuptable[50].alias = "far_range";
	lookuptable[50].etype = DOUBLE;
	lookuptable[50].offset = offsetof(struct PRM, far_range);

	lookuptable[51].name = "prf";
	lookuptable[51].alias = "PRF";
	lookuptable[51].etype = DOUBLE;
	lookuptable[51].offset = offsetof(struct PRM, prf);

	lookuptable[52].name = "xmi";
	lookuptable[52].alias = "I_mean";
	lookuptable[52].etype = DOUBLE;
	lookuptable[52].offset = offsetof(struct PRM, xmi);

	lookuptable[53].name = "xmq";
	lookuptable[53].alias = "Q_mean";
	lookuptable[53].etype = DOUBLE;
	lookuptable[53].offset = offsetof(struct PRM, xmq);

	lookuptable[54].name = "az_res";
	lookuptable[54].alias = "az_res";
	lookuptable[54].etype = DOUBLE;
	lookuptable[54].offset = offsetof(struct PRM, az_res);

	lookuptable[55].name = "fs";
	lookuptable[55].alias = "rng_samp_rate";
	lookuptable[55].etype = DOUBLE;
	lookuptable[55].offset = offsetof(struct PRM, fs);

	lookuptable[56].name = "chirp_slope";
	lookuptable[56].alias = "chirp_slope";
	lookuptable[56].etype = DOUBLE;
	lookuptable[56].offset = offsetof(struct PRM, chirp_slope);

	lookuptable[57].name = "pulsedur";
	lookuptable[57].alias = "pulsedur";
	lookuptable[57].etype = DOUBLE;
	lookuptable[57].offset = offsetof(struct PRM, pulsedur);

	lookuptable[58].name = "lambda";
	lookuptable[58].alias = "radar_wavelength";
	lookuptable[58].etype = DOUBLE;
	lookuptable[58].offset = offsetof(struct PRM, lambda);

	lookuptable[59].name = "rhww";
	lookuptable[59].alias = "rng_spec_wgt";
	lookuptable[59].etype = DOUBLE;
	lookuptable[59].offset = offsetof(struct PRM, rhww);

	lookuptable[60].name = "pctbw";
	lookuptable[60].alias = "rm_rng_band";
	lookuptable[60].etype = DOUBLE;
	lookuptable[60].offset = offsetof(struct PRM, pctbw);

	lookuptable[61].name = "pctbwaz";
	lookuptable[61].alias = "rm_az_band";
	lookuptable[61].etype = DOUBLE;
	lookuptable[61].offset = offsetof(struct PRM, pctbwaz);

	lookuptable[62].name = "fd1";
	lookuptable[62].alias = "fd1";
	lookuptable[62].etype = DOUBLE;
	lookuptable[62].offset = offsetof(struct PRM, fd1);

	lookuptable[63].name = "fdd1";
	lookuptable[63].alias = "fdd1";
	lookuptable[63].etype = DOUBLE;
	lookuptable[63].offset = offsetof(struct PRM, fdd1);

	lookuptable[64].name = "fddd1";
	lookuptable[64].alias = "fddd1";
	lookuptable[64].etype = DOUBLE;
	lookuptable[64].offset = offsetof(struct PRM, fddd1);

	lookuptable[65].name = "delr";
	lookuptable[65].alias = "delr";
	lookuptable[65].etype = DOUBLE;
	lookuptable[65].offset = offsetof(struct PRM, delr);

	lookuptable[66].name = "yaw";
	lookuptable[66].alias = "yaw";
	lookuptable[66].etype = DOUBLE;
	lookuptable[66].offset = offsetof(struct PRM, yaw);

	lookuptable[67].name = "SLC_scale";
	lookuptable[67].alias = "SLC_scale";
	lookuptable[67].etype = DOUBLE;
	lookuptable[67].offset = offsetof(struct PRM, SLC_scale);

	lookuptable[68].name = "sub_int_r";
	lookuptable[68].alias = "sub_int_r";
	lookuptable[68].etype = DOUBLE;
	lookuptable[68].offset = offsetof(struct PRM, sub_int_r);

	lookuptable[69].name = "sub_int_a";
	lookuptable[69].alias = "sub_int_a";
	lookuptable[69].etype = DOUBLE;
	lookuptable[69].offset = offsetof(struct PRM, sub_int_a);

	lookuptable[70].name = "sub_double";
	lookuptable[70].alias = "sub_double";
	lookuptable[70].etype = DOUBLE;
	lookuptable[70].offset = offsetof(struct PRM, sub_double);

	lookuptable[71].name = "stretch_r";
	lookuptable[71].alias = "stretch_r";
	lookuptable[71].etype = DOUBLE;
	lookuptable[71].offset = offsetof(struct PRM, stretch_r);

	lookuptable[72].name = "stretch_a";
	lookuptable[72].alias = "stretch_a";
	lookuptable[72].etype = DOUBLE;
	lookuptable[72].offset = offsetof(struct PRM, stretch_a);

	lookuptable[73].name = "a_stretch_r";
	lookuptable[73].alias = "a_stretch_r";
	lookuptable[73].etype = DOUBLE;
	lookuptable[73].offset = offsetof(struct PRM, a_stretch_r);

	lookuptable[74].name = "a_stretch_a";
	lookuptable[74].alias = "a_stretch_a";
	lookuptable[74].etype = DOUBLE;
	lookuptable[74].offset = offsetof(struct PRM, a_stretch_a);

	lookuptable[75].name = "baseline_start";
	lookuptable[75].alias = "baseline_start";
	lookuptable[75].etype = DOUBLE;
	lookuptable[75].offset = offsetof(struct PRM, baseline_start);

	lookuptable[76].name = "baseline_center";
	lookuptable[76].alias = "baseline_center";
	lookuptable[76].etype = DOUBLE;
	lookuptable[76].offset = offsetof(struct PRM, baseline_center);

	lookuptable[77].name = "baseline_end";
	lookuptable[77].alias = "baseline_end";
	lookuptable[77].etype = DOUBLE;
	lookuptable[77].offset = offsetof(struct PRM, baseline_end);

	lookuptable[78].name = "alpha_start";
	lookuptable[78].alias = "alpha_start";
	lookuptable[78].etype = DOUBLE;
	lookuptable[78].offset = offsetof(struct PRM, alpha_start);

	lookuptable[79].name = "alpha_center";
	lookuptable[79].alias = "alpha_center";
	lookuptable[79].etype = DOUBLE;
	lookuptable[79].offset = offsetof(struct PRM, alpha_center);

	lookuptable[80].name = "alpha_end";
	lookuptable[80].alias = "alpha_end";
	lookuptable[80].etype = DOUBLE;
	lookuptable[80].offset = offsetof(struct PRM, alpha_end);

	lookuptable[81].name = "bpara";
	lookuptable[81].alias = "B_parallel";
	lookuptable[81].etype = DOUBLE;
	lookuptable[81].offset = offsetof(struct PRM, bpara);

	lookuptable[82].name = "bperp";
	lookuptable[82].alias = "B_perpendicular";
	lookuptable[82].etype = DOUBLE;
	lookuptable[82].offset = offsetof(struct PRM, bperp);

	lookuptable[83].name = "B_offset_start";
	lookuptable[83].alias = "B_offset_start";
	lookuptable[83].etype = DOUBLE;
	lookuptable[83].offset = offsetof(struct PRM, B_offset_start);

	lookuptable[84].name = "B_offset_center";
	lookuptable[84].alias = "B_offset_center";
	lookuptable[84].etype = DOUBLE;
	lookuptable[84].offset = offsetof(struct PRM, B_offset_center);

	lookuptable[85].name = "B_offset_end";
	lookuptable[85].alias = "B_offset_end";
	lookuptable[85].etype = DOUBLE;
	lookuptable[85].offset = offsetof(struct PRM, B_offset_end);
}

int setvalue(struct PRM *prm, char *name, char *value) {
	struct PRMPRIMITIVE valuetype;
	int n = 0;
	int intvalue = 0;
	double doublevalue = 0.0;
	char *base;
	char *newbase;

	for (n = 0; n <= LOOKUPTABLE_SIZE; n++) {
		if (m_strcmp(lookuptable[n].name, name) == 0 || m_strcmp(lookuptable[n].alias, name) == 0) {
			valuetype.name = malloc(sizeof(char) * strlen(lookuptable[n].name));
			strcpy(valuetype.name, lookuptable[n].name);
			valuetype.etype = lookuptable[n].etype;
			valuetype.offset = lookuptable[n].offset;
			valuetype.size = lookuptable[n].size;
			break;
		}
	}
	/* Exit with a fatal error if variable isn't found */
	if (n == LOOKUPTABLE_SIZE) {
		printf("Fatal Error:\n");
		die("Variable name not found ", name);
	}

	/* We are all good and have valid values in valuetype struct */

	base = (char *)prm;
	newbase = (base + valuetype.offset);
	switch (valuetype.etype) {
	case CHAR:
		sprintf(newbase, "%s", value);
		break;

	case INT:
		intvalue = strtol(value, NULL, 10);
		/* An int cannot be copied with an assignment to the struct member here */
		memcpy(newbase, &intvalue, sizeof(int));
		break;

	case DOUBLE:
		doublevalue = atof(value);
		/* A double cannot be copied with an assignment to the struct member here */
		memcpy(newbase, &doublevalue, sizeof(double));
		break;
	}

	return 0;
}

char *get_PRM(char *prmfile, char *valuename) {
	struct PRMPRIMITIVE valuetype;
	int n = 0;
	char *szReturn = NULL;
	init_lookup_table();

	for (n = 0; n < LOOKUPTABLE_SIZE; n++) {
		if (m_strcmp(lookuptable[n].name, valuename) == 0 || m_strcmp(lookuptable[n].alias, valuename) == 0) {
			valuetype.name = strdup(lookuptable[n].name);
			valuetype.etype = lookuptable[n].etype;
			valuetype.offset = lookuptable[n].offset;
			valuetype.size = lookuptable[n].size;
			valuetype.alias = strdup(lookuptable[n].alias);
			break;
		}
	}

	/* Exit with a fatal error if variable isn't found */
	if (n == LOOKUPTABLE_SIZE) {
		printf("Fatal Error:\n");
		die("Variable name not found ", valuename);
	}

	/* We are all good and have a valid valuesname */

	szReturn = get_PRM_sub(prmfile, valuetype.name);

	if (szReturn == NULL)
		szReturn = get_PRM_sub(prmfile, valuetype.alias);

	free(valuetype.name);
	free(valuetype.alias);

	return szReturn;
}

char *get_PRM_sub(char *szFile, char *szValueName) {
	FILE *fp;
	char *line;
	int maxline = 0;
	LIST *current, *head;

	char *valuename;
	char *value;
	char *szRetVal = NULL;

	maxline = (findlongest(szFile, 1) + 1) * sizeof(char);
	line = (char *)malloc(maxline);
	valuename = (char *)malloc(maxline);
	value = (char *)malloc(maxline);

	head = current = NULL;
	fp = fopen(szFile, "r");

	while (fgets(line, maxline, fp)) {
		if (strchr(line, '=') != 0 && strchr(line, '#') == 0) {

			LIST *node = malloc(sizeof(LIST));

			char *token = strtok(line, "=");
			node->string = strdup(trimwhitespace(token));
			token = strtok(NULL, "=");
			node->value = strdup(trimwhitespace(token));
			node->next = NULL;

			if (head == NULL) {
				current = head = node;
			}
			else {
				current = current->next = node;
			}
		}
	}
	fclose(fp);
	for (current = head; current; current = current->next) {

		if (strcmp(trimwhitespace(szValueName), current->string) == 0) {
			free(line);
			free(valuename);
			free(value);
			if (strlen(trimwhitespace(current->value)))
				szRetVal = strdup(trimwhitespace(current->value));
			else
				szRetVal = NULL;
		}
	}
	// need free for each node
	for (current = head; current; current = current->next)
		free(current);

	return szRetVal;
}

int findlongest(char *szFile, int func) {

	FILE *file = fopen(szFile, "rb");
	if (file == NULL) {
		fprintf(stderr, "File open failed!\n");
		return EXIT_FAILURE;
	}

	size_t buffer_size = 8;
	char *line = (char *)malloc(buffer_size);
	char *longest = (char *)malloc(buffer_size);
	size_t longest_on = 0;
	size_t longest_length = 0;
	size_t num_lines = 0;

	if (line == NULL || longest == NULL) {
		fprintf(stderr, "Memory allocation failed!\n");
		return EXIT_FAILURE;
	}

	long pos_before = ftell(file);
	for (num_lines = 1; fgets(line, buffer_size, file); ++num_lines) {
		long pos_after = ftell(file);
		size_t line_length = pos_after - pos_before - 1;
		pos_before = pos_after - 1;
		while (line_length == buffer_size - 1 && line[line_length - 1] != '\n') {
			line = (char *)realloc(line, buffer_size * 2);
			longest = (char *)realloc(longest, buffer_size * 2);
			if (line == NULL || longest == NULL) {
				fprintf(stderr, "Memory allocation failed!\n");
				return EXIT_FAILURE;
			}
			if (fgets(&line[buffer_size - 1], buffer_size + 1, file) == NULL) {
				fprintf(stderr, "File read failed!\n");
				fclose(file);
				free(line);
				free(longest);
				return EXIT_FAILURE;
			}
			pos_after = ftell(file);
			line_length += pos_after - pos_before - 1;
			pos_before = pos_after - 1;
			buffer_size *= 2;
		}
		if (line_length > longest_length) {
			char *tmp = longest;
			longest = line;
			line = tmp;
			longest_on = num_lines;
			longest_length = line_length;
		}
	}

	if (longest_on) {
		switch (func) {
		case 1:
			return longest_length;
			break;

		case 2:
			return longest_on;
			break;

		default:
			return longest_length;
		}

		// printf("File had total of %lu lines.\n", num_lines);
		// printf("Longest line found on line %lu with %lu characters.\n",
		// longest_on, longest_length); printf("The line was: "); fwrite(longest, 1,
		// longest_length, stdout);
	}
	else {
		return EXIT_FAILURE;
		// fprintf(stderr, "The file had no readable lines!\n");
	}

	fclose(file);
	free(line);
	free(longest);
	return EXIT_SUCCESS;
}
