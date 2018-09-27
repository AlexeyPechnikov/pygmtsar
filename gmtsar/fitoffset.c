#define _GNU_SOURCE
#include "update_PRM.h"
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

struct valuelist {
	double value1;
	double value2;
	double value3;
	double value4;
	double value5;
	struct valuelist *next;
};

#ifndef MAX_PATH
#define MAX_PATH 1024
#endif

typedef struct valuelist VALUELIST;

VALUELIST *vhead;
int NPTS0;

int load_freq_xcorr(char *freq_xcorr_file);
char **str_split(char *a_str, const char a_delim);
void replace_multi_space_with_single_space(char *str);
int run_gmt_trend(char *file, int n_model_params, double *coefficient1, double *coefficient2, double *coefficient3);
int run_gmt_gmtinfo(char *file, double *coefficient1, double *coefficient2, double *coefficient3, double *coefficient4);
int file_exists(char *filename);

int main(int argc, char **argv) {
	double SNR = 20.0;
	int rng_p = 0;
	int azi_p = 0;
	VALUELIST *vcurrent;
	char s_r_xyz[MAX_PATH] = "";
	char s_a_xyz[MAX_PATH] = "";
	FILE *r_xyz = NULL;
	FILE *a_xyz = NULL;
	int n_r_xyz = -1;
	int n_a_xyz = -1;
	int NPTS = 0;
	double n_coeff1_r = 0.0;
	double n_coeff2_r = 0.0;
	double n_coeff3_r = 0.0;

	double n_coeff1_a = 0.0;
	double n_coeff2_a = 0.0;
	double n_coeff3_a = 0.0;

	double n_range1 = 0.0;
	double n_range2 = 0.0;
	double n_range3 = 0.0;
	double n_range4 = 0.0;

	double rshift = 0.0;
	double ashift = 0.0;

	int int_rshift = 0;
	int int_ashift = 0;
	double sub_int_r = 0.0;
	double sub_int_a = 0.0;
	double stretch_r = 0.0;
	double a_stretch_r = 0.0;
	double stretch_a = 0.0;
	double a_stretch_a = 0.0;

	char tmpstring[MAX_PATH];

	if (argc > 6 || argc == 1 || argv[1] == NULL || argv[2] == NULL || file_exists(argv[4]) == 0) {
		printf("Usage: fitoffset npar_rng npar_azi xcorr.dat prmfile.PRM [SNR]\n");
		return 1;
	}

	if (argc == 5 && atof(argv[4]) != 0.0)
		SNR = atof(argv[4]);

	rng_p = atoi(argv[1]);
	azi_p = atoi(argv[2]);

	load_freq_xcorr(argv[3]);

	/* Open a temporary file for r.xyz */
	strlcpy(s_r_xyz, "/tmp/r.xyz.XXXXXX", sizeof s_r_xyz);
	if ((n_r_xyz = mkstemp(s_r_xyz)) == -1 || (r_xyz = fdopen(n_r_xyz, "w+")) == NULL) {
		if (n_r_xyz != -1) {
			unlink(s_r_xyz);
			close(n_r_xyz);
		}
		fprintf(stderr, "%s: %s\n", s_r_xyz, strerror(errno));
		return (0);
	}
	/* Open a temporary file for a.xyz */
	strlcpy(s_a_xyz, "/tmp/a.xyz.XXXXXX", sizeof s_a_xyz);
	if ((n_a_xyz = mkstemp(s_a_xyz)) == -1 || (a_xyz = fdopen(n_a_xyz, "w+")) == NULL) {
		if (n_a_xyz != -1) {
			unlink(s_a_xyz);
			close(n_a_xyz);
		}
		fprintf(stderr, "%s: %s\n", s_a_xyz, strerror(errno));
		return (0);
	}

	for (vcurrent = vhead; vcurrent; vcurrent = vcurrent->next) {
		if (vcurrent->value5 > SNR) {
			/* Write the values to r.xyz */
			fprintf(r_xyz, "%.6f %.6f %.6f\n", vcurrent->value1, vcurrent->value3, vcurrent->value2);
			/* Increment NPTS by 1 for every line written */
			NPTS++;
		}
	}

	for (vcurrent = vhead; vcurrent; vcurrent = vcurrent->next) {
		if (vcurrent->value5 > SNR) {
			/* Write the values to a.xyz */
			fprintf(a_xyz, "%.6f %.6f %.6f\n", vcurrent->value1, vcurrent->value3, vcurrent->value4);
		}
	}

	if (NPTS < 8) {
		if (r_xyz)
			fclose(r_xyz);

		if (a_xyz)
			fclose(a_xyz);

		/* Remove the temporary files */
		unlink(s_r_xyz);
		unlink(s_a_xyz);

		printf("FAILED - not enough points to estimate parameters\ntry using a "
		       "lower SNR\nNPTS0 = %d NPTS = %d\n\n",
		       NPTS0, NPTS);
		return 1;
	}

	fflush(r_xyz);
	fflush(a_xyz);

	/* Close any open files */
	if (r_xyz)
		fclose(r_xyz);

	if (a_xyz)
		fclose(a_xyz);

	/* Run gmt to get coefficients */
	run_gmt_trend(s_r_xyz, rng_p, &n_coeff1_r, &n_coeff2_r, &n_coeff3_r);
	run_gmt_trend(s_a_xyz, azi_p, &n_coeff1_a, &n_coeff2_a, &n_coeff3_a);

	/* run gmt gmtinfo to get range and azimuth coefficients */
	run_gmt_gmtinfo(s_r_xyz, &n_range1, &n_range2, &n_range3, &n_range4);

	/* Calculate range and azimuth shifts */
	rshift = n_coeff1_r - n_coeff2_r * (n_range2 + n_range1) / (n_range2 - n_range1) -
	         n_coeff3_r * (n_range4 + n_range3) / (n_range4 - n_range3);
	ashift = n_coeff1_a - n_coeff2_a * (n_range2 + n_range1) / (n_range2 - n_range1) -
	         n_coeff3_a * (n_range4 + n_range3) / (n_range4 - n_range3);

	/* Calculate range and azimuth stretches */
	stretch_r = n_coeff2_r * 2.0 / (n_range2 - n_range1);
	a_stretch_r = n_coeff3_r * 2.0 / (n_range4 - n_range3);

	stretch_a = n_coeff2_a * 2.0 / (n_range2 - n_range1);
	a_stretch_a = n_coeff3_a * 2.0 / (n_range4 - n_range3);

	if (rshift >= 0) {
		int_rshift = trunc(rshift);
		sub_int_r = (rshift - (double)trunc(rshift));
	}
	else {
		int_rshift = floor(rshift);
		sub_int_r = 1 + (rshift - (double)trunc(rshift));
	}

	if (ashift >= 0) {
		int_ashift = trunc(ashift);
		sub_int_a = (ashift - (double)trunc(ashift));
	}
	else {
		int_ashift = floor(ashift);
		sub_int_a = 1 + (ashift - (double)trunc(ashift));
	}

#ifdef DEBUG
	printf("s_r_xyz %s\n", s_r_xyz);
	printf("%d lines\n", NPTS0);
	printf("SNR = %.3f\n", SNR);
	printf("n_coeff1_r = %.6f n_coeff2_r = %.6f n_coeff3_r = %.6f\n", n_coeff1_r, n_coeff2_r, n_coeff3_r);
	printf("n_coeff1_r = %.6f n_coeff2_a = %.6f n_coeff3_a = %.6f\n", n_coeff1_a, n_coeff2_a, n_coeff3_a);
	printf("n_range1 = %.0f n_range2 = %.0f n_range3 = %.0f n_range4 = %.0f\n\n", n_range1, n_range2, n_range3, n_range4);
	printf("Range Shift = %.4f\n", rshift);
	printf("Azimuth Shift = %.4f\n", ashift);

	printf("rshift = %d\n", int_rshift);
	printf("sub_int_r = %.4f\n", sub_int_r);
	printf("stretch_r = %.9f\n", stretch_r);
	printf("a_stretch_r = %.9f\n", a_stretch_r);

	printf("ashift = %d\n", int_ashift);
	printf("sub_int_a = %.4f\n", sub_int_a);
	printf("stretch_a = %.9f\n", stretch_a);
	printf("a_stretch_a = %.9f\n", a_stretch_a);
#endif

	/* Update rshift and associated params in PRM*/
	sprintf(tmpstring, "%d", int_rshift);
	update_PRM_sub(argv[4], "rshift", tmpstring);
	sprintf(tmpstring, "%.4f", sub_int_r);
	update_PRM_sub(argv[4], "sub_int_r", tmpstring);
	sprintf(tmpstring, "%.9f", stretch_r);
	update_PRM_sub(argv[4], "stretch_r", tmpstring);
	sprintf(tmpstring, "%.9f", a_stretch_r);
	update_PRM_sub(argv[4], "a_stretch_r", tmpstring);

	/* Update ashift and associated params in PRM*/
	sprintf(tmpstring, "%d", int_ashift);
	update_PRM_sub(argv[4], "ashift", tmpstring);
	sprintf(tmpstring, "%.4f", sub_int_a);
	update_PRM_sub(argv[4], "sub_int_a", tmpstring);
	sprintf(tmpstring, "%.9f", stretch_a);
	update_PRM_sub(argv[4], "stretch_a", tmpstring);
	sprintf(tmpstring, "%.9f", a_stretch_a);
	update_PRM_sub(argv[4], "a_stretch_a", tmpstring);

	/* Remove the temporary files */
	unlink(s_r_xyz);
	unlink(s_a_xyz);
	printf("fitoffset:\nUpdated %s with offset parameters calculated from %s\n", argv[4], argv[3]);
	return 0;
}

int load_freq_xcorr(char *freq_xcorr_file) {
	FILE *fp;
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	VALUELIST *vcurrent;
	char **tokens;

	vhead = vcurrent = NULL;
	NPTS0 = 0;
	fp = fopen(freq_xcorr_file, "r");

	if (fp == NULL)
		return EXIT_FAILURE;

	while ((read = getline(&line, &len, fp)) != -1) {
		NPTS0++;
		replace_multi_space_with_single_space(line);

		VALUELIST *vnode = malloc(sizeof(VALUELIST));
		tokens = str_split(trimwhitespace(line), ' ');

		if (tokens) {
			for (int i = 0; *(tokens + i); i++) {
				switch (i) {
				case 0:
					vnode->value1 = atof(*(tokens + i));
					break;

				case 1:
					vnode->value2 = atof(*(tokens + i));

					break;

				case 2:
					vnode->value3 = atof(*(tokens + i));

					break;

				case 3:
					vnode->value4 = atof(*(tokens + i));

					break;

				case 4:
					vnode->value5 = atof(*(tokens + i));
					break;
				}

				free(*(tokens + i));
			}
			free(tokens);
		}

		if (vhead == NULL) {
			vcurrent = vhead = vnode;
		}
		else {
			vcurrent = vcurrent->next = vnode;
		}
	}

	// need free for each node
	if (fp)
		fclose(fp);
	return EXIT_SUCCESS;
}

char **str_split(char *a_str, const char a_delim) {
	char **result = 0;
	size_t count = 0;
	char *tmp = a_str;
	char *last_comma = 0;
	char delim[2];
	delim[0] = a_delim;
	delim[1] = 0;

	/* Count how many elements will be extracted. */
	while (*tmp) {
		if (a_delim == *tmp) {
			count++;
			last_comma = tmp;
		}
		tmp++;
	}

	/* Add space for trailing token. */
	count += last_comma < (a_str + strlen(a_str) - 1);

	/* Add space for terminating null string so caller
	   knows where the list of returned strings ends. */
	count++;

	result = malloc(sizeof(char *) * count);

	if (result) {
		size_t idx = 0;
		char *token = strtok(a_str, delim);

		while (token) {
			assert(idx < count);
			*(result + idx++) = strdup(token);
			token = strtok(0, delim);
		}
		assert(idx == count - 1);
		*(result + idx) = 0;
	}

	return result;
}

void replace_multi_space_with_single_space(char *str) {
	char *dest = str; /* Destination to copy to */

	/* While we're not at the end of the string, loop... */
	while (*str != '\0') {
		/* Loop while the current character is a space, AND the next
		 * character is a space
		 */
		while (*str == ' ' && *(str + 1) == ' ')
			str++; /* Just skip to next character */

		/* Copy from the "source" string to the "destination" string,
		 * while advancing to the next character in both
		 */
		*dest++ = *str++;
	}

	/* Make sure the string is properly terminated */
	*dest = '\0';
}

/* Run gmt trend2d and collect the parameters */
int run_gmt_trend(char *file, int n_model_params, double *coefficient1, double *coefficient2, double *coefficient3) {

	char cmd_line[MAX_PATH] = "gmt trend2d /tmp/r.xyz.ZBbxXu -Fxyz -N\"2\"r -V >  /dev/null";
	// char buffer[1000];
	FILE *pipe;
	// int len;
	size_t glen;
	ssize_t read;
	// FILE *fp;
	char *line = NULL;
	char strfound[] = "trend2d: Model Coefficients:";
	char **tokens;

	sprintf(cmd_line, "gmt trend2d %s -Fxyz -N\"%d\"r -V  2>&1", file, n_model_params);

	pipe = popen(cmd_line, "r");

	if (NULL == pipe) {
		perror("pipe");
		return 1;
	}

	while ((read = getline(&line, &glen, pipe)) != -1) {
		if (strncmp(strfound, line, 28) == 0) {
			for (int i = 0; i < strlen(line); i++) {
				if (line[i] == '\t')
					line[i] = ' ';
			}
			tokens = str_split(trimwhitespace(line), ' ');
			if (tokens) {
				for (int i = 0; *(tokens + i); i++) {
					switch (i) {
					case 3:
						*coefficient1 = atof(*(tokens + i));
						break;
					case 4:
						*coefficient2 = atof(*(tokens + i));
						break;
					case 5:
						*coefficient3 = atof(*(tokens + i));
						break;
					}

					free(*(tokens + i));
				}
				free(tokens);
			}
		}
	}

	pclose(pipe);

	return 0;
}

/* Run gmt gmtinfo and collect the parameters */
int run_gmt_gmtinfo(char *file, double *coefficient1, double *coefficient2, double *coefficient3, double *coefficient4) {

	char cmd_line[MAX_PATH] = "";
	// char buffer[1000];
	FILE *pipe;
	// int len;
	size_t glen;
	ssize_t read;
	// FILE *fp;
	char *line = NULL;

	char **tokens;

	sprintf(cmd_line, "gmt gmtinfo  %s -C  2>&1", file);

	pipe = popen(cmd_line, "r");

	if (NULL == pipe) {
		perror("pipe");
		return 1;
	}

	while ((read = getline(&line, &glen, pipe)) != -1) {

		for (int i = 0; i < strlen(line); i++) {
			if (line[i] == '\t')
				line[i] = ' ';
		}
		tokens = str_split(trimwhitespace(line), ' ');
		if (tokens) {
			for (int i = 0; *(tokens + i); i++) {
				switch (i) {
				case 0:
					*coefficient1 = atof(*(tokens + i));
					break;
				case 1:
					*coefficient2 = atof(*(tokens + i));
					break;
				case 2:
					*coefficient3 = atof(*(tokens + i));
					break;
				case 3:
					*coefficient4 = atof(*(tokens + i));
					break;
				}

				free(*(tokens + i));
			}
			free(tokens);
		}
	}

	pclose(pipe);

	return 0;
}

int file_exists(char *filename) {
	FILE *file = NULL;
	file = fopen(filename, "r");
	if (file != NULL) {
		fclose(file);
		return 1;
	}
	return 0;
}
