/***************************************************************************
 * Creator:  Xiaohua(Eric) XU & David Sandwell                             *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  04/26/2015                                                    *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * DATE                                                                    *
 *                                                                         *
 ***************************************************************************/

#include "PRM.h"
#include "lib_functions.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double str2double(char *);
int strasign(char *, char *, int, int);
int strlocate(char *, int, int);

char *USAGE = "\n\nUsage: make_gaussian_filter name_of_PRM_file RNG_DEC "
              "AZI_DEC WAVELENGTH(m)\n"
              "\nExample: make_gaussian_filter "
              "IMG-HH-ALPSRP211830620-H1.0__A.PRM 2 4 200\n"
              "\nOutput: gauss_200\n";

int main(int argc, char **argv) {

	FILE *fid;
	int rng_dec, azi_dec, w, n_azi, n_rng;
	int i, j, idec, jdec;
	struct PRM prm;
	double azi_px_size, rng_px_size, x, y;
	double c_speed = 299792458.0;
	double sig_azi, sig_rng, a, rng, cost, cosa;
	double g[50][50];
	char out_name[100] = "gauss_";

	if (argc < 4)
		die(USAGE, "");
	rng_dec = (int)str2double(argv[2]);
	azi_dec = (int)str2double(argv[3]);
	w = (int)str2double(argv[4]);
	if (rng_dec < 1 || rng_dec > 4)
		die("Incorrect range decimation factor: \n", argv[2]);
	if (azi_dec < 1 || azi_dec > 4)
		die("Incorrect azimuth decimation factor: \n", argv[3]);
	if (w < 2 || w > 10000)
		die("Incorrect wavelength: \n", argv[4]);

	// get the prm
	get_prm(&prm, argv[1]);

	// compute the range and azimuth pixel size
	azi_px_size = prm.vel / sqrt(1 + prm.ht / prm.RE) / prm.prf; // real_vel/prf
	rng_px_size = c_speed / prm.fs / 2;

	// compute the cosine of the looking angle and the surface deviate angle
	a = prm.ht + prm.RE;
	prm.far_range = prm.near_range + rng_px_size * (double)prm.num_rng_bins;
	rng = (prm.near_range + prm.far_range) / 2;
	cost = (pow(a, 2.0) + pow(rng, 2.0) - pow(prm.RE, 2.0)) / 2 / a / rng;
	cosa = (pow(a, 2.0) + pow(prm.RE, 2.0) - pow(rng, 2.0)) / 2 / a / prm.RE;
	// fprintf(stderr,"cosa = %.9f, cost = %.9f\n", cosa ,cost);

	// compute the ground range pixel size
	rng_px_size = rng_px_size / sin(acos(cost) + acos(cosa));

	azi_px_size = azi_px_size * azi_dec;
	rng_px_size = rng_px_size * rng_dec;

	sig_azi = w / 5.3 / azi_px_size;
	sig_rng = w / 5.3 / rng_px_size;
	idec = floor(sig_azi / 2.);
	if (idec < 1)
		idec = 1;
	jdec = floor(sig_rng / 2.);
	if (jdec < 1)
		jdec = 1;
	fprintf(stdout, " %d %d \n", idec, jdec);

	n_azi = (int)(sig_azi * 4);
	n_rng = (int)(sig_rng * 4);

	if (n_azi % 2 == 0)
		n_azi = n_azi + 1;
	if (n_rng % 2 == 0)
		n_rng = n_rng + 1;

	strcat(out_name, argv[4]);

	if ((fid = fopen(out_name, "w")) == NULL)
		die("Couldn't open file: \n", out_name);
	fprintf(fid, "%d %d\n", n_rng, n_azi);
	for (i = 0; i < n_azi; i++) {
		for (j = 0; j < n_rng; j++) {
			x = (-(n_rng - 1) / 2 + j);
			y = (-(n_azi - 1) / 2 + i);
			g[i][j] = exp(-(x * x / sig_rng / sig_rng + y * y / sig_azi / sig_azi) / 2.0);
			fprintf(fid, "\t%.16e", g[i][j]);
		}
		fprintf(fid, "\n");
	}
	fclose(fid);
}

int strasign(char *str_out, char *str, int n1, int n2) {
	// asign n1-n2 of str to str_out
	int i;
	if (n1 > n2 || n2 > 199) {
		return (-1);
	}
	for (i = n1; i <= n2; i++) {
		str_out[i - n1] = str[i];
	}
	str_out[n2 - n1 + 1] = '\0';
	return (1);
}

int strlocate(char *str, int c, int n) {
	// locate the n-th c in str
	int i, j = 0;
	for (i = 0; i < strlen(str); i++) {
		if (str[i] == (char)c) {
			j++;
			if (j == n) {
				return (i);
			}
		}
	}
	return (-1);
}

double str2double(char *str) {
	int i, n, m;
	double value = 0.0, value1 = 0.0, value2 = 0.0, sgn = 1.0;
	char tmp1[100], tmp2[100], tmp[100], str_tmp[100];

	strasign(str_tmp, str, 0, strlen(str));

	// decide the sign
	if (str_tmp[0] == '-' || str_tmp[0] == '+') {
		if (str_tmp[0] == '-') {
			sgn = -1.0;
		}
		strasign(tmp, str_tmp, 1, strlen(str_tmp));
		strasign(str_tmp, tmp, 0, strlen(tmp));
	}

	// decide where it is sci form
	if (strlocate(str_tmp, 'e', 1) != -1 || strlocate(str_tmp, 'E', 1) != -1) {
		n = strlocate(str_tmp, 'e', 1);
		if (n == -1) {
			n = strlocate(str_tmp, 'E', 1);
		}
		strasign(tmp2, str_tmp, n + 1, strlen(str)); // exponential part
		strasign(tmp1, str_tmp, 0, n - 1);           // digits part
	}
	else {
		strasign(tmp1, str_tmp, 0, strlen(str_tmp)); // digits part
	}
	// decide whether it has fraction
	n = strlocate(tmp1, '.', 1);
	if (n != -1) {
		strasign(tmp, tmp1, 0, n - 1);
		m = strlen(tmp);
		for (i = 0; i < m; i++) {
			value1 = value1 + (double)((int)tmp[i] - 48) * pow(10.0, (double)(m - i - 1));
		}
		m = strlen(tmp1);
		strasign(tmp, tmp1, n + 1, m);
		m = strlen(tmp);
		for (i = 0; i < m; i++) {
			value2 = value2 + (double)((int)tmp[i] - 48) * pow(10.0, (double)(-i - 1));
		}
		// puts(tmp2);
		// fprintf(stderr,"%.12f    %.12f %.12f\n",value1,value2,str2double(tmp2));
		value = value1 + value2;
	}
	else {
		m = strlen(tmp1);
		for (i = 0; i < m; i++) {
			value = value + (double)((int)tmp1[i] - 48) * pow(10.0, (double)(m - i - 1));
		}
	}

	if (strlocate(str_tmp, 'e', 1) != -1 || strlocate(str_tmp, 'E', 1) != -1) {
		value = value * pow(10.0, str2double(tmp2));
	}

	return (value * sgn);
}
