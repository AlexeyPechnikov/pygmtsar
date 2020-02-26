/***************************************************************************
 * Creator:  David Sandwell                                                *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  09/23/2015                                                    *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * DATE   :  EX made the code capable of reading long vectors of orbits    *
 *                                                                         *
 ***************************************************************************/

#include "PRM.h"
#include "lib_defs.h"
#include "lib_functions.h"
#include "stateV.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// int pop_prm(struct PRM *, tree *, char *);
int pop_led_pre(tree *, state_vector *, double, double);
int write_orb(state_vector *sv, FILE *fp, int);
double yds2ydf(state_vector);

char *USAGE = "\n\nUsage: ext_orb_s1a.c name_of_prm name_of_orb_xml name_output\n"
              "\nExample: ext_orb_s1a.c S1A20140807.PRM "
              "S1A_OPER_AUX_POEORB_OPOD_20150615T155109_V20150525T225944_20150527T005944."
              "EOF output\n"
              "\nOutput: output.LED\n";

int main(int argc, char **argv) {

	FILE *XML_FILE, *INPUT_PRM, *OUTPUT_LED;
	char tmp_str[200];
	struct PRM prm;
	tree *xml_tree;
	state_vector *sv;
	int ch, n = 0, nc = 0, nlmx = 0;
	double t1, t2;

	if (argc < 3)
		die(USAGE, "");

	// find the number of lines and the maximum line length of the xml file
	if ((XML_FILE = fopen(argv[2], "r")) == NULL)
		die("Couldn't open xml file: \n", argv[2]);
	while (EOF != (ch = fgetc(XML_FILE))) {
		++nc;
		if (ch == '\n') {
			++n;
			if (nc > nlmx)
				nlmx = nc;
			nc = 0;
		}
	}
	xml_tree = (struct tree *)malloc(5 * n * sizeof(struct tree));
	sv = (struct state_vector *)malloc(n / 2 * sizeof(struct state_vector));
	fclose(XML_FILE);

	if ((XML_FILE = fopen(argv[2], "r")) == NULL)
		die("Couldn't open xml file: \n", argv[2]);
	get_tree(XML_FILE, xml_tree, 1);
	// show_tree(xml_tree,0,0);
	fclose(XML_FILE);

	// get the prm
	if ((INPUT_PRM = fopen(argv[1], "r")) == NULL)
		die("Couldn't open xml file: \n", argv[1]);
	get_sio_struct(INPUT_PRM, &prm);
	// fprintf(stderr,"%.12f  %.12f\n",prm.SC_clock_start,prm.SC_clock_stop);
    // revising t1 and t2 in case only one burst is used.
	t1 = prm.SC_clock_start - 10 * (prm.SC_clock_stop - prm.SC_clock_start);
	t2 = prm.SC_clock_stop + 10 * (prm.SC_clock_stop - prm.SC_clock_start);
    if (fabs(t2-t1) < 1/86400.0*100) {
        t1 = prm.SC_clock_start - 20 * (prm.SC_clock_stop - prm.SC_clock_start);
        t2 = prm.SC_clock_stop + 20 * (prm.SC_clock_stop - prm.SC_clock_start);
    }

	// generate the LED file
	n = pop_led_pre(xml_tree, sv, t1, t2);

	strcpy(tmp_str, argv[3]);
	strcat(tmp_str, ".LED");
	if ((OUTPUT_LED = fopen(tmp_str, "w")) == NULL)
		die("Couldn't open led file: \n", tmp_str);
	// tmp_str = search_tree()
	write_orb(sv, OUTPUT_LED, n);
	fclose(OUTPUT_LED);
	free(xml_tree);
	free(sv);
}

int write_orb(state_vector *sv, FILE *fp, int n) {
	int i;
	double dt;

	if (n <= 0) {
		fprintf(stderr, "NO orbit coverage in the selected file...\n");
		return (-1);
	}
	else {
		printf("Writing %d lines of precise orbit for the LED file...\n", n);
	}

	dt = round((sv[1].sec) * 1000.0) / 1000.0 - round((sv[0].sec) * 1000.0) / 1000.0;
	// printf("%f,%f\n",sv[1].sec,sv[0].sec);
	if (n <= 1)
		return (-1);
	fprintf(fp, "%d %d %d %.6lf %.3lf \n", n, sv[0].yr, sv[0].jd, sv[0].sec, dt);
	for (i = 0; i < n; i++) {
		fprintf(fp, "%d %d %.6lf %.6lf %.6lf %.6lf %.8lf %.8lf %.8lf \n", sv[i].yr, sv[i].jd, sv[i].sec, sv[i].x, sv[i].y,
		        sv[i].z, sv[i].vx, sv[i].vy, sv[i].vz);
	}
	return (1);
}

int pop_led_pre(tree *xml_tree, state_vector *sv, double t1, double t2) {
	int i, count, num = 0;
	char tmp_c[200], tmp_y[200];
	double tmp_d, tmp_t;
	state_vector tmp_sv;

	search_tree(xml_tree, "/Earth_Explorer_File/Data_Block/List_of_OSVs/", tmp_c, 3, 0, 1);
	count = (int)str2double(tmp_c);
	// printf("Reading %d lines from precise orbit...\n",count);
	for (i = 0; i < count; i++) {
		search_tree(xml_tree, "/Earth_Explorer_File/Data_Block/List_of_OSVs/OSV/UTC/", tmp_c, 2, 4, i + 1);
		tmp_d = str2double(tmp_c);
		search_tree(xml_tree, "/Earth_Explorer_File/Data_Block/List_of_OSVs/OSV/UTC/", tmp_c, 1, 4, i + 1);
		strasign(tmp_y, tmp_c, 4, 7);
		tmp_sv.yr = (int)(str2double(tmp_y));
		tmp_sv.jd = (int)(tmp_d - trunc(tmp_d / 1000.0) * 1000.0);
		tmp_sv.sec = (tmp_d - trunc(tmp_d)) * 86400;
		tmp_t = yds2ydf(tmp_sv);
		if (tmp_t < t1 || tmp_t > t2) {
			continue;
		}
		else {
			sv[num].yr = tmp_sv.yr;
			sv[num].jd = tmp_sv.jd;
			sv[num].sec = tmp_sv.sec;
			search_tree(xml_tree, "/Earth_Explorer_File/Data_Block/List_of_OSVs/OSV/X/", tmp_c, 1, 4, i + 1);
			sv[num].x = str2double(tmp_c);
			search_tree(xml_tree, "/Earth_Explorer_File/Data_Block/List_of_OSVs/OSV/Y/", tmp_c, 1, 4, i + 1);
			sv[num].y = str2double(tmp_c);
			search_tree(xml_tree, "/Earth_Explorer_File/Data_Block/List_of_OSVs/OSV/Z/", tmp_c, 1, 4, i + 1);
			sv[num].z = str2double(tmp_c);
			search_tree(xml_tree, "/Earth_Explorer_File/Data_Block/List_of_OSVs/OSV/VX/", tmp_c, 1, 4, i + 1);
			sv[num].vx = str2double(tmp_c);
			search_tree(xml_tree, "/Earth_Explorer_File/Data_Block/List_of_OSVs/OSV/VY/", tmp_c, 1, 4, i + 1);
			sv[num].vy = str2double(tmp_c);
			search_tree(xml_tree, "/Earth_Explorer_File/Data_Block/List_of_OSVs/OSV/VZ/", tmp_c, 1, 4, i + 1);
			sv[num].vz = str2double(tmp_c);
			num++;
		}
	}
	return (num);
}

double yds2ydf(state_vector s) {
	// convert year, day, sec to year, day, fraction of day
	double t = 0;
	t = t + (double)s.yr * 1000;
	t = t + (double)s.jd;
	t = t + (double)s.sec / 86400.;
	return (t);
}
