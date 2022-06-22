/***************************************************************************
 * Creator:  Xiaohua(Eric) XU                                              *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  10/25/2016                                                    *
 ***************************************************************************/
/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 ***************************************************************************/

#include "lib_defs.h"
#include "lib_functions.h"
#include "tiffio.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char *USAGE = "\nUsage: assemble_tops azi_1 azi_2 name_stem1 name_stem2 ... output_stem\n"
              "\nExample: assemble_tops 1685 9732 "
              "s1a-iw1-slc-vv-20150706t135900-20150706t135925-006691-008f28-001 "
              "s1a-iw1-slc-vv-20150706t135925-20150706t135950-006691-008f28-001 "
              "s1a-iw1-slc-vv-20150706t135900-20150706t135950-006691-008f28-001\n"
              "\nOutput:s1a-iw1-slc-vv-20150706t135900-20150706t135950-006691-008f28-001."
              "xml "
              "s1a-iw1-slc-vv-20150706t135900-20150706t135950-006691-008f28-001.tiff\n\n"
              "\nNote: output files are bursts that covers area between azi_1 and azi_2, "
              "set them to 0s to output all bursts\n";

int edit_tree(int, int, struct tree **, int, int, int *, int *);
int assemble_tifs(TIFF **, TIFF *, int, int, int, int);

int main(int argc, char **argv) {

	FILE *XML_FILE = NULL, *XML_OUTPUT = NULL;
	TIFF **tif, *tif_out;
	int na, nc, ncmx, nlmx, ii, nfiles, ch;
	char tmp_str[2000];
	struct tree **xml_tree;
	float azi_1, azi_2;
	int nb_start, nb_end, lpb;

	if (argc < 5) {
		fprintf(stderr, "%s", USAGE);
		die("Error: Incorrect input...", "\n");
	}

	azi_1 = atof(argv[1]);
	azi_2 = atof(argv[2]);
	if (azi_1 >= azi_2 && !(azi_1 == 0 && azi_2 == 0)) {
		die("azi_2 has to be larger than azi_1", "");
	}

	nfiles = argc - 4;
	null_MEM_STR();
	xml_tree = (struct tree **)malloc(nfiles * sizeof(struct tree *));
	TIFFSetWarningHandler(NULL);

	// figure out the sizes
	nlmx = 0;
	for (ii = 0; ii < nfiles; ii++) {
		na = 0;
		nc = 0;
		ncmx = 0;

		strcpy(tmp_str, argv[ii + 3]);
		strcat(tmp_str, ".xml");
		if ((XML_FILE = fopen(tmp_str, "r")) == NULL)
			die("Couldn't open xml file: \n", tmp_str);
		while (EOF != (ch = fgetc(XML_FILE))) {
			++nc;
			if (ch == '\n') {
				++na;
				if (nc > ncmx)
					ncmx = nc;
				nc = 0;
			}
			if (na > nlmx)
				nlmx = na;
		}
		fclose(XML_FILE);
	}
	// read in xmls and duplicate them in xml_tree[0]
	printf("Files to be assembled is %d (x %d lines max)\n", nfiles, nlmx);
	xml_tree[0] = (struct tree *)malloc(nlmx * nfiles * 5 * sizeof(struct tree));
	for (ii = 1; ii < nfiles; ii++)
		xml_tree[ii] = &xml_tree[0][nlmx * 5 * ii];
	for (ii = 0; ii < nfiles; ii++) {
		strcpy(tmp_str, argv[ii + 3]);
		strcat(tmp_str, ".xml");
		if ((XML_FILE = fopen(tmp_str, "r")) == NULL)
			die("Couldn't open xml file: \n", tmp_str);
		get_tree(XML_FILE, xml_tree[ii], 1);
		fclose(XML_FILE);
	}
	for (ii = 1; ii < nfiles; ii++) {
		xml_tree[ii] = (struct tree *)malloc(nlmx * 5 * sizeof(struct tree));
		strcpy(tmp_str, argv[ii + 3]);
		strcat(tmp_str, ".xml");
		if ((XML_FILE = fopen(tmp_str, "r")) == NULL)
			die("Couldn't open xml file: \n", tmp_str);
		get_tree(XML_FILE, xml_tree[ii], 1);
		fclose(XML_FILE);
	}

	// modify xml_tree[0] to get parameters from other xml_trees
	edit_tree(nfiles, nlmx, xml_tree, azi_1, azi_2, &nb_start, &nb_end);
	printf("Output burst from #%d to #%d ...\n", nb_start, nb_end);

	strcpy(tmp_str, argv[argc - 1]);
	strcat(tmp_str, ".xml");
	if ((XML_OUTPUT = fopen(tmp_str, "w")) == NULL)
		die("Couldn't open xml file: \n", tmp_str);
	// output the tree in the fromat of xml
	assemble_trees(nfiles, xml_tree, 0, 0, XML_OUTPUT);
	fclose(XML_OUTPUT);
	printf("XML file written...\n");

	tif = (TIFF **)malloc(nfiles * sizeof(TIFF *));
	for (ii = 0; ii < nfiles; ii++) {
		strcpy(tmp_str, argv[ii + 3]);
		strcat(tmp_str, ".tiff");
		if ((tif[ii] = TIFFOpen(tmp_str, "rb")) == NULL)
			die("Couldn't open tiff file: \n", tmp_str);
	}

	strcpy(tmp_str, argv[argc - 1]);
	strcat(tmp_str, ".tiff");
	if ((tif_out = TIFFOpen(tmp_str, "wb")) == NULL)
		die("Couldn't open tiff file: \n", tmp_str);
	search_tree(xml_tree[0], "/product/swathTiming/linesPerBurst/", tmp_str, 1, 0, 1);
	lpb = (int)str2double(tmp_str);

	assemble_tifs(tif, tif_out, nfiles, nb_start, nb_end, lpb);

	for (ii = 0; ii < nfiles; ii++)
		TIFFClose(tif[ii]);
	TIFFClose(tif_out);

	free(tif);
	for (ii = 1; ii < nfiles; ii++)
		free(xml_tree[ii]);
	free(xml_tree);
	return (1);
}

int assemble_tifs(TIFF **tif, TIFF *tif_out, int nfiles, int nb_start, int nb_end, int lpb) {

	int ii, jj;
	uint32 width, *height, height_all, ni = 0, nii, ni2 = 0;
	short *buf;
	uint16 s = 0;

	TIFFSetWarningHandler(NULL);
	TIFFGetField(tif[0], TIFFTAG_IMAGEWIDTH, &width);

	height_all = 0;
	height = (uint32 *)malloc(sizeof(uint32) * nfiles);
	for (ii = 0; ii < nfiles; ii++) {
		TIFFGetField(tif[ii], TIFFTAG_IMAGELENGTH, &height[ii]);
		height_all = height_all + height[ii];
	}

	TIFFSetField(tif_out, TIFFTAG_IMAGEWIDTH, width);
	TIFFSetField(tif_out, TIFFTAG_IMAGELENGTH, (nb_end - nb_start + 1) * lpb);
	TIFFSetField(tif_out, TIFFTAG_BITSPERSAMPLE, sizeof(short) * 8 * 2);
	TIFFSetField(tif_out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_COMPLEXINT);
	TIFFSetField(tif_out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);

	buf = (short *)_TIFFmalloc(TIFFScanlineSize(tif_out) * 2); // make some extra space in case the width differ for other images

	printf("Writing TIFF image Width(%d) X Height(%d)...\n", width, (nb_end - nb_start + 1) * lpb);

	for (ii = 0; ii < nfiles; ii++) {
		for (jj = 0; jj < height[ii]; jj++) {
			if (ni >= lpb * (nb_start - 1) && ni < lpb * nb_end) {
				nii = jj;
				TIFFReadScanline(tif[ii], buf, nii, s);
				if (TIFFFlushData(tif_out))
					TIFFWriteScanline(tif_out, buf, ni2, s);
				ni2++;
			}
			ni++;
		}
	}

	_TIFFfree(buf);
	free(height);

	return (1);
}

int add_index(struct tree *T, int ct, int ii) {

	int i;

	// printf("hahaha %d\n",ii);
	for (i = 0; i < ct; i++) {
		if (T[i].parent != -1)
			T[i].parent = T[i].parent + ii;
		if (T[i].sibl != -1)
			T[i].sibl = T[i].sibl + ii;
		if (T[i].sibr != -1)
			T[i].sibr = T[i].sibr + ii;
		if (T[i].firstchild != -1)
			T[i].firstchild = T[i].firstchild + ii;
	}

	// printf("hahaha\n");

	return (1);
}

double time2double(char *str) {
	double k;
	char s_name[200], s_out[200];

	cat_nums(s_name, str);
	str_date2JD(s_out, s_name);

	k = str2double(s_out);

	return (k);
}

int add_branch(struct tree **T, int nlmx, int qq, char *str, char *term, int mode) {

	double t1, t2;
	int ii1, jj1, kk1, ct, nn1, ii2, jj2, kk2, nn2;
	char tmp_c[60000];

	if (mode == 1) {
		// consecutive
		ii1 = search_tree(T[0], str, tmp_c, 3, 0, 1);
		nn1 = (int)str2double(tmp_c);
		kk1 = T[0][ii1].firstchild;
		for (jj1 = 1; jj1 < nn1; jj1++)
			kk1 = T[0][kk1].sibr;
		ii2 = search_tree(T[qq], str, tmp_c, 3, 0, 1);
		T[0][kk1].sibr = T[0][ii2 + nlmx * qq * 5].firstchild;
		nn2 = (int)str2double(tmp_c);
		kk2 = T[0][ii2 + nlmx * qq * 5].firstchild;
		for (jj2 = 1; jj2 < nn2; jj2++) {
			T[0][kk2].parent = ii1;
			kk2 = T[0][kk2].sibr;
		}
		T[0][kk2].parent = ii1;
		sprintf(tmp_c, "%s count=\"%d\"", term, nn1 + nn2);
		strcpy(T[0][ii1].name, tmp_c);
	}
	else if (mode == 2) {
		// in order of time, may have overlap
		ii1 = search_tree(T[0], str, tmp_c, 3, 0, 1);
		nn1 = (int)str2double(tmp_c);
		kk1 = T[0][ii1].firstchild;
		for (jj1 = 1; jj1 < nn1; jj1++)
			kk1 = T[0][kk1].sibr;
		t1 = time2double(T[0][T[0][T[0][kk1].firstchild].firstchild].name);
		ii2 = search_tree(T[qq], str, tmp_c, 3, 0, 1);
		nn2 = (int)str2double(tmp_c);
		// printf("%d %lf\n",nn1,t1);
		t2 = t1 - 1;
		ct = 0;
		kk2 = T[0][ii2 + nlmx * qq * 5].firstchild;
		t2 = time2double(T[0][T[0][T[0][kk2].firstchild].firstchild].name);
		while (t2 - 1e-6 / 86400.0 < t1) {
			kk2 = T[0][kk2].sibr;
            if (kk2 == -1) break;
			t2 = time2double(T[0][T[0][T[0][kk2].firstchild].firstchild].name);
			ct++;
		}

        if (kk2 != -1) { 
        // append the records from ct's record, otherwise don't do anything.
		    T[0][kk1].sibr = kk2;
		    kk2 = T[0][ii2 + nlmx * qq * 5].firstchild;
		    for (jj2 = 1; jj2 < nn2 - ct; jj2++) {
			    T[0][kk2].parent = ii1;
			    kk2 = T[0][kk2].sibr;
		    }
		    T[0][kk2].parent = ii1;
		    sprintf(tmp_c, "%s count=\"%d\"", term, nn1 + nn2 - ct);
		    strcpy(T[0][ii1].name, tmp_c);
        }
	}
	else if (mode == 3) {
		// same as mode 2, but time comes in second child
		ii1 = search_tree(T[0], str, tmp_c, 3, 0, 1);
		nn1 = (int)str2double(tmp_c);
		kk1 = T[0][ii1].firstchild;
		for (jj1 = 1; jj1 < nn1; jj1++)
			kk1 = T[0][kk1].sibr;
		t1 = time2double(T[0][T[0][T[0][T[0][kk1].firstchild].sibr].firstchild].name);
		ii2 = search_tree(T[qq], str, tmp_c, 3, 0, 1);
		nn2 = (int)str2double(tmp_c);
		// printf("%d %lf\n",nn1,t1);
		t2 = t1 - 1;
		ct = 0;
		kk2 = T[0][ii2 + nlmx * qq * 5].firstchild;
		t2 = time2double(T[0][T[0][T[0][T[0][kk2].firstchild].sibr].firstchild].name);
		while (t2 - 1e-6 / 86400.0 < t1) {
			kk2 = T[0][kk2].sibr;
            if (kk2 == -1) break;
			t2 = time2double(T[0][T[0][T[0][T[0][kk2].firstchild].sibr].firstchild].name);
			ct++;
		}
		// printf("%d %lf\n",nn2,t2);
        if (kk2 != -1) {
            // append the records from ct's record, otherwise don't do anything. 
		    T[0][kk1].sibr = kk2;
		    kk2 = T[0][ii2 + nlmx * qq * 5].firstchild;
		    for (jj2 = 1; jj2 < nn2 - ct; jj2++) {
			    T[0][kk2].parent = ii1;
			    kk2 = T[0][kk2].sibr;
		    }
		    T[0][kk2].parent = ii1;
		    sprintf(tmp_c, "%s count=\"%d\"", term, nn1 + nn2 - ct);
		    strcpy(T[0][ii1].name, tmp_c);
        }
	}

	return (1);
}

int edit_leaf(struct tree **T, int qq, char *str, int mode) {

	int ii, jj;
	double x1, x2;
	char tmp_c[200];

	if (mode == 1) {
		// add int numbers
		ii = search_tree(T[0], str, tmp_c, 1, 0, 1);
		x1 = str2double(tmp_c);
		jj = search_tree(T[qq], str, tmp_c, 1, 0, 1);
		x2 = str2double(tmp_c);
		sprintf(tmp_c, "%d", (int)(x1 + x2 + 1e-6));
		strcpy(T[0][T[0][ii].firstchild].name, tmp_c);
	}
	else if (mode == 2) {
		// replace with the second
		ii = search_tree(T[0], str, tmp_c, 1, 0, 1);
		jj = search_tree(T[qq], str, tmp_c, 1, 0, 1);
		strcpy(T[0][T[0][ii].firstchild].name, tmp_c);
	}

	return (1);
}

int edit_tree(int nfiles, int nlmx, struct tree **T, int azi_1, int azi_2, int *nb_start, int *nb_end) {

	/* note this is only editing firstchild and sibr, thus its not really getting
	 * a good tree structure */

	int ct, qq, i, ii, jj, lpb;
	char tmp_c[60000];
	double dta, t0, t1, t2;
	// char s_name[200],s_out[200];
	// double t1,t2;

	ct = 0;
	while (T[0][ct].sibr != -1 || T[0][ct].firstchild != -1) {
		if (T[0][ct].sibr != -1)
			ct = T[0][ct].sibr;
		else
			ct = T[0][ct].firstchild;
	}
	// printf("Original(first) xml has %d tree elements.\n",ct);

	if (nfiles > 1) {

		for (qq = 1; qq < nfiles; qq++) {
			ct = 0;
			while (T[qq][ct].sibr != -1 || T[qq][ct].firstchild != -1) {
				if (T[qq][ct].sibr != -1)
					ct = T[qq][ct].sibr;
				else
					ct = T[qq][ct].firstchild;
			}
			add_index(&T[0][nlmx * qq * 5], ct + 1, nlmx * qq * 5);
		}

		/* start editing the useful information to create one xml_tree */
		for (qq = 1; qq < nfiles; qq++) {
			/* start editing trunks of parameters */
			// burst List
			add_branch(T, nlmx, qq, "/product/swathTiming/burstList/", "burstList", 1);
			// orbit List
			add_branch(T, nlmx, qq, "/product/generalAnnotation/orbitList/", "orbitList", 2);
			// dc estimate
			add_branch(T, nlmx, qq, "/product/dopplerCentroid/dcEstimateList/", "dcEstimateList", 2);
			// azimuth fm rate
			add_branch(T, nlmx, qq, "/product/generalAnnotation/azimuthFmRateList/", "azimuthFmRateList", 2);
			// antenna pattern list
			add_branch(T, nlmx, qq, "/product/antennaPattern/antennaPatternList/", "antennaPatternList", 3);
			// geolocation grid, not used in our program, but anyway
			add_branch(T, nlmx, qq, "/product/geolocationGrid/geolocationGridPointList/", "geolocationGridPointList", 2);

			/* start editing other individual parameters */
			edit_leaf(T, qq, "/product/imageAnnotation/imageInformation/numberOfLines/", 1);
			edit_leaf(T, qq, "/product/adsHeader/stopTime/", 2);
		}
	}

	/* edit on how many burst to output */
	search_tree(T[0], "/product/swathTiming/linesPerBurst/", tmp_c, 1, 0, 1);
	lpb = (int)str2double(tmp_c);
	ii = search_tree(T[0], "/product/swathTiming/burstList/", tmp_c, 3, 0, 1);
	ct = (int)str2double(tmp_c);
	if (!(azi_1 == 0 && azi_2 == 0)) {
		printf("Editing burst output based on input azimuth references...\n");
		search_tree(T[0], "/product/imageAnnotation/imageInformation/azimuthTimeInterval/", tmp_c, 1, 0, 1);
		dta = str2double(tmp_c);
		search_tree(T[0], "/product/swathTiming/burstList/burst/azimuthTime/", tmp_c, 2, 4, 1);
		t0 = str2double(tmp_c);
		t2 = t0;
		t1 = t0;
		*nb_start = -1;
		*nb_end = -1;
		for (i = 2; i <= ct; i++) {
			search_tree(T[0], "/product/swathTiming/burstList/burst/azimuthTime/", tmp_c, 2, 4, i);
			t2 = str2double(tmp_c);
			if ((t2 - t0) * 86400.0 / dta >= azi_1 && (t1 - t0) * 86400.0 / dta <= azi_1)
				*nb_start = i - 1;
			if ((t2 - t0) * 86400.0 / dta >= azi_2 && (t1 - t0) * 86400.0 / dta <= azi_2)
				*nb_end = i - 1;
			t1 = t2;
		}

		if (*nb_start == -1)
			*nb_start = 1;
		if (*nb_end == -1)
			*nb_end = ct;

		ii = search_tree(T[0], "/product/imageAnnotation/imageInformation/numberOfLines/", tmp_c, 1, 0, 1);
		sprintf(tmp_c, "%d", (int)(-(*nb_start) + (*nb_end) + 1) * lpb);
		strcpy(T[0][T[0][ii].firstchild].name, tmp_c);
		ii = search_tree(T[0], "/product/imageAnnotation/imageInformation/productFirstLineUtcTime/", tmp_c, 1, 0, 1);
		search_tree(T[0], "/product/swathTiming/burstList/burst/azimuthTime/", tmp_c, 1, 4, *nb_start);
		strcpy(T[0][T[0][ii].firstchild].name, tmp_c);
		ii = search_tree(T[0], "/product/adsHeader/startTime/", tmp_c, 1, 0, 1);
		search_tree(T[0], "/product/swathTiming/burstList/burst/azimuthTime/", tmp_c, 1, 4, *nb_start);
		strcpy(T[0][T[0][ii].firstchild].name, tmp_c);

		if (*nb_end != ct) {
			ii = search_tree(T[0], "/product/imageAnnotation/imageInformation/productLastLineUtcTime/", tmp_c, 1, 0, 1);
			search_tree(T[0], "/product/swathTiming/burstList/burst/azimuthTime/", tmp_c, 1, 4, *nb_start);
			search_tree(T[0], "/product/swathTiming/burstList/burst/azimuthTime/", tmp_c, 1, 4, *nb_end + 1);
			strcpy(T[0][T[0][ii].firstchild].name, tmp_c);
			ii = search_tree(T[0], "/product/adsHeader/stopTime/", tmp_c, 1, 0, 1);
			search_tree(T[0], "/product/swathTiming/burstList/burst/azimuthTime/", tmp_c, 1, 4, *nb_end + 1);
			strcpy(T[0][T[0][ii].firstchild].name, tmp_c);
		}
		ii = search_tree(T[0], "/product/swathTiming/burstList/burst/", tmp_c, 1, 4, *nb_end);
		T[0][ii].sibr = -1;
		ii = search_tree(T[0], "/product/swathTiming/burstList/", tmp_c, 3, 0, 1);
		jj = search_tree(T[0], "/product/swathTiming/burstList/burst/", tmp_c, 1, 4, *nb_start);
		T[0][ii].firstchild = jj;

		ii = search_tree(T[0], "/product/swathTiming/burstList/", tmp_c, 1, 0, 1);
		sprintf(tmp_c, "%s count=\"%d\"", "burstList", *nb_end - *nb_start + 1);

		strcpy(T[0][ii].name, tmp_c);
	}
	else {
		*nb_start = 1;
		*nb_end = ct;
	}

	return (1);
}
