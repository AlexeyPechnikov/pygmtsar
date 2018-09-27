/***************************************************************************
 * Creator:  Xiaohua(Eric) XU                                              *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  07/09/2015                                                    *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * DATE                                                                    *
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

int pop_prm(struct PRM *, tree *, char *);
int pop_led(tree *, state_vector *);
int write_orb(state_vector *sv, FILE *fp, int);
int write_slc(FILE *, FILE *, int, int);

static int is_big_endian() {
	union {
		long l;
		char c[sizeof(long)];
	} u;
	u.l = 1;
	return (u.c[sizeof(long) - 1] == 1 ? 1 : -1);
}

#ifndef __CYGWIN__
static inline unsigned short bswap_16(unsigned short x) { return (x >> 8) | (x << 8); }

static inline unsigned int bswap_32(unsigned int x) { return (bswap_16(x & 0xffff) << 16) | (bswap_16(x >> 16)); }

/*static inline unsigned long long bswap_64(unsigned long long x) {
    return (((unsigned long long)bswap_32(x&0xffffffffull))<<32) |
(bswap_32(x>>32));
}*/
#endif

char *USAGE = "\n\nUsage: make_slc_tsx name_of_xml_file name_of_image_file name_output\n"
              "\nExample: make_slc_s1a "
              "TSX1_SAR__SSC______SM_S_SRA_20120615T162057_20120615T162105.xml "
              "IMAGE_HH_SRA_strip_007.cos TSX_HH_20120615\n"
              "\nOutput: TSX_HH_20120615.SLC TSX_HH_20120615.PRM TSX_HH_20120615.LED\n";

int main(int argc, char **argv) {

	FILE *XML_FILE, *OUTPUT_PRM, *OUTPUT_SLC, *OUTPUT_LED, *INPUT_SLC;
	char tmp_str[200];
	struct PRM prm;
	tree *xml_tree;
	state_vector sv[200];
	int rows, cols;
	int ch, n = 0, nc = 0, nlmx = 0;

	if (argc < 4)
		die(USAGE, "");

	// find the number of lines and the maximum line length of the xml file
	if ((XML_FILE = fopen(argv[1], "r")) == NULL)
		die("Couldn't open xml file: \n", argv[1]);
	while (EOF != (ch = fgetc(XML_FILE))) {
		++nc;
		if (ch == '\n') {
			++n;
			if (nc > nlmx)
				nlmx = nc;
			nc = 0;
		}
	}
	xml_tree = (struct tree *)malloc(n * 5 * sizeof(struct tree));
	fclose(XML_FILE);

	// generate the xml tree
	if ((XML_FILE = fopen(argv[1], "r")) == NULL)
		die("Couldn't open xml file: \n", argv[1]);
	get_tree(XML_FILE, xml_tree, 0);
	fclose(XML_FILE);

	// show_tree(xml_tree,0,0);

	// initiate the prm
	null_sio_struct(&prm);

	// generate the PRM file
	pop_prm(&prm, xml_tree, argv[3]);

	strcpy(tmp_str, argv[3]);
	strcat(tmp_str, ".PRM");
	if ((OUTPUT_PRM = fopen(tmp_str, "w")) == NULL)
		die("Couldn't open prm file: \n", tmp_str);
	put_sio_struct(prm, OUTPUT_PRM);
	fclose(OUTPUT_PRM);

	// generate the LED file
	n = pop_led(xml_tree, sv);

	strcpy(tmp_str, argv[3]);
	strcat(tmp_str, ".LED");
	if ((OUTPUT_LED = fopen(tmp_str, "w")) == NULL)
		die("Couldn't open led file: \n", tmp_str);
	write_orb(sv, OUTPUT_LED, n);
	fclose(OUTPUT_LED);

	// generate the SLC file
	if ((INPUT_SLC = fopen(argv[2], "rb")) == NULL)
		die("Couldn't open data file: \n", argv[2]);

	strcpy(tmp_str, argv[3]);
	strcat(tmp_str, ".SLC");
	if ((OUTPUT_SLC = fopen(tmp_str, "wb")) == NULL)
		die("Couldn't open tiff file: \n", tmp_str);
	search_tree(xml_tree, "/level1Product/productInfo/imageDataInfo/imageRaster/numberOfColumns/", tmp_str, 1, 0, 1);
	cols = (int)str2double(tmp_str);
	search_tree(xml_tree, "/level1Product/productInfo/imageDataInfo/imageRaster/numberOfRows/", tmp_str, 1, 0, 1);
	rows = (int)str2double(tmp_str);

	write_slc(INPUT_SLC, OUTPUT_SLC, rows, cols);
	fclose(INPUT_SLC);
	fclose(OUTPUT_SLC);

	// TIFFClose(TIFF_FILE);
	// fclose(OUTPUT_SLC);
}

int write_slc(FILE *input, FILE *slc, int rows, int cols) {

	int i, j, tj, k, tk, x;
	short *buf = malloc(sizeof(unsigned short) * (cols + 2) * 2);
	int bib, rsri, rs, as, bi, rtnb, tnl, asri, asfv, aslv, rsfv, rslv;

	i = is_big_endian();
	if (i == 1) {
		printf("System is Big Endian...\n");
	}
	else {
		printf("System is Little Endian...\n");
	}

	printf("Writing SLC..Image Size: %d X %d...\n", cols, rows);
	// fread(buf,sizeof(short),(cols+2)*2,input);
	j = 0;
	tj = rows;
	while (j < tj) {
		// first line
		fread(&bib, sizeof(int), 1, input);
		fread(&rsri, sizeof(int), 1, input);
		fread(&rs, sizeof(int), 1, input);
		fread(&as, sizeof(int), 1, input);
		fread(&bi, sizeof(int), 1, input);
		fread(&rtnb, sizeof(int), 1, input);
		fread(&tnl, sizeof(int), 1, input);
		fread(buf, sizeof(short), (cols - 5) * 2, input);
		// second line
		fread(buf, sizeof(short), 4, input);
		fread(&asri, sizeof(int), 1, input);
		fread(buf, sizeof(short), (cols - 1) * 2, input);
		// third line
		fread(buf, sizeof(short), 4, input);
		fread(&asfv, sizeof(int), 1, input);
		fread(buf, sizeof(short), (cols - 1) * 2, input);
		// fourth line
		fread(buf, sizeof(short), 4, input);
		fread(&aslv, sizeof(int), 1, input);
		fread(buf, sizeof(short), (cols - 1) * 2, input);

		if (i != 1) {
			bib = bswap_32(bib);
			rsri = bswap_32(rsri);
			rs = bswap_32(rs);
			as = bswap_32(as);
			bi = bswap_32(bi);
			rtnb = bswap_32(rtnb);
			tnl = bswap_32(tnl);
			asri = bswap_32(asri);
			asfv = bswap_32(asfv);
			aslv = bswap_32(aslv);
		}

		// printf("Burst Info: \n\tBytes in Burst: %u\n\tRange Samples:
		// %u\n\tAzimuth Samples: %u\n\tRange Total Number of Bytes: %u\n\tTotal
		// Number of Lines: %u\n",bib,rs,as,rtnb,tnl); printf("ASRI: %u    ASFV: %u
		// ASLV: %u\n",asri,asfv,aslv);

		// printf("Writing %u Bytes to SLC...\n",(tnl-4)*(rtnb-8));
		if (i != 1) {
			printf("Swaping Bytes...\n");
		}
		tk = tnl - 4;
		for (k = 0; k < tk; k++) {
			fread(&rsfv, sizeof(int), 1, input);
			fread(&rslv, sizeof(int), 1, input);
			fread(buf, sizeof(short), cols * 2, input);
			if (i != 1) {
				// printf("RSFV: %u    RSLV: %u\n",bswap_32(rsfv),bswap_32(rslv));
				for (x = 0; x < cols * 2; x++) {
					buf[x] = (short)bswap_16(buf[x]);
				}
			}
			fwrite(buf, sizeof(short), cols * 2, slc);
		}
		j = j + tk;
	}
	free(buf);
	return (1);
}

int write_orb(state_vector *sv, FILE *fp, int n) {
	int i;
	double dt;

	dt = trunc((sv[1].sec) * 1e4) / 1e4 - trunc((sv[0].sec) * 1e4) / 1e4;
	if (n <= 1)
		return (-1);
	fprintf(fp, "%d %d %d %.3lf %.3lf \n", n, sv[0].yr, sv[0].jd, sv[0].sec, dt);
	for (i = 0; i < n; i++) {
		fprintf(fp, "%d %d %.3lf %.6lf %.6lf %.6lf %.8lf %.8lf %.8lf \n", sv[i].yr, sv[i].jd, sv[i].sec, sv[i].x, sv[i].y,
		        sv[i].z, sv[i].vx, sv[i].vy, sv[i].vz);
	}
	return (1);
}

int pop_led(tree *xml_tree, state_vector *sv) {
	int i, count;
	char tmp_c[200];
	double tmp_d;

	search_tree(xml_tree, "/level1Product/platform/orbit/orbitHeader/numStateVectors/", tmp_c, 1, 0, 1);
	count = (int)str2double(tmp_c);
	for (i = 0; i < count; i++) {
		search_tree(xml_tree, "/level1Product/platform/orbit/stateVec/timeUTC/", tmp_c, 2, 4, i + 1);
		tmp_d = str2double(tmp_c);
		search_tree(xml_tree, "/level1Product/platform/orbit/stateVec/timeUTC/", tmp_c, 1, 4, i + 1);
		tmp_c[4] = '\0';
		sv[i].yr = (int)(str2double(tmp_c));
		sv[i].jd = (int)(tmp_d - trunc(tmp_d / 1000.0) * 1000.0);
		sv[i].sec = (tmp_d - trunc(tmp_d)) * 86400;
		search_tree(xml_tree, "/level1Product/platform/orbit/stateVec/posX/", tmp_c, 1, 4, i + 1);
		sv[i].x = str2double(tmp_c);
		search_tree(xml_tree, "/level1Product/platform/orbit/stateVec/posY/", tmp_c, 1, 4, i + 1);
		sv[i].y = str2double(tmp_c);
		search_tree(xml_tree, "/level1Product/platform/orbit/stateVec/posZ/", tmp_c, 1, 4, i + 1);
		sv[i].z = str2double(tmp_c);
		search_tree(xml_tree, "/level1Product/platform/orbit/stateVec/velX/", tmp_c, 1, 4, i + 1);
		sv[i].vx = str2double(tmp_c);
		search_tree(xml_tree, "/level1Product/platform/orbit/stateVec/velY/", tmp_c, 1, 4, i + 1);
		sv[i].vy = str2double(tmp_c);
		search_tree(xml_tree, "/level1Product/platform/orbit/stateVec/velZ/", tmp_c, 1, 4, i + 1);
		sv[i].vz = str2double(tmp_c);
	}
	printf("%d Lines Written for Orbit...\n", count);
	return (count);
}

int pop_prm(struct PRM *prm, tree *xml_tree, char *file_name) {
	char tmp_c[200];
	double tmp_d;
	int tmp_i;
	double c_speed = 299792458.0;

	// define some of the variables
	prm->first_line = 1;
	prm->st_rng_bin = 1;

	search_tree(xml_tree, "/level1Product/processing/processingParameter/rangeLooks/", tmp_c, 1, 0, 1);
	prm->nlooks = (int)str2double(tmp_c);
	prm->rshift = 0;
	prm->ashift = 0;
	prm->sub_int_r = 0.0;
	prm->sub_int_a = 0.0;
	prm->stretch_r = 0.0;
	prm->stretch_a = 0.0;
	prm->a_stretch_r = 0.0;
	prm->a_stretch_a = 0.0;
	prm->first_sample = 1;
	strasign(prm->dtype, "a", 0, 0);

	search_tree(xml_tree, "/level1Product/productInfo/imageDataInfo/imageRaster/rowSpacing/", tmp_c, 1, 0, 1);
	prm->fs = 1 / str2double(tmp_c); // rng_samp_rate
	prm->SC_identity = 7;            /* (1)-ERS1 (2)-ERS2 (3)-Radarsat (4)-Envisat (5)-ALOS
	                                    (6)-  (7)-TSX (8)-CSK (9)-RS2 (10) Sentinel-1a*/

	search_tree(xml_tree, "/level1Product/instrument/radarParameters/centerFrequency/", tmp_c, 1, 0, 1);
	prm->lambda = c_speed / str2double(tmp_c);

	search_tree(xml_tree,
	            "/level1Product/processing/processingParameter/rangeCompression/"
	            "chirps/referenceChirp/pulseLength/",
	            tmp_c, 1, 0, 1);
	tmp_d = str2double(tmp_c);
	search_tree(xml_tree,
	            "/level1Product/processing/processingParameter/rangeCompression/"
	            "chirps/referenceChirp/pulseBandwidth/",
	            tmp_c, 1, 0, 1);
	prm->chirp_slope = str2double(tmp_c) / tmp_d * pow(10.0, 9.0);
	search_tree(xml_tree,
	            "/level1Product/processing/processingParameter/rangeCompression/"
	            "chirps/referenceChirp/chirpSlope/",
	            tmp_c, 1, 0, 1);
	if (strcmp(tmp_c, "DOWN") == 0) {
		prm->chirp_slope = -1.0 * prm->chirp_slope;
	}

	prm->pulsedur = tmp_d / pow(10.0, 9.0);

	// search_tree(xml_tree,"/product/qualityInformation/qualityDataList/qualityData/imageQuality/imageStatistics/outputDataMean/re/",tmp_c,1,0,1);
	prm->xmi = 0.0; // str2double(tmp_c); //I_mean

	// search_tree(xml_tree,"/product/qualityInformation/qualityDataList/qualityData/imageQuality/imageStatistics/outputDataMean/im/",tmp_c,1,0,1);
	prm->xmq = 0.0; // str2double(tmp_c); //Q_mean

	search_tree(xml_tree, "/level1Product/productSpecific/complexImageInfo/commonPRF/", tmp_c, 1, 0, 1);
	prm->prf = str2double(tmp_c);

	search_tree(xml_tree, "/level1Product/productInfo/sceneInfo/rangeTime/firstPixel/", tmp_c, 1, 0, 1);
	prm->near_range = str2double(tmp_c) * c_speed / 2;
	prm->ra = 6378137.00; // equatorial_radius
	prm->rc = 6356752.31; // polar_radius

	search_tree(xml_tree, "/level1Product/productInfo/missionInfo/orbitDirection/", tmp_c, 1, 0, 1);
	strasign(prm->orbdir, tmp_c, 0, 0);

	search_tree(xml_tree, "/level1Product/productInfo/acquisitionInfo/lookDirection/", tmp_c, 1, 0, 1);
	strasign(prm->lookdir, tmp_c, 0, 0);

	strcpy(tmp_c, file_name);
	strcat(tmp_c, ".raw");
	strcpy(prm->input_file, tmp_c);

	strcpy(tmp_c, file_name);
	strcat(tmp_c, ".LED");
	strcpy(prm->led_file, tmp_c);

	strcpy(tmp_c, file_name);
	strcat(tmp_c, ".SLC");
	strcpy(prm->SLC_file, tmp_c);

	prm->SLC_scale = 1.0;

	search_tree(xml_tree, "/level1Product/productInfo/sceneInfo/start/timeUTC/", tmp_c, 2, 0, 1);
	prm->clock_start = str2double(tmp_c);
	search_tree(xml_tree, "/level1Product/productInfo/sceneInfo/start/timeUTC/", tmp_c, 1, 0, 1);
	tmp_c[4] = '\0';
	prm->SC_clock_start = prm->clock_start + 1000. * str2double(tmp_c);

	strasign(prm->iqflip, "n", 0, 0); // Flip_iq
	strasign(prm->deskew, "n", 0, 0); // deskew
	strasign(prm->offset_video, "n", 0, 0);

	search_tree(xml_tree, "/level1Product/productInfo/imageDataInfo/imageRaster/numberOfColumns/", tmp_c, 1, 0, 1);
	tmp_i = (int)str2double(tmp_c);
	prm->bytes_per_line = tmp_i * 4;
	prm->good_bytes = prm->bytes_per_line;
	prm->caltone = 0.0;
	prm->pctbwaz = 0.0;            // rm_az_band
	prm->pctbw = 0.2;              // rm_rng_band
	prm->rhww = 1.0;               // rng_spec_wgt
	strasign(prm->srm, "0", 0, 0); // scnd_rng_mig
	search_tree(xml_tree, "/level1Product/calibration/nominalGeometricPerformance/azimuthRes/", tmp_c, 1, 0, 1);
	prm->az_res = str2double(tmp_c);

	/*search_tree(xml_tree,"/level1Product/productInfo/acquisitionInfo/lookDirection/",tmp_c,1,0,1);
	prm->antenna_side = 1;
	if (strcmp(tmp_c,"R")==0){
	    prm->antenna_side = -1;
	}
	*/
	prm->fdd1 = 0.0;
	prm->fddd1 = 0.0;

	search_tree(xml_tree, "/level1Product/productInfo/imageDataInfo/imageRaster/numberOfRows/", tmp_c, 1, 0, 1);
	tmp_i = (int)str2double(tmp_c);
	prm->num_lines = tmp_i - tmp_i % 4;

	// search_tree(xml_tree,"/product/adsHeader/stopTime/",tmp_c,2,0,1);
	prm->SC_clock_stop = prm->SC_clock_start + prm->num_lines / prm->prf / 86400;
	prm->clock_stop = prm->clock_start + prm->num_lines / prm->prf / 86400;

	prm->nrows = prm->num_lines;
	prm->num_valid_az = prm->num_lines;
	prm->num_patches = 1;
	prm->num_rng_bins = prm->bytes_per_line / 4;
	prm->chirp_ext = 0;

	printf("PRM set for Image File...\n");
	return (1);
}
