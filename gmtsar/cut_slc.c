/*  $Id cut_slc.c 2018-08 Xiaohua XU$                                     */
/***************************************************************************
 * cut coregistered slc to a smaller area                                   *
 **************************************************************************/
/***************************************************************************
 * Creator:  Xiaohua Xu                                                    *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  07/11/2018                                                    *
 **************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 * DATE                                                                    *
 *                                                                         *
 ***************************************************************************/

#include "gmtsar.h"
#include <stdio.h>
#include <stdlib.h>

char *USAGE = "\nUsage: "
              "cut_slc stem.PRM new_stem range [prm_only]\n"
              "    stem.PRM      - PRM file for coregistered image to be cut\n"
              "    new_stem      - stem for newly generated SLC and PRM file\n"
              "    range         - range to cut the SLC e.g. 500/1500/900/3600\n"
              "\n"
              "Output: new_stem.PRM new_stem.SLC (old .LED file can still be used)\n\n";

int get_range(char *str, int *xl, int *xh, int *yl, int *yh) {

	int ii = 0, jj = 0, kk = 0, rr[4];
	char c;
	char tmp_c[128];

	c = str[ii];
	while (c != '\0') {
		if (c != '/') {
			tmp_c[jj] = c;
			jj++;
		}
		else if (c == '/') {
			tmp_c[jj] = '\0';
			rr[kk] = atoi(tmp_c);
			jj = 0;
			kk++;
		}
		ii++;
		c = str[ii];
	}
	tmp_c[jj] = c;
	rr[kk] = atoi(tmp_c);
	*xl = rr[0];
	*xh = rr[1];
	*yl = rr[2];
	*yh = rr[3];

	return (1);
}

int main(int argc, char **argv) {

	int ii, jj, xl, xh, yl, yh, nl, nr;
	short *buf_in, *buf_out;
	struct PRM pp;
	FILE *SLC_in, *SLC_out, *PRM_out;
	char str[1024];
	double dr;
    int do_slc = 1;

	if (argc < 4 || argc > 5)
		die(USAGE, "");
    if (argc == 5)
        do_slc = 0;

	get_prm(&pp, argv[1]);
	// xl = 0; xh = pp.num_rng_bins-1;
	// yl = 0; yh = pp.num_patches*num_valid_az-1;
	get_range(argv[3], &xl, &xh, &yl, &yh);
	nl = pp.num_patches * pp.num_valid_az;
	nr = pp.num_rng_bins;

	// printf("%d %d %d %d %d %d
	// %d\n",xl,xh,yl,yh,pp.num_rng_bins,pp.num_patches,pp.num_valid_az);

	//if (xl < 0 || xl > xh || xh > nr - 1 || yl < 0 || yl > yh || yh > nl - 1)
	if (xl < 0 || xl > xh || xh > nr || yl < 0 || yl > yh || yh > nl)
		die("wrong range ", "");

    if (do_slc == 1)
	    if ((SLC_in = fopen(pp.SLC_file, "rb")) == NULL)
		    die("Can't open SLCfile", pp.SLC_file);

	strcpy(str, argv[2]);
	// str[strlen(str)-3] = '\0';
	strcat(str, ".SLC");

    if (do_slc == 1)
	    if ((SLC_out = fopen(str, "wb")) == NULL)
		    die("Can't open SLCfile", str);

    if (do_slc == 1)
	    buf_in = (short *)malloc(pp.num_rng_bins * 2 * sizeof(short));

	pp.num_lines = (yh - yl) + 1;
	if (pp.num_lines % 4 != 0) {
		printf("Number of lines is set to a multiple of 4\n");
		pp.num_lines = pp.num_lines - pp.num_lines % 4;
	}


	pp.SC_clock_start = pp.SC_clock_start + ((double)xl*pp.stretch_a + (double)yl*(1+pp.a_stretch_a) + (double)(pp.nrows-pp.num_valid_az)/2.0) / pp.prf / 86400.0;
	pp.clock_start = pp.clock_start + ((double)xl*pp.stretch_a + (double)yl*(1+pp.a_stretch_a) + (double)(pp.nrows-pp.num_valid_az)/2.0) / pp.prf / 86400.0;
	pp.SC_clock_stop = pp.SC_clock_start + pp.num_lines / pp.prf / 86400.0;
	pp.clock_stop = pp.clock_start + pp.num_lines / pp.prf / 86400.0;

	pp.num_patches = 1;
	pp.num_valid_az = pp.num_lines;
	pp.nrows = pp.num_lines;

	dr = SOL / (2.0 * pp.fs);
	pp.near_range = pp.near_range + dr * xl * (1 + pp.stretch_r) + dr * yl * pp.a_stretch_r;

	pp.num_rng_bins = xh - xl + 1;
	if (pp.num_rng_bins % 4 != 0) {
		printf("Number of range pixels is set to a multiple of 4\n");
		pp.num_rng_bins = pp.num_rng_bins - pp.num_rng_bins % 4;
	}
	pp.bytes_per_line = pp.num_rng_bins * 4;
	pp.good_bytes = pp.bytes_per_line;

	strcpy(str, argv[2]);
	strcat(str, ".PRM");
	if ((PRM_out = fopen(str, "w")) == NULL)
		die("can't open prm file", str);
	put_sio_struct(pp, PRM_out);
	fclose(PRM_out);
	printf("New PRM file written ...\n");


    if (do_slc == 1) {
	    buf_out = (short *)malloc(pp.num_rng_bins * 2 * sizeof(short));

	    printf("Writing image (%d x %d) ...\n", pp.num_rng_bins, pp.num_lines);
	    // printf("Hahahaha %d %d %d %d\n",nr,nl,pp.num_rng_bins,pp.num_lines);
	    for (ii = 0; ii < nl; ii++) {
		    // printf("%d ",ii);
		    fread(buf_in, sizeof(short), nr * 2, SLC_in);
		    if (ii >= yl && ii <= yl + pp.num_lines - 1) {
			    for (jj = xl; jj <= xh; jj++) {
				    buf_out[jj * 2 - xl * 2] = buf_in[jj * 2];
				    buf_out[jj * 2 - xl * 2 + 1] = buf_in[jj * 2 + 1];
			    }
			    fwrite(buf_out, sizeof(short), pp.num_rng_bins * 2, SLC_out);
		    }
	    }
    
	fclose(SLC_out);
	fclose(SLC_in);
    free(buf_out);
    free(buf_in);
    }

	return (1);
}
