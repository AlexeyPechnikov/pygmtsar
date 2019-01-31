//
//  split_aperture.c
//  
//
//  Created by Xiaohua Xu on 1/30/19.
//
//  Program used for creating along-track multi-aperture interferogram
//
/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * DATE                                                                    *
 ***************************************************************************/

#include "gmt.h"
#include "gmtsar.h"
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/mman.h>
#include <sys/types.h>

int write_column_slc(short *, short *, int, int, int);
int read_column_slc(short *, short *, int, int, int, int);

char *USAGE = "nUSAGE: split_aperture prm_file\n\n"
              "   program used to split azimuth spectrum.\n\n"
              "   SLCs are bandpassed output is SLCF (forward) SLCB (backward)\n\n";

int main(int argc, char **argv) {
    if (argc != 2)
        die("",USAGE);
    
    FILE *SLC, *SLC_F, *SLC_B;
    struct PRM p1;
    int nl,num_rng_bins,i,j,fin, fout_f, fout_b, nffti;
    char str_f[1024],str_b[1024];
    short *buf, *buf_f, *buf_b;
    short *in, *out_f, *out_b;
    float *fbuf, *fbuf_f, *fbuf_b;
    size_t st_size;
    void *API = NULL;
    
    if ((API = GMT_Create_Session(argv[0], 0U, 0U, NULL)) == NULL)
        return EXIT_FAILURE;
    
    // read in PRM file
    get_prm(&p1,argv[1]);
    nl = p1.num_lines;
    num_rng_bins = p1.num_rng_bins;
    nffti = find_fft_length(nl);
    
    // open all SLCs
    if ((SLC = fopen(p1.SLC_file, "rb")) == NULL)
        die("Can't open ", p1.SLC_file);
    
    strcpy(str_f,p1.SLC_file);
    strcat(str_f,"F");
    if ((SLC_F = fopen(str_f, "wb")) == NULL)
        die("Can't open ", str_f);
    
    strcpy(str_b,p1.SLC_file);
    strcat(str_b,"B");
    if ((SLC_B = fopen(str_b, "wb")) == NULL)
        die("Can't open ", str_b);
    
    // malloc buf
    buf = (short *)malloc(nffti * 2 * sizeof(short));
    buf_f = (short *)malloc(nl * 2 * sizeof(short));
    buf_b = (short *)malloc(nl * 2 * sizeof(short));
    fbuf = (float *)malloc(nffti * 2 * sizeof(float));
    fbuf_f = (float *)malloc(nffti * 2 * sizeof(float));
    fbuf_b = (float *)malloc(nffti * 2 * sizeof(float));

    
    // write original image
    for (i=0;i<nl;i++) {
        fread(buf, sizeof(short),num_rng_bins*2, SLC);
        fwrite((void *)buf, sizeof(short),num_rng_bins*2,SLC_F);
        fwrite((void *)buf, sizeof(short),num_rng_bins*2,SLC_B);
    }
    
    fclose(SLC_F);
    fclose(SLC_B);
    fclose(SLC);
    
    // mmap the files
    if ((fin = open(p1.SLC_file, O_RDONLY)) < 0)
        die("can't open %s for reading", p1.SLC_file);
    if ((fout_f = open(str_f, O_RDWR)) < 0)
        die("can't open %s for reading", str_f);
    if ((fout_b = open(str_b, O_RDWR)) < 0)
        die("can't open %s for reading", str_b);
    

    st_size = (size_t)4 * (size_t)nl * (size_t)num_rng_bins;
    
    if ((in = mmap(0, st_size, PROT_READ, MAP_SHARED, fin, 0)) == MAP_FAILED)
        die("mmap error for input", " ");
    if ((out_f = mmap(0, st_size, PROT_READ, MAP_SHARED, fout_f, 0)) == MAP_FAILED)
        die("mmap error for input", " ");
    if ((out_b = mmap(0, st_size, PROT_READ, MAP_SHARED, fout_b, 0)) == MAP_FAILED)
        die("mmap error for input", " ");
    for (j=0;j<num_rng_bins;j++) {

        read_column_slc(in, buf, j, nl, num_rng_bins, nffti);

        for (i=0;i<nffti;i++) {
            fbuf[2*i] = buf[2*i];
            fbuf[2*i+1] = buf[2*i+1];
        }
        GMT_FFT_1D(API, (float *)fbuf, nffti, GMT_FFT_FWD, GMT_FFT_COMPLEX);
        for (i=0;i<nffti;i++) {
            if (i < nffti/2) {
                fbuf_f[2*i] = fbuf[2*i];
                fbuf_f[2*i+1] = fbuf[2*i+1];
                fbuf_b[2*i] = 0.0;
                fbuf_b[2*i+1] = 0.0;
            }
            else {
                fbuf_b[2*i] = fbuf[2*i];
                fbuf_b[2*i+1] = fbuf[2*i+1];
                fbuf_f[2*i] = 0.0;
                fbuf_f[2*i+1] = 0.0;
            }
        }
        GMT_FFT_1D(API, (float *)fbuf_f, nffti, GMT_FFT_INV, GMT_FFT_COMPLEX);
        GMT_FFT_1D(API, (float *)fbuf_b, nffti, GMT_FFT_INV, GMT_FFT_COMPLEX);
        for (i=0;i<nl;i++) {
            buf_f[2*i] = (int)round(fbuf_f[2*i]);
            buf_f[2*i+1] = (int)round(fbuf_f[2*i+1]);
            buf_b[2*i] = (int)round(fbuf_b[2*i]);
            buf_b[2*i+1] = (int)round(fbuf_b[2*i+1]);
        }

        write_column_slc(out_f, buf_f, j, nl, num_rng_bins);
        write_column_slc(out_b, buf_b, j, nl, num_rng_bins);

    }
    
    free(buf);free(buf_b);free(buf_f);
    free(fbuf);free(fbuf_b);free(fbuf_f);
    close(fin);close(fout_b);close(fout_f);
    
    return(1);
}

int read_column_slc(short *slc, short *buf, int col, int nl, int num_rng_bins, int nffti) {
    // read slc to buf
    short *pt;
    int i;
    
    for (i=0;i<nl;i++) {
        pt = slc; pt++;
        pt = slc + (size_t)(i*2*num_rng_bins) + (size_t)(col*2);
        buf[2*i] = pt[0];
        buf[2*i+1] = pt[1];
    }
    for (i=nl;i<nffti;i++) {
        buf[2*i] = 0;
        buf[2*i+1] = 0;
    }
    return(1);
}

int write_column_slc(short *slc, short *buf, int col, int nl, int num_rng_bins) {
    // write buf to slc
    short *pt;
    int i;
    
    for (i=0;i<nl;i++) {
        pt = slc; pt++;
        pt = slc + (size_t)(i*2*num_rng_bins) + (size_t)(col*2);
        pt[0] = buf[2*i];
        pt[1] = buf[2*i+1];
    }
    return(1);
}






