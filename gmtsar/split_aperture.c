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
fcomplex fmean(short *, int);

char *USAGE = "\nUSAGE: split_aperture prm_file\n\n"
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
    fcomplex *fbuf, *fbuf_f, *fbuf_b,fm;
    size_t st_size;
    void *API = NULL;
    double f_f=0.0, f_b=0.0, w_f=0.0, w_b=0.0, prf;
    
    if ((API = GMT_Create_Session(argv[0], 0U, 0U, NULL)) == NULL)
        return EXIT_FAILURE;
    
    // read in PRM file
    get_prm(&p1,argv[1]);
    nl = p1.num_lines;
    num_rng_bins = p1.num_rng_bins;
    prf = p1.prf;
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
    fbuf = (fcomplex *)malloc(nffti * sizeof(fcomplex));
    fbuf_f = (fcomplex *)malloc(nffti * sizeof(fcomplex));
    fbuf_b = (fcomplex *)malloc(nffti * sizeof(fcomplex));
    
    // write original image
    fprintf(stderr,"Duplicating original SLCs (%d X %d) ...\n",nl,num_rng_bins);
    for (i=0;i<nl;i++) {
        //if (i%1000 == 0) fprintf(stderr,"%d ",i);
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
        die("can't open %s for writing", str_f);
    if ((fout_b = open(str_b, O_RDWR)) < 0)
        die("can't open %s for writing", str_b);
    
    fprintf(stderr,"Number of FFT %d ...\n",nffti);
    st_size = (size_t)4 * (size_t)nl * (size_t)num_rng_bins;
    
    if ((in = mmap(0, st_size, PROT_READ, MAP_SHARED, fin, 0)) == MAP_FAILED)
        die("mmap error for input", " ");
    if ((out_f = mmap(0, st_size, PROT_WRITE, MAP_SHARED, fout_f, 0)) == MAP_FAILED)
        die("mmap error for input", " ");
    if ((out_b = mmap(0, st_size, PROT_WRITE, MAP_SHARED, fout_b, 0)) == MAP_FAILED)
        die("mmap error for input", " ");

    fprintf(stderr,"Working on column ");

    for (j=0;j<num_rng_bins;j++) {
        if (j%1000 == 0) fprintf(stderr,"%d ",j);
        read_column_slc(in, buf, j, nl, num_rng_bins, nffti);
        fm = fmean(buf, num_rng_bins);
        //fm.r = 0.0; fm.i = 0.0;

        for (i=0;i<nffti;i++) {
            fbuf[i].r = (float)buf[2*i] - fm.r;
            fbuf[i].i = (float)buf[2*i+1] - fm.i;
        }
        GMT_FFT_1D(API, (float *)fbuf, nffti, GMT_FFT_FWD, GMT_FFT_COMPLEX);
        for (i=0;i<nffti;i++) {
            if (i < nffti/2) {
            //if (i < nffti) {
                fbuf_f[i].r = fbuf[i].r;
                fbuf_f[i].i = fbuf[i].i;
                fbuf_b[i].r = 0.0;
                fbuf_b[i].i = 0.0;
                //fbuf_b[i].r = fbuf[i].r;
                //fbuf_b[i].i = fbuf[i].i;
                w_f += (fbuf[i].r*fbuf[i].r+fbuf[i].i*fbuf[i].i);
                f_f += (fbuf[i].r*fbuf[i].r+fbuf[i].i*fbuf[i].i)*prf/nffti*i;
            }
            else {
                fbuf_b[i].r = fbuf[i].r;
                fbuf_b[i].i = fbuf[i].i;
                fbuf_f[i].r = 0.0;
                fbuf_f[i].i = 0.0;
                w_b += (fbuf[i].r*fbuf[i].r+fbuf[i].i*fbuf[i].i);
                f_b += (fbuf[i].r*fbuf[i].r+fbuf[i].i*fbuf[i].i)*prf/nffti*(i-nffti);
            }
        }
        GMT_FFT_1D(API, (float *)fbuf_f, nffti, GMT_FFT_INV, GMT_FFT_COMPLEX);
        GMT_FFT_1D(API, (float *)fbuf_b, nffti, GMT_FFT_INV, GMT_FFT_COMPLEX);
        for (i=0;i<nl;i++) {
            buf_f[2*i] = (int)round(fbuf_f[i].r);
            buf_f[2*i+1] = (int)round(fbuf_f[i].i);
            buf_b[2*i] = (int)round(fbuf_b[i].r);
            buf_b[2*i+1] = (int)round(fbuf_b[i].i);
        }
        write_column_slc(out_f, buf_f, j, nl, num_rng_bins);
        write_column_slc(out_b, buf_b, j, nl, num_rng_bins);
    }
    fprintf(stderr,"...\n");
    printf("average_spectrum_frequency_forward  = %.6f\n",f_f/w_f);
    printf("average_spectrum_frequency_backward = %.6f\n",f_b/w_b);
    printf("average_spectrum_frequency_separation = %.6f\n",f_f/w_f-f_b/w_b);
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
        //fprintf(stderr,"%d ",i);
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

fcomplex fmean(short *c, int N) {
    
    fcomplex m;
    int i;
    
    m.r = 0.0;
    m.i = 0.0;
    for (i=0;i<N;i++) {
        m.r += (float)c[2*i];
        m.i += (float)c[2*i+1];
    }
    m.r = m.r/(float)N;
    m.i = m.i/(float)N;
    
    return(m);
}




