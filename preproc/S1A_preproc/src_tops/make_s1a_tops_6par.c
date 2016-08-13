/***************************************************************************
 * Creator:  David Sandwell and Xiaohua(Eric) XU                           *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  04/06/2015                                                    *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * Date   :  08/21/15 DTS begin modification for precision shifting        *
 * Date   :  10/24/15 DTS  Added the deramp and reramp                     *
 * Date   :  12/29/15 DTS added 3-parameter phase shift and reramp         *
 * Date   :  01/14/15 EXU added range stretch and line interpolator        *
 *                                                                         *
 ***************************************************************************/
#include"tiffio.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"PRM.h"
#include"lib_functions.h"
#include"stateV.h"
#include"xmlC.h"
#include"lib_defs.h"
#include "gmtsar.h"

typedef struct burst_bounds{
    int SL;
    int SC;
    int SH;
    int EL;
    int EC;
    int EH;
}burst_bounds;

int pop_led(struct tree *, struct state_vector *);
int write_orb(struct state_vector *sv, FILE *fp, int);
int pop_burst(struct PRM *, struct tree *, struct burst_bounds *, char *);
int shift_write_slcs(void *, struct PRM *, struct tree *,burst_bounds *, TIFF *, FILE *, FILE *, FILE *, int, double, double, double, double, double ,double);
int dramp_dmod(struct tree *, int, fcomplex *, int , int, double, double, double);
int azi_shift (void *, fcomplex *, int , int , fcomplex *, int ,  double, double);
int rng_shift (void *, fcomplex *, int , int , fcomplex *, int ,  double );
int rng_interp (fcomplex *, fcomplex *, int, int, int, int, double, double, double);
double sinc_kernel (double);

char *USAGE = "\nUsage: make_s1a_tops_6par xml_file tiff_file output mode rng_shift azi_shift stretch_a a_stretch_a stretch_r a_stretch_r\n"
              "         xml_file    - name of xml file \n" 
              "         tiff_file   - name of tiff file \n" 
              "         output      - stem name of output files .PRM, .LED, .SLC \n" 
              "         mode        - (0) no SLC; (1) center SLC; (2) high SLCH and low SLCL \n" 
              "         rng_shift   - fractional part of desired range shift \n" 
              "         azi_shift   - fractional part of desired azimuth shift \n" 
              "         stretch_a   - additional azimuth shift in range dir \n" 
              "         a_stretch_a - additional azimuth shift in azimuth dir \n" 
              "         stretch_r   - additional range shift in range dir \n"
              "         a_stretch_r - additional range shift in azimuth dir \n"
"\nExample: make_s1a_tops_6par s1a-s1-slc-vv-20140807.xml s1a-s1-slc-vv-20140807.tiff S1A20140807 1 0.4596 .9109 2.18065e-06 -3.63255e-06 1.37343e-05 -3.13352e-05 \n"
"\nOutput: S1A20140807.PRM S1A20140807.LED S1A20140807.SLC\n";

int main(int argc, char **argv){
    
    FILE *XML_FILE,*OUTPUT_PRM,*OUTPUT_LED;
    FILE *OUTPUT_SLCL=NULL,*OUTPUT_SLCC=NULL,*OUTPUT_SLCH=NULL;
    TIFF *TIFF_FILE;
    char tmp_str[200];
    struct PRM prm;
    struct tree *xml_tree;
    struct state_vector sv[400];
    struct burst_bounds bb[200];
    int ch, n=0, nc=0, nlmx=0, imode=0;
    double rng, azi, stretch_a, a_stretch_a, stretch_r, a_stretch_r;

    // Begin: Initializing new GMT5 session 
    void    *API = NULL; // GMT API control structure 
    if ((API = GMT_Create_Session (argv[0], 0U, 0U, NULL)) == NULL) return EXIT_FAILURE;
    
    if (argc < 11) die (USAGE,"");

    // find the number of lines and the maximum line length of the xml file
    if ((XML_FILE = fopen(argv[1],"r")) == NULL) die("Couldn't open xml file: \n",argv[1]);
    while (EOF != (ch=fgetc(XML_FILE))) {
        ++nc;
        if (ch == '\n') {
           ++n;
           if(nc > nlmx) nlmx = nc;
           nc=0;
        }
    }
    fclose(XML_FILE);
    xml_tree = (struct tree *)malloc(n*5*sizeof(struct tree));
    
    // generate the xml tree
    if ((XML_FILE = fopen(argv[1],"r")) == NULL) die("Couldn't open xml file: \n",argv[1]);
    get_tree(XML_FILE,xml_tree,1);
    fclose(XML_FILE);

    //show_tree(xml_tree,0,0);
    
    // generate the LED file
    n = pop_led(xml_tree,sv);
    strcpy(tmp_str,argv[3]);
    strcat(tmp_str,".LED");
    if ((OUTPUT_LED = fopen(tmp_str,"w")) == NULL) die ("Couldn't open led file: \n",tmp_str);
    write_orb(sv,OUTPUT_LED,n);
    fclose(OUTPUT_LED);
    
    // initiate the prm
    null_sio_struct(&prm);

    // analyze the burst and generate the PRM
    pop_burst(&prm,xml_tree,bb,argv[3]);
    
    // open the TIFF file and the three SLC files SLCL-low  SLCC-center and SLCH-high
    TIFFSetWarningHandler(NULL);
    if ((TIFF_FILE = TIFFOpen(argv[2], "rb")) == NULL) die ("Couldn't open tiff file: \n",argv[2]);

    // open output files depending on the imode
    imode = atoi(argv[4]);
    if(imode == 1) {
       strcpy(tmp_str,argv[3]);
       strcat(tmp_str,".SLC");
       if ((OUTPUT_SLCC = fopen(tmp_str,"wb")) == NULL) die ("Couldn't open slc(C) file: \n",tmp_str);
    }
    else if(imode == 2){
       strcpy(tmp_str,argv[3]);
       strcat(tmp_str,".SLCL");
       if ((OUTPUT_SLCL = fopen(tmp_str,"wb")) == NULL) die ("Couldn't open slc(L) file: \n",tmp_str);
       strcpy(tmp_str,argv[3]);
       strcat(tmp_str,".SLCH");
       if ((OUTPUT_SLCH = fopen(tmp_str,"wb")) == NULL) die ("Couldn't open slc(H) file: \n",tmp_str);
     }
     else {
     // don't open any files
     }

    /* apply range and azimuth shifts to each burst and write the three SLC files SLCL SLC and SLCH depending on imode */
    /* the negative signs are needed for all but the azimuth stretch versus range */
    rng = atof(argv[5]);
    azi = -atof(argv[6]);
    stretch_r = atof(argv[9]);
    a_stretch_r = atof(argv[10]);
    stretch_a = atof(argv[7]);
    a_stretch_a = -atof(argv[8]);
    prm.sub_int_r = rng;
    prm.sub_int_a = -azi;
    prm.stretch_a = stretch_a;
    prm.a_stretch_a = -a_stretch_a;
    prm.stretch_r = stretch_r;
    prm.a_stretch_r = a_stretch_r;

    shift_write_slcs(API,&prm,xml_tree,bb,TIFF_FILE,OUTPUT_SLCL,OUTPUT_SLCC,OUTPUT_SLCH,imode,rng,azi,stretch_a,a_stretch_a,stretch_r,a_stretch_r);
    
    TIFFClose(TIFF_FILE);
    if(imode == 2) fclose(OUTPUT_SLCL);
    if(imode == 1) fclose(OUTPUT_SLCC);
    if(imode == 2) fclose(OUTPUT_SLCH);
    
    strcpy(tmp_str,argv[3]);
    strcat(tmp_str,".PRM");
    if ((OUTPUT_PRM = fopen(tmp_str,"w")) == NULL) die ("Couldn't open prm file: \n",tmp_str);
    put_sio_struct(prm, OUTPUT_PRM);
    fclose(OUTPUT_PRM);
    free(xml_tree);
    if (GMT_Destroy_Session (API)) return EXIT_FAILURE;     /* Remove the GMT machinery */
    return (EXIT_SUCCESS);
}

int write_orb(state_vector *sv, FILE *fp, int n){
    int i;
    double dt;
    
    dt = trunc((sv[1].sec)*100.0)/100.0-trunc((sv[0].sec)*100.0)/100.0;
    if(n<=1) return(-1);
    fprintf(fp,"%d %d %d %.6lf %lf \n",n,sv[0].yr,sv[0].jd,sv[0].sec,dt);
    for(i=0;i<n;i++){
        fprintf(fp,"%d %d %.6lf %.6lf %.6lf %.6lf %.8lf %.8lf %.8lf \n",sv[i].yr,sv[i].jd,sv[i].sec,sv[i].x,sv[i].y,sv[i].z,sv[i].vx,sv[i].vy,sv[i].vz);
    }
    return(1);
}

int pop_led(tree *xml_tree,state_vector *sv){
    int i,count;
    char tmp_c[200];
    double tmp_d;
    
    search_tree(xml_tree,"/product/generalAnnotation/orbitList/",tmp_c,3,0,1);
    count = (int)str2double(tmp_c);
    for (i=0;i<count;i++){
        search_tree(xml_tree,"/product/generalAnnotation/orbitList/orbit/time/",tmp_c,2,4,i+1);
        tmp_d = str2double(tmp_c);
        search_tree(xml_tree,"/product/generalAnnotation/orbitList/orbit/time/",tmp_c,1,4,i+1);
        tmp_c[4] = '\0';
        sv[i].yr = (int)(str2double(tmp_c));
        sv[i].jd = (int)(tmp_d - trunc(tmp_d/1000.0)*1000.0);
        sv[i].sec = (tmp_d - trunc(tmp_d))*86400;
        search_tree(xml_tree,"/product/generalAnnotation/orbitList/orbit/position/x/",tmp_c,1,4,i+1);
        sv[i].x = str2double(tmp_c);
        search_tree(xml_tree,"/product/generalAnnotation/orbitList/orbit/position/y/",tmp_c,1,4,i+1);
        sv[i].y = str2double(tmp_c);
        search_tree(xml_tree,"/product/generalAnnotation/orbitList/orbit/position/z/",tmp_c,1,4,i+1);
        sv[i].z = str2double(tmp_c);
        search_tree(xml_tree,"/product/generalAnnotation/orbitList/orbit/velocity/x/",tmp_c,1,4,i+1);
        sv[i].vx = str2double(tmp_c);
        search_tree(xml_tree,"/product/generalAnnotation/orbitList/orbit/velocity/y/",tmp_c,1,4,i+1);
        sv[i].vy = str2double(tmp_c);
        search_tree(xml_tree,"/product/generalAnnotation/orbitList/orbit/velocity/z/",tmp_c,1,4,i+1);
        sv[i].vz = str2double(tmp_c);
    }
    printf("Writing %d lines for orbit...\n",count);
    return(count);
}

int pop_burst(struct PRM *prm, tree *xml_tree, struct burst_bounds *bb, char *file_name){
    
    char tmp_c[200],tmp_cc[60000];
    double tmp_d,dt,t[100];
    int i,j,k,nl=0,nlf,ntl,count,lpb,tmp_i,flag,flag0;
    int k_start, kC;
    int *kF, *ksa, *ksr, *kea, *ker, *kover;
    double t0=-1.,time;
    char *cflag,*cflag_orig;
    
    // define some of the variables
    prm->first_line = 1;
    prm->st_rng_bin = 1;
    search_tree(xml_tree,"/product/imageAnnotation/processingInformation/swathProcParamsList/swathProcParams/rangeProcessing/numberOfLooks/",tmp_c,1,0,1);
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
    strasign(prm->dtype,"a",0,0);
    search_tree(xml_tree,"/product/generalAnnotation/productInformation/rangeSamplingRate/",tmp_c,1,0,1);
    prm->fs = str2double(tmp_c);//rng_samp_rate
    prm->SC_identity = 10; /* (1)-ERS1 (2)-ERS2 (3)-Radarsat (4)-Envisat (5)-ALOS (6)-  (7)-TSX (8)-CSK (9)-RS2 (10) Sentinel-1a*/
    search_tree(xml_tree,"/product/generalAnnotation/productInformation/radarFrequency/",tmp_c,1,0,1);
    prm->lambda = SOL/str2double(tmp_c);
    search_tree(xml_tree,"/product/generalAnnotation/downlinkInformationList/downlinkInformation/downlinkValues/txPulseLength/",tmp_c,1,0,1);
    tmp_d = str2double(tmp_c);
    search_tree(xml_tree,"/product/imageAnnotation/processingInformation/swathProcParamsList/swathProcParams/rangeProcessing/lookBandwidth/",tmp_c,1,0,1);
    prm->chirp_slope = str2double(tmp_c)/tmp_d;
    prm->pulsedur = tmp_d;
    search_tree(xml_tree,"/product/qualityInformation/qualityDataList/qualityData/imageQuality/imageStatistics/outputDataMean/re/",tmp_c,1,0,1);
    //prm->xmi = str2double(tmp_c); //I_mean
    prm->xmi = 0.0; //I_mean
    search_tree(xml_tree,"/product/qualityInformation/qualityDataList/qualityData/imageQuality/imageStatistics/outputDataMean/im/",tmp_c,1,0,1);
    //prm->xmq = str2double(tmp_c); //Q_mean
    prm->xmq = 0.0; //Q_mean
    search_tree(xml_tree,"/product/imageAnnotation/imageInformation/azimuthTimeInterval/",tmp_c,1,0,1);
    dt = str2double(tmp_c);
    prm->prf = 1/dt;
    search_tree(xml_tree,"/product/imageAnnotation/imageInformation/slantRangeTime/",tmp_c,1,0,1);
    prm->near_range = str2double(tmp_c)*SOL/2;
    prm->ra = 6378137.00; //equatorial_radius
    prm->rc = 6356752.31; //polar_radius
    search_tree(xml_tree,"/product/generalAnnotation/productInformation/pass/",tmp_c,1,0,1);
    strasign(prm->orbdir,tmp_c,0,0);
    strasign(prm->lookdir,"R",0,0);
    strcpy(tmp_c,file_name);
    strcat(tmp_c,".raw");
    strcpy(prm->input_file,tmp_c);
    strcpy(tmp_c,file_name);
    strcat(tmp_c,".LED");
    strcpy(prm->led_file,tmp_c);
    strcpy(tmp_c,file_name);
    strcat(tmp_c,".SLC");
    strcpy(prm->SLC_file,tmp_c);
    prm->SLC_scale = 1.0;
    strasign(prm->iqflip,"n",0,0); //Flip_iq
    strasign(prm->deskew,"n",0,0); //deskew
    strasign(prm->offset_video,"n",0,0);
    search_tree(xml_tree,"/product/imageAnnotation/imageInformation/numberOfSamples/",tmp_c,1,0,1);
    tmp_i = (int)str2double(tmp_c) - (int)str2double(tmp_c)%4;
    prm->bytes_per_line = tmp_i*4;
    prm->good_bytes = prm->bytes_per_line;
    prm->num_rng_bins = prm->bytes_per_line/4;
    prm->caltone = 0.0;
    prm->pctbwaz = 0.0; //rm_az_band
    prm->pctbw = 0.2; //rm_rng_band
    prm->rhww = 1.0; //rng_spec_wgt
    strasign(prm->srm,"0",0,0); //scnd_rng_mig
    prm->az_res = 0.0;
    //prm->antenna_side = -1;
    prm->fdd1 = 0.0;
    prm->fddd1 = 0.0;
    
    // manipulate the lines
    search_tree(xml_tree,"/product/swathTiming/burstList/",tmp_c,3,0,1);
    count = (int)str2double(tmp_c);
    //count = 1;
    search_tree(xml_tree,"/product/imageAnnotation/imageInformation/numberOfLines/",tmp_c,1,0,1);
    ntl = (int)str2double(tmp_c);
    search_tree(xml_tree,"/product/swathTiming/linesPerBurst/",tmp_c,1,4,0);
    lpb = (int)str2double(tmp_c);
    nlf = count*lpb;
    // allocate some memory
    kF = (int *)malloc(nlf*sizeof(int));
    ksa = (int *)malloc((count+1)*sizeof(int));
    ksr = (int *)malloc((count+1)*sizeof(int));
    kea = (int *)malloc((count+1)*sizeof(int));
    ker = (int *)malloc((count+1)*sizeof(int));
    kover = (int *)malloc((count+1)*sizeof(int));
    cflag_orig = (char *)malloc(sizeof(char)*80*( lpb + 1 ) );
    cflag = cflag_orig;
    
    search_tree(xml_tree,"/product/imageAnnotation/imageInformation/productFirstLineUtcTime/",tmp_c,2,0,1);
    prm->clock_start = str2double(tmp_c);
    search_tree(xml_tree,"/product/imageAnnotation/imageInformation/productFirstLineUtcTime/",tmp_c,1,0,1);
    tmp_c[4] = '\0';
    prm->SC_clock_start = prm->clock_start + 1000.*str2double(tmp_c);
    
    k = 0;
    flag0 = -1;
    for (i=1;i<=count;i++){
        search_tree(xml_tree,"/product/swathTiming/burstList/burst/azimuthAnxTime/",tmp_c,1,4,i);
        t[i] = str2double(tmp_c);
        search_tree(xml_tree,"/product/swathTiming/burstList/burst/firstValidSample/",tmp_cc,1,4,i);
	strcpy(cflag,tmp_cc);
        for (j=0;j<lpb;j++){
            flag = (int)strtol(cflag,&cflag,10);
            time = t[i] + (double)j/prm->prf;
            kF[k] = -1;
        // don't use flagged data
            if(flag >= 0) {
		if(t0 < 0.) {
                t0=time;
                k_start = k;
                }
       //  the kF array uniquely maps the line on the input file into the SLC line on output
              kF[k] = (int)((time - t0)*prm->prf + .1);
              nl = kF[k];
            }
       // calculate the start and end index for each burst
            if(flag > flag0) {
              ksa[i] = kF[k];
              ksr[i] = j;
            }
            if(flag < flag0) {
              kea[i] = kF[k-1];
              ker[i] = j-1;
            }
            flag0 = flag;
            k++;
        }
    }
    // calculate the burst overlap.  The first one is zero.
    kover[0] = 0;
    for (i=1;i<=count;i++){
        kover[i] = 0;
        if(i < count) kover[i] = kea[i] - ksa[i+1] +1;
    }
    // calculate the approximate center shift
    kC = kover[1]/2;
    // calculate the relative start S and end E lines for each burst in the L and H configurations
    for (i=1;i<=count;i++){
        bb[i].SL = ksr[i];
        bb[i].SH = ksr[i] + kover[i-1];
        bb[i].EL = ker[i] - kover[i];
        bb[i].EH = ker[i];
    }
    // make the center start end locations
    for (i=1;i<count;i++){
        bb[i].EC = bb[i].EL + kC;
        bb[i+1].SC = bb[i+1].SL + kC;
    }
    bb[1].SC = bb[1].SL;
    bb[count].EC = bb[i].EL;

    // make sure the number of lines is divisible by 4
    prm->num_lines = nl-nl%4;
    // advance the start time to account for the zero lines at the start of the first burst
    prm->SC_clock_start = prm->SC_clock_start + (double)(k_start)/prm->prf/86400.;
    prm->SC_clock_stop = prm->SC_clock_start + prm->num_lines/prm->prf/86400;
    prm->clock_start = prm->clock_start + (double)(k_start)/prm->prf/86400.;
    prm->clock_stop = prm->clock_start + prm->num_lines/prm->prf/86400;
    prm->nrows = prm->num_lines;
    prm->num_valid_az = prm->num_lines;
    prm->num_patches = 1;
    prm->chirp_ext = 0;
    free(kF);
    free(ksa);
    free(ksr);
    free(kea);
    free(ker);
    free(kover);
    free(cflag_orig); /* freeing this memory causes a crash!!!! */
    return(1);
}
 
int shift_write_slcs (void *API, struct PRM *prm, tree *xml_tree, burst_bounds *bb, TIFF *tif, FILE *slcl, FILE *slcc, FILE *slch, int imode, double rng, double azi, double stretch_a, double a_stretch_a, double stretch_r, double a_stretch_r){

    uint16 s=0;
    uint16 *buf;
    uint32 it;
    short *tmp, *brst;
    float *rtmp;
    int ii,jj,nl,k,k2,kk;
    int count, lpb, nlf, width2, nclip=0;
    uint32 width,height,widthi;
    char tmp_c[200];
    fcomplex *cbrst, *cramp;
    fcomplex *fft_vec_rng, *fft_vec_azi;
    int ranfft_rng,ranfft_azi;
    double rtest, itest, azi_new;
    int cl;

    // get the burst information and compare with the TIFF file
    search_tree(xml_tree,"/product/swathTiming/burstList/",tmp_c,3,0,1);
    count = (int)str2double(tmp_c);
    search_tree(xml_tree,"/product/swathTiming/linesPerBurst/",tmp_c,1,4,0);
    lpb = (int)str2double(tmp_c);
    nlf = count*lpb;

    // get the width the file and make width divisible by 4
    TIFFGetField(tif,TIFFTAG_IMAGEWIDTH,&widthi);
    TIFFGetField(tif,TIFFTAG_IMAGELENGTH,&height);
    width = widthi - widthi%4;
    width2 = 2*width;

    if(nlf != height) {
       fprintf(stderr,"XML and TIFF file disagree on image height %d %d \n", nlf, height);
    }

    // allocate memory for the burst and other arrays
    brst = (short *)malloc(lpb*width2*sizeof(short));
    cbrst = (fcomplex *) malloc(lpb*width*sizeof(fcomplex));
    cramp = (fcomplex *) malloc(lpb*width*sizeof(fcomplex));
    buf = (uint16 *)_TIFFmalloc(TIFFScanlineSize(tif));
    tmp = (short *)malloc(width*2*sizeof(short));
    rtmp = (float *)malloc(width*2*sizeof(float));
    ranfft_rng = fft_bins(width);
    ranfft_azi = fft_bins(lpb);
    fft_vec_rng = (fcomplex *) malloc(ranfft_rng*sizeof(fcomplex));
    fft_vec_azi = (fcomplex *) malloc(ranfft_azi*sizeof(fcomplex));
    nl = prm-> num_lines;
 
    // don't do anything if imode = 0
    if(imode == 0) return(1);
 
    printf("Writing SLC..Image Size: %d X %d...\n",width,nl);
    it = 0;

    // loop over the bursts
    cl = 0;

    for (kk=1;kk<=count;kk++){
    fprintf(stderr,"burst # %d \n",kk);

    // read a burst into memory
        for (ii=0;ii<lpb;ii++){
            TIFFReadScanline(tif, buf, it, s);
            for (jj=0;jj<width2;jj++){
                k = ii*width2+jj;
                brst[k]=(short)buf[jj];
	    }
            it++;
        }

    // load the burst into a complex float array
	for (ii=0;ii<lpb;ii++){
	    for (jj=0;jj<width;jj++){
                k = ii*width+jj;
                k2 = ii*width2+jj*2;
                cbrst[k].r= (float) brst[k2];
                cbrst[k].i= (float) brst[k2+1];
            }
        }

   // don't need this code if there are no shifts
    if(azi != 0. || rng != 0.){

    // shift the burst in range if != 0.
       if(rng != 0.){
           rng_interp(cbrst,fft_vec_rng,cl,bb[kk].SC,bb[kk].EC,width,rng,stretch_r,a_stretch_r);
           //rng_shift (API,cbrst,lpb,width,fft_vec_rng,ranfft_rng,rng);
       }

    // complute the complex deramp_demod array for each burst with no shift
       dramp_dmod (xml_tree,kk,cramp,lpb,width,0.,0.,0.);

    // apply the dramp_dmod  
	for (ii=0;ii<lpb;ii++){
	    for (jj=0;jj<width;jj++){
                k = ii*width+jj;
                cbrst[k] = Cmul(cbrst[k],cramp[k]);
            }
        }

    // shift the burst cramp in azimuth if azi != 0.
       if(azi != 0.){
           azi_new=azi+(cl+(bb[kk].EC-bb[kk].SC+1)/2.)*a_stretch_a;
           azi_shift (API,cbrst,lpb,width,fft_vec_azi,ranfft_azi,azi_new,stretch_a);
       }

    // now make a shifted version of the dramp_dmod
       dramp_dmod (xml_tree,kk,cramp,lpb,width,0.0,azi_new,stretch_a);

    // reramp the slc
	for (ii=0;ii<lpb;ii++){
	    for (jj=0;jj<width;jj++){
                k = ii*width+jj;
                cramp[k].i = -cramp[k].i;
                cbrst[k] = Cmul(cbrst[k],cramp[k]);
            }
        }

     }
  
    // unload the float complex array into a short burst array, multiply by 2 and clip if needed
 	for (ii=0;ii<lpb;ii++){
	    for (jj=0;jj<width;jj++){
                k = ii*width+jj;
                k2 = ii*width2+jj*2;
		rtest = 2.*cbrst[k].r;
		itest = 2.*cbrst[k].i;
		brst[k2]   = (short) clipi2(rtest);
                brst[k2+1] = (short) clipi2(itest);
                if((int)(rtest) > I2MAX || (int)(itest) > I2MAX ){
                  nclip++;
                }
            }
        }

    // write the burst in L, C, and H configurations
          for (ii=0;ii<lpb;ii++){
            // write low
            if(ii >= bb[kk].SL && ii <= bb[kk].EL){
            for (jj=0;jj<width2;jj++){
                k = ii*width2+jj;
                tmp[jj] = brst[k];
            }
            if(imode == 2) fwrite(tmp,sizeof(short),width*2,slcl);
            }
    
            // write center
            if(ii >= bb[kk].SC && ii <= bb[kk].EC){
            for (jj=0;jj<width2;jj++){
                k = ii*width2+jj;
                tmp[jj] = brst[k];
            }
            if(imode == 1) fwrite(tmp,sizeof(short),width*2,slcc);
            cl++;
            }

            // write high
            if(ii >= bb[kk].SH && ii <= bb[kk].EH){
            for (jj=0;jj<width2;jj++){
                k = ii*width2+jj;
                tmp[jj] = brst[k];
            }
            if(imode == 2) fwrite(tmp,sizeof(short),width*2,slch);
            }
        }
    }

    fprintf(stderr,"number of points clipped to short int %d \n", nclip);
    _TIFFfree(buf);
    free(fft_vec_rng);
    free(fft_vec_azi);
    free(brst);
    free(cbrst);
    free(cramp);
    free(tmp);
    free(rtmp);
    return(1);
}

int dramp_dmod (tree *xml_tree, int nb, fcomplex *cramp, int lpb, int width, double rng, double azi, double stretch_a){

/*  this is a routine to apply an azimuth shift to a burst of TOPS data */
    int ii,jj,k;
    char tmp_c[200];
    double kpsi,fc,dta,dts=0.,ts0,tau0,vx,vy,vz;
    double ks,vtot,ka,taus,phase,pramp,pmod;
    double fnc[3],fka[3];
    double *eta,*etaref,*kt,*fnct;
    double azi_rng;
    double t_brst, t1, t2;

    // get all the parameters needed for the remod_deramp
    search_tree(xml_tree,"/product/generalAnnotation/productInformation/radarFrequency/",tmp_c,1,0,1);
    fc = str2double(tmp_c);
    search_tree(xml_tree,"/product/generalAnnotation/productInformation/azimuthSteeringRate/",tmp_c,1,0,1);
    kpsi=PI*str2double(tmp_c)/180.;
    search_tree(xml_tree,"/product/imageAnnotation/imageInformation/azimuthTimeInterval/",tmp_c,1,0,1);
    dta = str2double(tmp_c);
    search_tree(xml_tree,"/product/generalAnnotation/productInformation/rangeSamplingRate/",tmp_c,1,0,1);
    if(str2double(tmp_c) != 0.) dts = 1./str2double(tmp_c);
    search_tree(xml_tree,"/product/imageAnnotation/imageInformation/slantRangeTime/",tmp_c,1,0,1);
    ts0 = str2double(tmp_c);
    search_tree(xml_tree,"/product/generalAnnotation/azimuthFmRateList/azimuthFmRate/t0/",tmp_c,1,0,1);
    tau0 = str2double(tmp_c);
    
    search_tree(xml_tree,"/product/swathTiming/burstList/burst/azimuthTime/",tmp_c,2,4,nb);
    t_brst = str2double(tmp_c)+dta*(double)lpb/2.0/86400.0;

    // first find the parameters that depends on burst, find the neareast one.
    search_tree(xml_tree,"/product/dopplerCentroid/dcEstimateList/",tmp_c,3,0,1); 
    jj = (int)str2double(tmp_c);
    t2 = 0.0;
    ii = 0;
    while(t2<t_brst && ii < jj){
        t1 = t2;
	ii++;
	search_tree(xml_tree,"/product/dopplerCentroid/dcEstimateList/dcEstimate/azimuthTime/",tmp_c,2,4,ii);
	t2 = str2double(tmp_c);
    }
    if (t_brst-t1 < t2-t_brst) jj = ii-1;
    else jj = ii;
    //fprintf(stderr,"finding the %d DcPolynomial\n",jj);

    search_tree(xml_tree,"/product/dopplerCentroid/dcEstimateList/dcEstimate/dataDcPolynomial/",tmp_c,1,4,jj);
    str2dbs(fnc,tmp_c);


    search_tree(xml_tree,"/product/generalAnnotation/azimuthFmRateList/",tmp_c,3,0,1);
    jj = (int)str2double(tmp_c);
    t2 = 0.0;
    ii = 0;
    while(t2<t_brst && ii < jj){
	t1 = t2;
	ii++;
        search_tree(xml_tree,"/product/generalAnnotation/azimuthFmRateList/azimuthFmRate/azimuthTime/",tmp_c,2,4,ii);
	t2 = str2double(tmp_c);
    }
    if (t_brst-t1 < t2-t_brst) jj = ii-1;
    else jj = ii;
    //fprintf(stderr,"finding the %d AziFmRate\n",jj);

    ii = search_tree(xml_tree,"/product/generalAnnotation/azimuthFmRateList/azimuthFmRate/t0/",tmp_c,1,0,1);
    if (xml_tree[xml_tree[ii].sibr].sibr < 0){
        search_tree(xml_tree, "/product/generalAnnotation/azimuthFmRateList/azimuthFmRate/azimuthFmRatePolynomial/", tmp_c,1,4,jj);
        str2dbs(fka,tmp_c);
    }
    else{
        search_tree(xml_tree,"/product/generalAnnotation/azimuthFmRateList/azimuthFmRate/c0/",tmp_c,1,4,jj);
        fka[0] = str2double(tmp_c);
        search_tree(xml_tree,"/product/generalAnnotation/azimuthFmRateList/azimuthFmRate/c1/",tmp_c,1,4,jj);
        fka[1] = str2double(tmp_c);
        search_tree(xml_tree,"/product/generalAnnotation/azimuthFmRateList/azimuthFmRate/c2/",tmp_c,1,4,jj);
        fka[2] = str2double(tmp_c);
    }

    // find the velocity by linearly interpolating between 2 neareast point
    search_tree(xml_tree,"/product/generalAnnotation/orbitList/",tmp_c,3,0,1);
    jj = (int)str2double(tmp_c);

    t2 = 0.0;
    ii = 0;
    while(t2 < t_brst && ii < jj){
        t1 = t2;
        ii++;
        search_tree(xml_tree,"/product/generalAnnotation/orbitList/orbit/time/",tmp_c,2,4,ii);
        t2 = str2double(tmp_c);
    }

    search_tree(xml_tree,"/product/generalAnnotation/orbitList/orbit/velocity/x/",tmp_c,1,4,ii-1);
    vx = str2double(tmp_c);
    search_tree(xml_tree,"/product/generalAnnotation/orbitList/orbit/velocity/x/",tmp_c,1,4,ii);
    vx = (str2double(tmp_c)*(t2-t_brst)+vx*(t_brst-t1))/(t2-t1);
    search_tree(xml_tree,"/product/generalAnnotation/orbitList/orbit/velocity/y/",tmp_c,1,4,ii-1);
    vy = str2double(tmp_c);
    search_tree(xml_tree,"/product/generalAnnotation/orbitList/orbit/velocity/y/",tmp_c,1,4,ii);
    vy = (str2double(tmp_c)*(t2-t_brst)+vy*(t_brst-t1))/(t2-t1);
    search_tree(xml_tree,"/product/generalAnnotation/orbitList/orbit/velocity/z/",tmp_c,1,4,ii-1);
    vz = str2double(tmp_c);
    search_tree(xml_tree,"/product/generalAnnotation/orbitList/orbit/velocity/z/",tmp_c,1,4,ii);
    vz = (str2double(tmp_c)*(t2-t_brst)+vz*(t_brst-t1))/(t2-t1);

    // malloc the memory
    eta    = (double *) malloc(lpb*sizeof(double));
    etaref = (double *) malloc(width*sizeof(double));
    kt     = (double *) malloc(width*sizeof(double));
    fnct   = (double *) malloc(width*sizeof(double));
    // malloc the memory
    eta    = (double *) malloc(lpb*sizeof(double));
    etaref = (double *) malloc(width*sizeof(double));
    kt     = (double *) malloc(width*sizeof(double));
    fnct   = (double *) malloc(width*sizeof(double));

    // compute velocity and then Doppler rate
    vtot = sqrt(vx*vx+vy*vy+vz*vz);
    ks = 2.*vtot*fc*kpsi/SOL;

    for (ii=0;ii<lpb;ii++){
        //eta[ii] = ((double)(ii-lpb/2)-azi)*dta;
        eta[ii] = ((double)ii-(double)lpb/2.+.5-azi)*dta;
    }
    
    for (jj=0;jj<width;jj++){
    // This is silly because ts0 = tau0 in the test I did
        taus=ts0+((double)jj-rng)*dts-tau0;
        ka=fka[0]+fka[1]*taus+fka[2]*taus*taus;
        kt[jj]=ka*ks/(ka-ks);
        fnct[jj]=fnc[0]+fnc[1]*taus+fnc[2]*taus*taus;
        etaref[jj]=-fnct[jj]/ka+fnc[0]/fka[0];
    }

    for (ii=0;ii<lpb;ii++){
        for (jj=0;jj<width;jj++){
            k = ii*width+jj;
            /* include the phase shift with range */
            azi_rng = jj*stretch_a;
            pramp = -M_PI*kt[jj]*pow((eta[ii] + azi_rng*dta - etaref[jj]),2);
            pmod  = -2.*M_PI*fnct[jj]*eta[ii];
            phase = pramp + pmod;
            cramp[k]=Cexp(phase);
        }
    }

    return(1);
}

int azi_shift (void *API, fcomplex *cbrst, int lpb, int width, fcomplex *fft_vec, int ranfft,  double azi, double stretch_a){

/*  this is a routine to apply an azimuth shift to a burst of TOPS data */
    int ii,jj,k;
    double sumr, sumi, azi_new;

       // loop over columns of the burst
       for (jj=0;jj<width;jj++){
           azi_new = azi + jj*stretch_a;
           // compute the mean value
           for (ii=0;ii<lpb;ii++){
               k = ii*width+jj;
               sumr=sumr+cbrst[k].r;
               sumi=sumi+cbrst[k].i;
           }
           sumr=sumr/(double)lpb;
           sumi=sumi/(double)lpb;
           // load the complex vector and append the mean
           for (ii=0;ii<ranfft;ii++){
               if(ii < lpb){
               k = ii*width+jj;
                   fft_vec[ii].r = cbrst[k].r;
                   fft_vec[ii].i = cbrst[k].i;
               }
               else{
                   fft_vec[ii].r = sumr;
                   fft_vec[ii].i = sumi;
               }
           }
           // apply the shift theorem
           shift(API,ranfft,fft_vec,azi_new);
           // unload the complex vector 
           for (ii=0;ii<lpb;ii++){
               k = ii*width+jj;
                   cbrst[k].r = fft_vec[ii].r;
                   cbrst[k].i = fft_vec[ii].i;
               }
       }
    return(1);
}

int rng_shift (void *API, fcomplex *cbrst, int lpb, int width, fcomplex *fft_vec, int ranfft,  double rng){

/*  this is a routine to apply a range shift to a burst of TOPS data */
    int ii,jj,k;

       // loop over rows of the burst
       for (ii=0;ii<lpb;ii++){
           // load the complex vector
           for (jj=0;jj<ranfft;jj++){
               if(jj < width){
               k = ii*width+jj;
                   fft_vec[jj].r = cbrst[k].r;
                   fft_vec[jj].i = cbrst[k].i;
               }
               else{
                   fft_vec[jj].r = 0.0f;
                   fft_vec[jj].i = 0.0f;
               }
           }
        // apply the shift theorem
        shift(API,ranfft,fft_vec,rng);
        // unload the complex vector
        for (jj=0;jj<width;jj++){
            k = ii*width+jj;
            cbrst[k].r = fft_vec[jj].r;
            cbrst[k].i = fft_vec[jj].i;
        }
       }
    return(1);
}

int rng_interp (fcomplex *cbrst, fcomplex *ctmp, int cl, int lstart, int lend, int width, double rng, double stretch_r, double a_stretch_r){

/* this is a routine to shift and stretch along range for TOPS data */
    
    int ii,jj,ns2,j,k,kk;
    double ra,w,wsum;  
 
    ns2=NS/2-1;

    for (ii=lstart;ii<=lend;ii++){
        // read in the line
        for (jj=0;jj<width;jj++){
            kk = ii*width+jj;
            ctmp[jj].i = cbrst[kk].i;
            ctmp[jj].r = cbrst[kk].r;
        }
        // compute the output
        for (jj=0;jj<width;jj++){
            ra = (double)jj + ( rng + (double)jj*stretch_r + (cl+ii-lstart)*a_stretch_r );
            j = (int)ra;
            // make sure the point is safe for interpolation, use original input if can not be interpolated
            if ( j-ns2 < 0 || j+ns2+1 >= width) continue;
            
            kk = ii*width+jj;
            cbrst[kk].i = 0.0;
            cbrst[kk].r = 0.0;
            wsum = 0.0;
            for (k=j-ns2;k<=j+ns2+1;k++) {
                w = sinc_kernel(fabs(ra-(double)k));
                cbrst[kk].i += w * ctmp[k].i;
                cbrst[kk].r += w * ctmp[k].r;
                wsum += w;
            }
            cbrst[kk].i = cbrst[kk].i/wsum;
            cbrst[kk].r = cbrst[kk].r/wsum;
        }
    }
    return(1);
}


double sinc_kernel (double x) {

    double arg, f;
        
        arg = fabs(PI*x);
        if(arg > 0.){
            f=sin(arg)/arg;
        }
        else {
            f=1.;
        }
        return (f);
}

