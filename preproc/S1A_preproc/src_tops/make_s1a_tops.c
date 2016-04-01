/***************************************************************************
 * Creator:  Xiaohua(Eric) XU and David Sandwell                           *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  02/01/2016                                                    *
 ***************************************************************************/

/***************************************************************************
 * Adapted from previous make_slc_s1a_tops code, start with point-by-point *
 * precise shift for aligning two images                                   *
 ***************************************************************************/

/***************************************************************************
 * Modification history of previous make_slc_s1a_tops code:                *
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
double dramp_dmod(struct tree *, int, fcomplex *, int , int, int, struct GMT_GRID *, struct GMT_GRID *, int,int);
double shift_write_slc(void *, struct PRM *, struct tree *, burst_bounds *, int, TIFF *, FILE *, FILE *, FILE *, char *, char *, char *);
int shift_burst(fcomplex *, int , int , int , struct GMT_GRID *, struct GMT_GRID *,int);
void fbisinc (double *, fcomplex *, int, int, fcomplex *);


char *USAGE =  "\nUsage: make_slc_s1a_tops xml_file tiff_file output mode dr.grd da.grd [aux_file]\n"
              "         xml_file    - name of xml file \n"
              "         tiff_file   - name of tiff file \n"
              "         output      - stem name of output files .PRM, .LED, .SLC \n"
              "         mode        - (0) no SLC; (1) center SLC; (2) high SLCH and low SLCL \n"
              "         dr.grd      - range shift table to be read in \n"
              "         da.grd      - azimuth shift table to be read in \n"
              "         aux_file    - used to compute elevation antenna pattern\n"
"\nExample: make_slc_s1a_tops s1a-s1-slc-vv-20140807.xml s1a-s1-slc-vv-20140807.tiff S1A20140807 1 dr.grd da.grd [aux_file]\n"
"\n         make_slc_s1a_tops s1a-s1-slc-vv-20140807.xml s1a-s1-slc-vv-20140807.tiff S1A20140807 1 [aux_file]\n"
"\nOutput: mode 1: S1A20140807.PRM S1A20140807.LED S1A20140807.SLC\n"
"\n        mode 2: S1A20140807.PRM S1A20140807.LED S1A20140807.SLCH S1A20140807.SLCL S1A20140807.BB"
"\nNote: if dr and da are not given, SLCs will be written with no shifts.\n"
"\n      aux_file is needed for acquisitions acquired before Mar 2015.\n";


int main(int argc, char **argv){

    FILE *XML_FILE,*OUTPUT_PRM,*OUTPUT_LED;
    FILE *OUTPUT_SLCL=NULL,*OUTPUT_SLCC=NULL,*OUTPUT_SLCH=NULL,*BB=NULL;
    TIFF *TIFF_FILE;
    char tmp_str[200],rshifts[200],ashifts[200],aux_file[200];
    struct PRM prm;
    struct tree *xml_tree;
    struct state_vector sv[400];
    struct burst_bounds bb[200];
    int ch, n=0, nc=0, nlmx=0, imode=0;
    double spec_sep = 0.0, dta = 0.0;

    // Begin: Initializing new GMT5 session
    void    *API = NULL; // GMT API control structure
    if ((API = GMT_Create_Session (argv[0], 0U, 0U, NULL)) == NULL) return EXIT_FAILURE;

    if (argc == 5) {
        rshifts[0] = '\0';
        ashifts[0] = '\0';
        aux_file[0] = '\0';
    }
    else if (argc == 6) {
        rshifts[0] = '\0';
        ashifts[0] = '\0';
        strcpy(aux_file,argv[7]);
    }
    else if (argc == 7) {
        strcpy(rshifts,argv[5]);
        strcpy(ashifts,argv[6]);
        aux_file[0] = '\0';
    }
    else if (argc == 8) {
        strcpy(rshifts,argv[5]);
        strcpy(ashifts,argv[6]);
        strcpy(aux_file,argv[7]);
    }
    else{
        die (USAGE,"");
    }
    
    imode = atoi(argv[4]);
    
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
       strcpy(tmp_str,argv[3]);
       strcat(tmp_str,".BB");
       if ((BB = fopen(tmp_str,"w")) == NULL) die ("Couldn't open burst boundary file: \n",tmp_str);
    }
    else {
     // don't open any files
    }

    /* apply range and azimuth shifts to each burst and write the three SLC files SLCL SLC and SLCH depending on imode */
    if (imode == 1 || imode == 2) {
        spec_sep = shift_write_slc(API,&prm,xml_tree,bb,imode,TIFF_FILE,OUTPUT_SLCL,OUTPUT_SLCC,OUTPUT_SLCH,rshifts,ashifts,aux_file);
    }
    /* shift applied */
    
    // write the boundary file if mode == 2
    if (imode == 2) {
        search_tree(xml_tree,"/product/swathTiming/burstList/",tmp_str,3,0,1);
        nc = (int)str2double(tmp_str);
        search_tree(xml_tree,"/product/imageAnnotation/imageInformation/azimuthTimeInterval/",tmp_str,1,0,1);
        dta = str2double(tmp_str);
        search_tree(xml_tree,"/product/imageAnnotation/imageInformation/numberOfSamples/",tmp_str,1,0,1);
        fprintf(BB,"%d %d %.6f %.12f\n",nc,(int)str2double(tmp_str) - (int)str2double(tmp_str)%4,spec_sep,dta);
        for (n=1;n<=nc;n++) {
            fprintf(BB,"%d %d %d %d %d %d\n",bb[n].SL,bb[n].SC,bb[n].SH,bb[n].EL,bb[n].EC,bb[n].EH);
        }
    }
 
    TIFFClose(TIFF_FILE);
    if(imode == 2) fclose(OUTPUT_SLCL);
    if(imode == 1) fclose(OUTPUT_SLCC);
    if(imode == 2) fclose(OUTPUT_SLCH);
    if(imode == 2) fclose(BB);

    strcpy(tmp_str,argv[3]);
    strcat(tmp_str,".PRM");
    if ((OUTPUT_PRM = fopen(tmp_str,"w")) == NULL) die ("Couldn't open prm file: \n",tmp_str);
    put_sio_struct(prm, OUTPUT_PRM);
    fclose(OUTPUT_PRM);
    free(xml_tree);
    if (GMT_Destroy_Session (API)) return EXIT_FAILURE;     /* Remove the GMT machinery */
    return (EXIT_SUCCESS);

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
    free(cflag_orig); 
    return(1);
}


double dramp_dmod (struct tree *xml_tree, int nb, fcomplex *cramp, int lpb, int width, int al_start, struct GMT_GRID *R, struct GMT_GRID *A, int bshift, int imode){

/*  this is a routine to apply an azimuth shift to a burst of TOPS data */
    int ii,jj,k;
    char tmp_c[200];
    double kpsi,fc,dta,dts=0.,ts0,tau0,vx,vy,vz;
    double ks,vtot,ka,taus,phase,pramp,pmod;
    double fnc[3],fka[3];
    double *eta,*etaref,*kt,*fnct;
    double t_brst, t1, t2;
    double dr,da;
    double sum_spec_sep=0.0;
    
    // get all the parameters needed for the remod_deramp
    // find the constant parameters
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

    // compute velocity and then Doppler rate
    vtot = sqrt(vx*vx+vy*vy+vz*vz);
    ks = 2.*vtot*fc*kpsi/SOL;

    if (imode == 0) {
        for (ii=0;ii<lpb;ii++){
            eta[ii] = ((double)ii-(double)lpb/2.+.5)*dta;
        }

        for (jj=0;jj<width;jj++){
        // This is silly because ts0 = tau0 in the test I did
            taus=ts0+((double)jj)*dts-tau0;
            ka=fka[0]+fka[1]*taus+fka[2]*taus*taus;
            kt[jj]=ka*ks/(ka-ks);
            fnct[jj]=fnc[0]+fnc[1]*taus+fnc[2]*taus*taus;
            etaref[jj]=-fnct[jj]/ka+fnc[0]/fka[0];
        }

        for (ii=0;ii<lpb;ii++){
            for (jj=0;jj<width;jj++){
                k = ii*width+jj;
                pramp = -M_PI*kt[jj]*pow((eta[ii] - etaref[jj]),2);
                pmod  = -2.*M_PI*fnct[jj]*eta[ii];
                phase = pramp + pmod;
                cramp[k]=Cexp(phase);
            }
        }
    }
    else if (imode == 1){
        for (ii=0;ii<lpb;ii++){
            for (jj=0;jj<width;jj++){
                k = ii*width+jj;
                if (floor((al_start+ii)/R->header->inc[GMT_Y]+0.5) < 0 || (int)floor((al_start+ii)/R->header->inc[GMT_Y]+0.5) >= R->header->ny){
                    cramp[k].r=1;
                    cramp[k].i=0;
                }                  
                else{
                    dr = R->data[(int)(jj/R->header->inc[GMT_X]+0.5)+R->header->nx*(int)((al_start+ii)/R->header->inc[GMT_Y]+0.5)];
                    da = A->data[(int)(jj/A->header->inc[GMT_X]+0.5)+A->header->nx*(int)((al_start+ii)/A->header->inc[GMT_Y]+0.5)] - (double)bshift;
                    
                    eta[0] = ((double)ii - (double)lpb/2.+.5 + da)*dta;
                    taus=ts0+((double)jj + dr)*dts-tau0;
                    ka=fka[0]+fka[1]*taus+fka[2]*taus*taus;
                    kt[0]=ka*ks/(ka-ks);
                    fnct[0]=fnc[0]+fnc[1]*taus+fnc[2]*taus*taus;
                    etaref[0]=-fnct[0]/ka+fnc[0]/fka[0];
                    pramp = -M_PI*kt[0]*pow((eta[0] - etaref[0]),2);
                    pmod  = -2.*M_PI*fnct[0]*eta[0];
                    phase = pramp + pmod;
                    cramp[k]=Cexp(phase);
                }
            }
        }
    }
    else if (imode == 2){
        for (jj=0;jj<width;jj++){
            taus=ts0+((double)jj)*dts-tau0;
            ka=fka[0]+fka[1]*taus+fka[2]*taus*taus;
            kt[jj]=ka*ks/(ka-ks);
        }

        for (ii=0;ii<lpb;ii++) {
            for (jj=0;jj<width;jj++){
                sum_spec_sep += kt[jj];
            }
        }
    }

    return(sum_spec_sep);
}


double shift_write_slc(void *API,struct PRM *prm,struct tree *xml_tree,struct burst_bounds *bb,int imode,TIFF *tif,FILE *slcl,FILE *slcc,FILE *slch,char *dr_table,char *da_table, char *aux_file) {

    uint16 s=0;
    uint16 *buf;
    uint32 it;
    short *tmp, *brst;
    float *rtmp;
    int ii,jj,nl,k,k2,kk;
    int count, lpb, nlf, width2, nclip=0;
    uint32 width, height, widthi;
    char tmp_c[200];
    fcomplex *cbrst, *cramp;
    float rtest,itest;
    int al_start,cl;
    int bshift=0;
    double sum_spec_sep = 0.0, spec_sep = 0.0,dta;

    struct GMT_GRID *R = NULL, *A = NULL;

    if (dr_table[0] != '\0' && da_table[0] != '\0') {
        fprintf(stderr,"Reading in range and azimuth shifts table...\n");
        if ((R = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, dr_table, NULL)) == NULL) die("cannot open range shift tables",dr_table);
        if ((A = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, da_table, NULL)) == NULL) die("cannot open azimuth shift tables",dr_table);    
        if (R->header->inc[GMT_X] != A->header->inc[GMT_X]) die("shift table size does not match","");
        if (R->header->inc[GMT_Y] != A->header->inc[GMT_Y]) die("shift table size does not match","");
        if (GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, dr_table, R) == NULL) return EXIT_FAILURE;
        if (GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, da_table, A) == NULL) return EXIT_FAILURE;
    }

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
    nl = prm-> num_lines;
    
    if (imode == 1) printf("Writing SLC..Image Size: %d X %d...\n",width,nl);
    else if (imode == 2) printf("Writing SLCL & SLCH..\n");
    it = 0;

    // fix the shift when there are burst with no overlap
    cl = 0;

    if (dr_table[0] != '\0' && da_table[0] != '\0') {
        for (kk=1;kk<=count;kk++){
            cl = cl + (bb[kk].EL-bb[kk].SL+1);
            //fprintf(stderr,"SC:%d, EC:%d, SH:%d, EH:%d, SL:%d, EL:%d\n",bb[kk].SC,bb[kk].EC,bb[kk].SH,bb[kk].EH,bb[kk].SL,bb[kk].EL);
            if (floor(A->data[0]/(double)cl+0.2) == 1) {
                bshift = cl;
                fprintf(stderr,"Image has %d bursts mismatch...making corrections to azimuth shifts (%d)...\n",kk,bshift);
            }
            else if(ceil(A->data[0]/(double)cl-0.2) == -1) {
                bshift = -cl;
                fprintf(stderr,"Image has %d bursts mismatch...making corrections to azimuth shifts (%d)...\n",kk,bshift);
            }
        }
    
        if (bshift != 0) {
            prm->ashift = bshift;
            fprintf(stderr,"Updated ashift is %d\n",bshift);
        }
        cl = 0;
    }

    // loop over the bursts
    for (kk=1;kk<=count;kk++){
        fprintf(stderr,"burst # %d \n",kk);
        if (imode == 2) {
            sum_spec_sep = dramp_dmod(xml_tree,kk,cramp,lpb,width,al_start,R,A,bshift,2);
            spec_sep += sum_spec_sep;
        }
    
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

        // do not shift anything if dr and da is not given
        if (dr_table[0] != '\0' && da_table[0] != '\0') {
            // complute the complex deramp_demod array for each burst with no shift
            al_start = cl-bb[kk].SC;
            
	    dramp_dmod(xml_tree,kk,cramp,lpb,width,al_start,R,A,bshift,0);

            // apply the dramp_dmod
            for (ii=0;ii<lpb;ii++){
                for (jj=0;jj<width;jj++){
                    k = ii*width+jj;
                    cbrst[k] = Cmul(cbrst[k],cramp[k]);
                }
            }  
       
            // shift the burst with the given table
            shift_burst(cbrst,al_start,lpb,width,R,A,bshift); 

            dramp_dmod(xml_tree,kk,cramp,lpb,width,al_start,R,A,bshift,1);   

            //prm->rshift = -11;
            //prm->ashift = -15;
            // reramp the slc
            for (ii=0;ii<lpb;ii++){
                for (jj=0;jj<width;jj++){
                    k = ii*width+jj;
                    cramp[k].i = -cramp[k].i;
                    cbrst[k] = Cmul(cbrst[k],cramp[k]);
                }
            }
        }

        // compute the elevation antenna pattern (EAP) change if aux_file is specified
/*
        if (aux_file[0] != '\0') {
            compute_eap(cramp,lpb,width,aux_file);
            for (ii=0;ii<lpb;ii++){
                for (jj=0;jj<width;jj++){
                    k = ii*width+jj;
                    cbrst[k] = Cmul(cbrst[k],cramp[k]);
                }
            }
        }
*/
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
            if(kk > 1 && ii >= bb[kk].SL && ii <= bb[kk].SH-1){
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
            if(kk < count && ii >= bb[kk].EL+1 && ii <= bb[kk].EH){
                for (jj=0;jj<width2;jj++){
                    k = ii*width2+jj; 
                    tmp[jj] = brst[k];
                }
                if(imode == 2) fwrite(tmp,sizeof(short),width*2,slch);
            }
        }
    }

    search_tree(xml_tree,"/product/imageAnnotation/imageInformation/azimuthTimeInterval/",tmp_c,1,0,1);
    dta = str2double(tmp_c);
    if(imode == 2) {
        spec_sep = spec_sep*dta/count/width;
    }

    fprintf(stderr,"number of points clipped to short int %d \n", nclip);
    _TIFFfree(buf);
    free(brst);
    free(cbrst);
    free(cramp);
    free(tmp);
    free(rtmp);
    return(spec_sep);
    
}

/*
int compute_eap(fcomplex *cbrst, tree *xml_tree char *aux_file, int mode) {
    // 
    //  mode 1=S1_HH, 2=S1_HV, 3=S1_VV, 4=S1_VH, 5=S2_HH, 6=S2_HV, 7=S2_VV, 8=S2_VH, 
    //       9=S3_HH, 10=S3_HV, 11=S3_VV, 12=S3_VH,
    //

    FILE *xml;
    struct tree *xml_aux;
    char *tmp_str[200], *str;
    int nc=0,n=0,nlmx=0;
    double *d;

    if ((xml = fopen(aux_file,"r")) == NULL) {
        die("Couldn't open xml file: \n",aux_file);
    }

    while (EOF != (ch=fgetc(XML_FILE))) {
        ++nc;
        if (ch == '\n') {
           ++n;
           if(nc > nlmx) nlmx = nc;
           nc=0;
        }
    }
    fclose(xml);

    if ((xml = fopen(aux_file,"r")) == NULL) {
        die("Couldn't open xml file: \n",aux_file);
    }
    xml_tree = (struct tree *)malloc(n*5*sizeof(struct tree));
    get_tree(xml,xml_aux,1);
    fclose(xml);

    search_tree(xml_aux,"/auxiliaryCalibration/calibrationParamsList/calibrationParams/elevationAntennaPattern/values/",tmp_str,3,3,mode);
    n = (int)str2double(tmp_str);    
    str = (char *)malloc(80*n*sizeof(char));
    d = (double *)malloc(n*2*sizeof(double));

    search_tree(xml_aux,"/auxiliaryCalibration/calibrationParamsList/calibrationParams/elevationAntennaPattern/values/",str,1,3,mode);
    str2dbs(d,str);
     

aaa
}
*/

int shift_burst(fcomplex *cbrst, int al_start, int lpb, int width, struct GMT_GRID *R, struct GMT_GRID *A, int bshift){

    int ii,jj,k;
    double ras[2];
    fcomplex *cbrst2;
    double incx,incy;
    int kr,ka;
    
    incx = R->header->inc[GMT_X];
    incy = R->header->inc[GMT_Y];

    cbrst2 = (fcomplex *) malloc(lpb*width*sizeof(fcomplex));
    for (ii=0;ii<lpb;ii++){
        for (jj=0;jj<width;jj++){
            k = ii*width+jj;
            cbrst2[k].r = cbrst[k].r;
            cbrst2[k].i = cbrst[k].i;
        }
    }

    for (ii=0;ii<lpb;ii++){
        for (jj=0;jj<width;jj++){
            k = ii*width+jj;
            if(floor((al_start+ii)/incy+0.5) < 0 || (int)floor((al_start+ii)/incy+0.5) >= R->header->ny){
                cbrst[k].r = 0;
		cbrst[k].i = 0;
            }
            else{
                kr = (int)floor(jj/incx+0.5)+R->header->nx*(int)floor((al_start+ii)/incy+0.5);
                ka = (int)floor(jj/incx+0.5)+A->header->nx*(int)floor((al_start+ii)/incy+0.5);
                ras[0] = (double)jj+R->data[kr];
                ras[1] = (double)ii+A->data[ka]  - (double)bshift;
                fbisinc(ras,cbrst2,lpb,width,&cbrst[k]);
            }
        }
    }
 
    free(cbrst2);
    return(1);
}


void sinc_one(double *, double *, double , double , double *);

void fbisinc (double *ras, fcomplex *s_in, int ydims, int xdims, fcomplex *sout)
{
        double dr, da, ns2=NS/2-1;
        double rdata[NS*NS], idata[NS*NS], cz[2];
        int i, j, k, kk;
        int i0, j0;
        int nclip;

        /* compute the residual offsets */
        nclip = 0;
        j0 = (int)floor(ras[0]);
        i0 = (int)floor(ras[1]);
        dr = ras[0] - (double)j0;
        da = ras[1] - (double)i0;
        if(dr < 0. || dr > 1. || da < 0. || da > 1) fprintf(stderr," dr or da out of bounds %f %f \n",dr,da);

        /* make sure all 4 corners are within the bounds of the slave array */

        if((i0-ns2) < 0 || (i0+ns2+1) >= ydims || (j0-ns2) < 0 || (j0+ns2+1) >= xdims) {
          sout->r = 0;
          sout->i = 0;
        }
        else {

        /* safe to do the interpolation */

        for (i=0; i<NS; i++){
                for (j=0; j<NS; j++){
                k = i*NS + j;
                kk = xdims*(i0-ns2+i) + (j0-ns2+j);
                rdata[k] = (double)s_in[kk].r;
                idata[k] = (double)s_in[kk].i;
                }
        }

        /* interpolate the real and imaginary data */

        sinc_one(rdata, idata, dr, da, cz);

        if((int)fabs(cz[0]) > I2MAX) nclip = nclip + 1;
        sout->r = (float)cz[0];
        if((int)fabs(cz[1]) > I2MAX) nclip = nclip + 1;
        sout->i = (float)cz[1];
        }
        //if(nclip > 0) fprintf(stderr," %d integers were clipped \n",nclip);
}


double sinc_kernel( double );

void sinc_one(double *rdata, double *idata, double x, double y, double *cz)
{
        int i, j, ij, ns2 = NS/2-1 ;
        double wx[NS], wy[NS];
        double arg, w, wsum, rsum, isum;

        for(i = 0; i < NS ; i++){
                arg = fabs(x + ns2 - i);
                wx[i] = sinc_kernel(arg);
                arg = fabs(y + ns2 - i);
                wy[i] = sinc_kernel(arg);
        }

        rsum = isum = wsum = 0.0;
        ij = 0;
        for (j = 0; j < NS; j++) {
                for (i = 0; i < NS; i++) {
                w = wx[i] * wy[j];
                rsum += rdata[ij+i] * w;
                isum += idata[ij+i] * w;
                wsum += w;
                }
                ij += NS;
        }
        if(wsum <= 0.0) printf(" error wsum is zero \n");
        cz[0] = rsum/wsum;
        cz[1] = isum/wsum;
}


#define PI 3.1415926535897932

double sinc_kernel (double x)
{
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



