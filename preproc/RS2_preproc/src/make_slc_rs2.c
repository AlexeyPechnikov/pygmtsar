/***************************************************************************
 * Creator:  Xiaohua(Eric) XU                                              *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  04/06/2015                                                    *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * DATE                                                                    *
 *                                                                         *
 ***************************************************************************/

#include "tiffio.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "PRM.h"
#include "lib_functions.h"
#include "stateV.h"
#include "xmlC.h"
#include "lib_defs.h"

int pop_prm(struct PRM *, tree *, char *);
int pop_led(tree *,state_vector *);
int write_orb(state_vector *sv, FILE *fp, int);
int write_slc(TIFF *, FILE *, char *);


char *USAGE = "\n\nUsage: make_slc_rs2 name_of_xml_file name_of_tiff_file name_output\n"
"\nExample: make_slc_rs2 product.xml imagery_HH.tif RS220110515\n"
"\nOutput: RS220110515.SLC RS220110515.PRM RS220110515.LED\n";

int main(int argc, char **argv){
    
    FILE *XML_FILE,*OUTPUT_PRM,*OUTPUT_SLC,*OUTPUT_LED;
    TIFF *TIFF_FILE;
    char tmp_str[200];
    struct PRM prm;
    tree *xml_tree;
    state_vector sv[200];
    int ch, n=0, nc=0, nlmx=0;
    
    if (argc < 4) die (USAGE,"");
    
    if ((XML_FILE = fopen(argv[1],"rb")) == NULL) die("Couldn't open xml file: \n",argv[1]);

    // find the number of lines and the maximum line length of the xml file
    while (EOF != (ch=fgetc(XML_FILE))) {
        ++nc;
        if (ch == '\n') {
           ++n;
           if(nc > nlmx) nlmx = nc;
           nc=0;
        }
    }
    fprintf(stderr,"%d %d \n",n,nlmx);
    xml_tree = (struct tree *)malloc(n*5*sizeof(struct tree));
    fclose(XML_FILE);
    
    // generate the xml tree
    if ((XML_FILE = fopen(argv[1],"r")) == NULL) die("Couldn't open xml file: \n",argv[1]);
    get_tree(XML_FILE,xml_tree,1);
    fclose(XML_FILE);

    //show_tree(xml_tree,0,0);
    
    // initiate the prm
    null_sio_struct(&prm);

    // generate the PRM file
    pop_prm(&prm,xml_tree,argv[3]);
    strcpy(tmp_str,argv[3]);
    strcat(tmp_str,".PRM");
    if ((OUTPUT_PRM = fopen(tmp_str,"w")) == NULL) die ("Couldn't open prm file: \n",tmp_str);
    put_sio_struct(prm, OUTPUT_PRM);
    fclose(OUTPUT_PRM);
    
    // generate the LED file
    n = pop_led(xml_tree,sv);
    
    strcpy(tmp_str,argv[3]);
    strcat(tmp_str,".LED");
    if ((OUTPUT_LED = fopen(tmp_str,"w")) == NULL) die ("Couldn't open led file: \n",tmp_str);
    write_orb(sv,OUTPUT_LED,n);
    fclose(OUTPUT_LED);
    
    // generate the SLC file
    TIFFSetWarningHandler(NULL);
    if ((TIFF_FILE = TIFFOpen(argv[2], "r")) == NULL) die ("Couldn't open tiff file: \n",argv[2]);
    
    strcpy(tmp_str,argv[3]);
    strcat(tmp_str,".SLC");
    if ((OUTPUT_SLC = fopen(tmp_str,"wb")) == NULL) die ("Couldn't open slc file: \n",tmp_str);
    write_slc(TIFF_FILE,OUTPUT_SLC,prm.orbdir);
    
    TIFFClose(TIFF_FILE);
    fclose(OUTPUT_SLC);
}

int write_slc(TIFF *tif, FILE *slc, char *orbdir){

    uint32 width,height,widthi;
    int i,j;
    uint16 s=0,nsamples;
    uint16 *buf;
    short *tmp;
    
    // get the width and the height of the file, make width dividable by 4
    TIFFGetField(tif,TIFFTAG_IMAGEWIDTH,&widthi);
    TIFFGetField(tif,TIFFTAG_IMAGELENGTH,&height);
    width = widthi - widthi%4;
    //printf("%d %d \n",width,height);
    
    buf = (uint16 *)_TIFFmalloc(TIFFScanlineSize(tif));
    tmp = (short *)malloc(width*2*sizeof(short));
    printf("Writing SLC..Image Size: %d X %d...\n",width,height);
    
    TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &nsamples);
    //printf("Total Number of Samples: %d\n",nsamples);
    // For Ascending Condition...
    if (strcmp(orbdir,"A") == 0){
        printf("Fliping upside down for Ascending Image...\n");
        for (i=height-1;i>=0;i--){
            TIFFReadScanline(tif, buf, i, s);
            for (j=0;j<width*2;j++){
                tmp[j] = (short)buf[j];
            }
            fwrite(tmp,sizeof(short),width*2,slc);
        }
    }
    // For Descending Condition...
    else{
        printf("Fliping leftside right for Descending Image...\n");
        for (i=0;i<height;i++){
            TIFFReadScanline(tif, buf, i, s);
            for (j=0;j<width;j++){
                tmp[(width-1-j)*2] = (short)buf[j*2];
                tmp[(width-1-j)*2+1] = (short)buf[j*2+1];
            }
            fwrite(tmp,sizeof(short),width*2,slc);
        }
    }
    _TIFFfree(buf);
    free(tmp);
    return(1);
}

int write_orb(state_vector *sv, FILE *fp, int n){
    int i;
    double dt;
    
    dt = trunc((sv[1].sec)*1000.0)/1000.0-trunc((sv[0].sec)*1000.0)/1000.0;
    if(n<=1) return(-1);
    fprintf(fp,"%d %d %d %.3lf %lf \n",n,sv[0].yr,sv[0].jd,sv[0].sec,dt);
    for(i=0;i<n;i++){
        fprintf(fp,"%d %d %.3lf %.6lf %.6lf %.6lf %.8lf %.8lf %.8lf \n",sv[i].yr,sv[i].jd,sv[i].sec,sv[i].x,sv[i].y,sv[i].z,sv[i].vx,sv[i].vy,sv[i].vz);
    }
    return(1);
}

int pop_led(tree *xml_tree,state_vector *sv){
    int i,ct;
    char tmp_c[200];
    double tmp_d;
    
    //search_tree(xml_tree,"/product/generalAnnotation/orbitList/",tmp_c,3,0,1);
    //count = 1;//(int)str2double(tmp_c);
    
    i = 1;
    while (1){
        ct = search_tree(xml_tree,"/product/sourceAttributes/orbitAndAttitude/orbitInformation/stateVector/",tmp_c,1,5,i);
        search_tree(xml_tree,"/product/sourceAttributes/orbitAndAttitude/orbitInformation/stateVector/timeStamp/",tmp_c,2,5,i);
        tmp_d = str2double(tmp_c);
        search_tree(xml_tree,"/product/sourceAttributes/orbitAndAttitude/orbitInformation/stateVector/timeStamp/",tmp_c,1,5,i);
        tmp_c[4] = '\0';
        sv[i-1].yr = (int)(str2double(tmp_c));
        sv[i-1].jd = (int)(tmp_d - trunc(tmp_d/1000.0)*1000.0);
        sv[i-1].sec = (tmp_d - trunc(tmp_d))*86400;
        search_tree(xml_tree,"/product/sourceAttributes/orbitAndAttitude/orbitInformation/stateVector/xPosition/",tmp_c,1,5,i);
        sv[i-1].x = str2double(tmp_c);
        search_tree(xml_tree,"/product/sourceAttributes/orbitAndAttitude/orbitInformation/stateVector/yPosition/",tmp_c,1,5,i);
        sv[i-1].y = str2double(tmp_c);
        search_tree(xml_tree,"/product/sourceAttributes/orbitAndAttitude/orbitInformation/stateVector/zPosition/",tmp_c,1,5,i);
        sv[i-1].z = str2double(tmp_c);
        search_tree(xml_tree,"/product/sourceAttributes/orbitAndAttitude/orbitInformation/stateVector/xVelocity/",tmp_c,1,5,i);
        sv[i-1].vx = str2double(tmp_c);
        search_tree(xml_tree,"/product/sourceAttributes/orbitAndAttitude/orbitInformation/stateVector/yVelocity/",tmp_c,1,5,i);
        sv[i-1].vy = str2double(tmp_c);
        search_tree(xml_tree,"/product/sourceAttributes/orbitAndAttitude/orbitInformation/stateVector/zVelocity/",tmp_c,1,5,i);
        sv[i-1].vz = str2double(tmp_c);
        //fprintf(stderr,"Hahahahahaahah...\n");
        if (xml_tree[ct].sibr == -1 || strcmp(xml_tree[xml_tree[ct].sibr].name,"stateVector") != 0) {
            break;
        }
        else{
            i = i+1;
        }
        
    }
    printf("%d Lines Written for Orbit...\n",i);
    return(i);
}

int pop_prm(struct PRM *prm,tree *xml_tree,char *file_name){
    char tmp_c[2000];
    double tmp_d;
    int tmp_i;
    double c_speed = 299792458.0;
    
    // define some of the variables
    prm->first_line = 1;
    prm->st_rng_bin = 1;
    search_tree(xml_tree,"/product/imageGenerationParameters/sarProcessingInformation/numberOfRangeLooks/",tmp_c,1,0,1);
    fprintf(stderr,"%s \n",tmp_c);
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
    
    search_tree(xml_tree,"/product/imageAttributes/rasterAttributes/sampledPixelSpacing/",tmp_c,1,0,1);
    prm->fs = c_speed/2.0/str2double(tmp_c);//rng_samp_rate
    prm->SC_identity = 9; /* (1)-ERS1 (2)-ERS2 (3)-Radarsat (4)-Envisat (5)-ALOS (6)-  (7)-TSX (8)-CSK (9)-RS2 (10) Sentinel-1a*/
    
    search_tree(xml_tree,"/product/sourceAttributes/radarParameters/radarCenterFrequency/",tmp_c,1,0,1);
    prm->lambda = c_speed/str2double(tmp_c);
    
    search_tree(xml_tree,"/product/sourceAttributes/radarParameters/pulseLength/",tmp_c,1,0,1);
    tmp_d = str2double(tmp_c);
    search_tree(xml_tree,"/product/sourceAttributes/radarParameters/pulseBandwidth/",tmp_c,1,0,1);
    prm->chirp_slope = str2double(tmp_c)/tmp_d;
    prm->pulsedur = tmp_d;
    
    //search_tree(xml_tree,"/product/qualityInformation/qualityDataList/qualityData/imageQuality/imageStatistics/outputDataMean/re/",tmp_c,1,0,1);
    prm->xmi = 0; //str2double(tmp_c); //I_mean
    
    //search_tree(xml_tree,"/product/qualityInformation/qualityDataList/qualityData/imageQuality/imageStatistics/outputDataMean/im/",tmp_c,1,0,1);
    prm->xmq = 0; //str2double(tmp_c); //Q_mean
    
    search_tree(xml_tree,"/product/sourceAttributes/radarParameters/pulseRepetitionFrequency/",tmp_c,1,0,1);
    prm->prf = str2double(tmp_c)/2.0;
    search_tree(xml_tree,"/product/imageGenerationParameters/slantRangeToGroundRange/slantRangeTimeToFirstRangeSample/",tmp_c,1,0,1);
    prm->near_range = str2double(tmp_c)*c_speed/2.0;
    prm->ra = 6378137.00; //equatorial_radius
    prm->rc = 6356752.31; //polar_radius
    
    search_tree(xml_tree,"/product/sourceAttributes/orbitAndAttitude/orbitInformation/passDirection/",tmp_c,1,0,1);
    strasign(prm->orbdir,tmp_c,0,0);
    
    search_tree(xml_tree,"/product/sourceAttributes/radarParameters/antennaPointing/",tmp_c,1,0,1);
    strasign(prm->lookdir,tmp_c,0,0);

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
    if (strcmp(prm->orbdir,"D") == 0){
        search_tree(xml_tree,"/product/imageGenerationParameters/sarProcessingInformation/zeroDopplerTimeFirstLine/",tmp_c,2,0,1);
    }
    else{
        search_tree(xml_tree,"/product/imageGenerationParameters/sarProcessingInformation/zeroDopplerTimeLastLine/",tmp_c,2,0,1);
    }
    prm->clock_start = str2double(tmp_c);
    search_tree(xml_tree,"/product/imageGenerationParameters/sarProcessingInformation/zeroDopplerTimeLastLine/",tmp_c,1,0,1);
    tmp_c[4] = '\0';
    prm->SC_clock_start = prm->clock_start + 1000.*str2double(tmp_c);
    
    strasign(prm->iqflip,"n",0,0); //Flip_iq
    strasign(prm->deskew,"n",0,0); //deskew
    strasign(prm->offset_video,"n",0,0);
    
    search_tree(xml_tree,"/product/imageAttributes/rasterAttributes/numberOfSamplesPerLine/",tmp_c,1,0,1);
    //tmp_i = (int)str2double(tmp_c) - (int)str2double(tmp_c)%4;
    prm->bytes_per_line = (int)str2double(tmp_c)*4;//tmp_i*4;
    prm->good_bytes = prm->bytes_per_line;
    prm->caltone = 0.0;
    prm->pctbwaz = 0.0; //rm_az_band
    prm->pctbw = 0.2; //rm_rng_band
    prm->rhww = 1.0; //rng_spec_wgt
    strasign(prm->srm,"0",0,0); //scnd_rng_mig
    prm->az_res = 0.0;
    //prm.antenna_side = -1;
    prm->fdd1 = 0.0;
    prm->fddd1 = 0.0;
    
    search_tree(xml_tree,"/product/imageAttributes/rasterAttributes/numberOfLines/",tmp_c,1,0,1);
    tmp_i = (int)str2double(tmp_c);
    prm->num_lines = tmp_i - tmp_i%4;
    
    //search_tree(xml_tree,"/product/adsHeader/stopTime/",tmp_c,2,0,1);
    prm->SC_clock_stop = prm->SC_clock_start + prm->num_lines/prm->prf/86400;
    prm->clock_stop = prm->clock_start + prm->num_lines/prm->prf/86400;
    
    prm->nrows = prm->num_lines;
    prm->num_valid_az = prm->num_lines;
    prm->num_patches = 1;
    prm->num_rng_bins = prm->bytes_per_line/4;
    prm->chirp_ext	= 0;
    
    printf("PRM set for Image File...\n");
    return(1);
}

