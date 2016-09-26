/***************************************************************************
 * Creator:  Xiaohua(Eric) XU                                              *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  02/01/2016                                                    *
 ***************************************************************************/

/***************************************************************************
 * apply spectral diversity to estimate residual shift after               *
 * coregistration                                                          *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "xmlC.h"
#include "PRM.h"
#include "lib_functions.h"
#include "lib_defs.h"
#include "gmtsar.h"

typedef struct burst_bounds{
    int SL;
    int SC;
    int SH;
    int EL;
    int EC;
    int EH;
    int SLi;
    int SCi;
    int SHi;
    int ELi;
    int ECi;
    int EHi;
    int S;
    int E;
}burst_bounds;

char *USAGE = "\n Usage: spectral_diversity master_stem slave_stem bshfit filter\n"
              "\n Example: spectral_diversity S1A20150322_F1 S1A20150415_F1 0 gauss5x5\n"
              "\n Output: resitual_shift = 0.001234  [with a file ddphase]\n"
              "\n Note: make sure stem.SLCH stem.SLCL stem.BB exist \n";

int main(int argc, char **argv){

    FILE *MF, *MB, *SF, *SB, *BBM, *BBS, *FILTER ,*OUTP;
    char tmp_str[200];
    struct burst_bounds bbm[200],bbs[200];
    double tmp_d[200];
    int kkm,kks,splm,spls,spl,nbm,nbs,ii,jj,nlm,nls,ntls,ntlm,ntl;
    int bshift=0,llm,lls,nboff = 0;
    short *mf, *mb, *sf, *sb;
    float fmi,fmr,bmi,bmr,fsi,fsr,bsi,bsr;
    float i1,r1,i2,r2,*real,*imag,*filter,filtin,*freal,*fimag,filtdat,fsum=0.0,famp1,famp2;
    float *amp1,*amp2,*corr;
    double phase,isum=0.0,rsum=0.0; 
    int yarr, xarr, zz, *zz_r;
    double spec_sep,dta;

    if (argc < 5) die("",USAGE);

    // open corresponding files
    strcpy(tmp_str,argv[1]);
    strcat(tmp_str,".SLCH");
    if ((MF = fopen(tmp_str,"rb")) == NULL) die ("Couldn't open slch file: \n",tmp_str);
    strcpy(tmp_str,argv[1]);
    strcat(tmp_str,".SLCL");
    if ((MB = fopen(tmp_str,"rb")) == NULL) die ("Couldn't open slcl file: \n",tmp_str);
    strcpy(tmp_str,argv[2]);
    strcat(tmp_str,".SLCH");
    if ((SF = fopen(tmp_str,"rb")) == NULL) die ("Couldn't open slch file: \n",tmp_str);
    strcpy(tmp_str,argv[2]);
    strcat(tmp_str,".SLCL");
    if ((SB = fopen(tmp_str,"rb")) == NULL) die ("Couldn't open slcl file: \n",tmp_str);
    strcpy(tmp_str,argv[1]);
    strcat(tmp_str,".BB");
    if ((BBM = fopen(tmp_str,"rb")) == NULL) die ("Couldn't open bb file: \n",tmp_str);
    strcpy(tmp_str,argv[2]);
    strcat(tmp_str,".BB");
    if ((BBS = fopen(tmp_str,"rb")) == NULL) die ("Couldn't open bb file: \n",tmp_str);
    strcpy(tmp_str,argv[4]);
    if ((FILTER = fopen(tmp_str,"rb")) == NULL) die ("Couldn't open filter file: \n",tmp_str);
    bshift = (int)str2double(argv[3]);
    nboff = floor((double)bshift/1400.0+0.5);
    if (nboff != 0) {
        printf("Image has %d burst offset\n",nboff);
    }

    strcpy(tmp_str,"ddphase");
    if ((OUTP = fopen(tmp_str,"w")) == NULL) die ("Couldn't open output file: \n",tmp_str);
    //bshift = (int)str2double(argv[3]);
    bbm[0].SLi = bbm[0].SCi = bbm[0].SHi = bbm[0].ELi = bbm[0].ECi = bbm[0].EHi = -1;
    bbs[0].SLi = bbs[0].SCi = bbs[0].SHi = bbs[0].ELi = bbs[0].ECi = bbs[0].EHi = -1-bshift;
    bbm[0].S = bbm[0].E = -1;
    bbs[0].S = bbs[0].E = -1;

    // get parameters
    fgets(tmp_str,200*sizeof(char),BBM);
    tmp_str[strlen(tmp_str)-1]=' ';
    str2dbs(tmp_d,tmp_str);
    nbm = (int)tmp_d[0];
    splm = (int)tmp_d[1];
    nlm = 0;
    ntlm = 0;

    for (ii=1;ii<=nbm;ii++){
        fgets(tmp_str,200*sizeof(char),BBM);
        //tmp_str[strlen(tmp_str)-1]=' ';
        str2dbs(tmp_d,tmp_str);
        bbm[ii].SL=(int)tmp_d[0];
        bbm[ii].SC=(int)tmp_d[1];
        bbm[ii].SH=(int)tmp_d[2];
        bbm[ii].EL=(int)tmp_d[3];
        bbm[ii].EC=(int)tmp_d[4];
        bbm[ii].EH=(int)tmp_d[5];
        //fprintf(stderr,"%d %d %d %d %d %d\n",bbm[ii].SL,bbm[ii].SC,bbm[ii].SH,bbm[ii].EL,bbm[ii].EC,bbm[ii].EH);
        bbm[ii].SLi = bbm[ii-1].ELi + 1;
        bbm[ii].SCi = bbm[ii-1].ECi + 1;
        bbm[ii].SHi = bbm[ii-1].EHi + 1;
        bbm[ii].ELi = bbm[ii].SLi + bbm[ii].EL - bbm[ii].SL;
        bbm[ii].ECi = bbm[ii].SCi + bbm[ii].EC - bbm[ii].SC;
        bbm[ii].EHi = bbm[ii].SHi + bbm[ii].EH - bbm[ii].SH;
        bbm[ii].S = bbm[ii-1].E + 1;
        bbm[ii].E = bbm[ii].S + bbm[ii].EH - (bbm[ii].EL+1);
        //fprintf(stderr,"%d %d %d %d %d %d %d %d\n",bbm[ii].SLi,bbm[ii].SCi,bbm[ii].SHi,bbm[ii].ELi,bbm[ii].ECi,bbm[ii].EHi,bbm[ii].S,bbm[ii].E);
        nlm = nlm + bbm[ii].EH - (bbm[ii].EL+1)+1;
        ntlm = ntlm + bbm[ii].EC - bbm[ii].SC + 1;
    }
    //fprintf(stderr,"Slave Image Size %d x %d \n",nlm,splm);

    fgets(tmp_str,200*sizeof(char),BBS);
    tmp_str[strlen(tmp_str)-1]=' ';
    str2dbs(tmp_d,tmp_str);
    nbs = (int)tmp_d[0];
    spls = (int)tmp_d[1];
    spec_sep = tmp_d[2];
    dta = tmp_d[3];
    nls = 0;
    ntls = 0;

    for (ii=1;ii<=nbs;ii++){
        fgets(tmp_str,200*sizeof(char),BBS);
        //tmp_str[strlen(tmp_str)-1]=' ';
        str2dbs(tmp_d,tmp_str);
        bbs[ii].SL=(int)tmp_d[0];
        bbs[ii].SC=(int)tmp_d[1];
        bbs[ii].SH=(int)tmp_d[2];
        bbs[ii].EL=(int)tmp_d[3];
        bbs[ii].EC=(int)tmp_d[4];
        bbs[ii].EH=(int)tmp_d[5];
        //fprintf(stderr,"%d %d %d %d %d %d\n",bbs[ii].SL,bbs[ii].SC,bbs[ii].SH,bbs[ii].EL,bbs[ii].EC,bbs[ii].EH);
        bbs[ii].SLi = bbs[ii-1].ELi + 1;
        bbs[ii].SCi = bbs[ii-1].ECi + 1;
        bbs[ii].SHi = bbs[ii-1].EHi + 1;
        bbs[ii].ELi = bbs[ii].SLi + bbs[ii].EL - bbs[ii].SL;
        bbs[ii].ECi = bbs[ii].SCi + bbs[ii].EC - bbs[ii].SC;
        bbs[ii].EHi = bbs[ii].SHi + bbs[ii].EH - bbs[ii].SH;
        bbs[ii].S = bbs[ii-1].E + 1;
        bbs[ii].E = bbs[ii].S + bbs[ii].EH - (bbs[ii].EL+1);
        //fprintf(stderr,"%d %d %d %d %d %d %d %d\n",bbs[ii].SLi,bbs[ii].SCi,bbs[ii].SHi,bbs[ii].ELi,bbs[ii].ECi,bbs[ii].EHi,bbs[ii].S,bbs[ii].E);
        nls = nls + (bbs[ii].SH-1) - bbs[ii].SL + 1;
        ntls = ntls + bbs[ii].EC - bbs[ii].SC + 1;
    }    

    //fprintf(stderr,"Slave Image Size %d x %d \n",nls,spls);
    

    // malloc memory for images
    mf = (short *)malloc(nlm*splm*2*sizeof(short));
    mb = (short *)malloc(nlm*splm*2*sizeof(short));
    sf = (short *)malloc(nls*spls*2*sizeof(short));
    sb = (short *)malloc(nls*spls*2*sizeof(short));
    fread(mf,nlm*splm*2,sizeof(short),MF);
    fread(mb,nlm*splm*2,sizeof(short),MB);
    fread(sf,nls*spls*2,sizeof(short),SF);
    fread(sb,nls*spls*2,sizeof(short),SB);

    kkm = 1; 
    kks = kkm+nboff;
    ntl = ntlm;
    if (ntl>=ntls) ntl = ntls;

    if (bshift > ntl) {
        fprintf(stderr,"Images does not overlap, returning 0 to res_shift\n");
        printf("residual_phase =  %.6f\n   isum = %.2g   rsum = %.2g\n",0.0,0.0,0.0);
        printf("spectral_spectrationXdta = %.6f\n",spec_sep*dta);
        printf("residual_shift = %.12f\n",0.0);
        
    }


    spl = splm;
    if (spl>=spls) spl = spls;
    real = (float *)malloc(ntl*spl*sizeof(float));
    imag = (float *)malloc(ntl*spl*sizeof(float));
    freal = (float *)malloc(ntl*spl*sizeof(float));
    fimag = (float *)malloc(ntl*spl*sizeof(float));
    zz_r = (int *)malloc(ntl*sizeof(int));
    amp1 = (float *)malloc(ntl*spl*sizeof(float));
    amp2 = (float *)malloc(ntl*spl*sizeof(float));
    
    corr = (float *)malloc(ntl*spl*sizeof(float));
    // compute sum real and sum imagenary
    //fprintf(stderr,"Some pars: ntl: %d %d %d, spl: %d %d %d\n",ntlm,ntls,ntl,splm,spls,spl);
    zz = 0;
    while(kkm < 1 || kks < 1){
        kks++;
        kkm++;
    }
    if(kkm!=kks) printf("starting bursts are %d for master and %d for slave\n",kkm,kks);
    //printf("Working on burst %d (master, zz = %d)...\n",kkm,zz);
    for(ii=0;ii<ntl;ii++){
        if(ii>=bbm[kkm].ELi+1 && ii<=bbm[kkm].EHi && ii>=bbs[kks].ELi+1 && ii<=bbs[kks].EHi){
            //fprintf(stderr,"working on Line %d...\n",ii);
            llm = ii - (bbm[kkm].ELi+1) + bbm[kkm].S;
            lls = ii - (bbs[kks].ELi+1) + bbs[kks].S;
            //if (llm%100 == 0) fprintf(stderr,"master_line %d, slave_line %d\n",llm,lls);
            for(jj=0;jj<spl;jj++){
                fmr = (float)mf[(llm*splm+jj)*2];
                fmi = (float)mf[(llm*splm+jj)*2+1];
                bmr = (float)mb[(llm*splm+jj)*2];
                bmi = (float)mb[(llm*splm+jj)*2+1];
                fsr = (float)sf[(lls*spls+jj)*2];
                fsi = (float)sf[(lls*spls+jj)*2+1];
                bsr = (float)sb[(lls*spls+jj)*2];
                bsi = (float)sb[(lls*spls+jj)*2+1];               
                
                r1 = fmr*fsr+fmi*fsi;
                i1 = fmi*fsr-fmr*fsi;
                r2 = bmr*bsr+bmi*bsi;
                i2 = bmi*bsr-bmr*bsi;

                amp1[zz*spl+jj] = r1*r1+i1*i1;
                amp2[zz*spl+jj] = r2*r2+i2*i2;
                 
                real[zz*spl+jj] = r1*r2+i1*i2;
                imag[zz*spl+jj] = i1*r2-r1*i2;
                //phase = sqrt(real*real+imag*imag)/sqrt((r1*r1+i1*i1)*(r2*r2+i2*i2));
                //phase = atan2(imag[zz*spl+jj],real[zz*spl+jj]);
            }
            zz_r[zz] = ii;
            zz++;
        }
        if(ii>bbm[kkm].EHi && ii>bbs[kks].EHi) {
            kkm++;kks++;
            //printf("Working on burst %d (master, zz = %d)...\n",kkm,zz);
            //fprintf(stderr,"Computing next burst %d...\n",kk);
        }

    }
    // filter the phase
    if (fscanf(FILTER,"%d%d", &xarr, &yarr) != 2 || xarr < 1 || yarr < 1 || (xarr & 1) == 0 || (yarr & 1) == 0) die("filter incomplete","");
    if (( filter = (float *)malloc(sizeof(float)*xarr*yarr)) == NULL) die("memory allocation","");
    for (ii=0;ii<xarr*yarr;ii++) {
        if (fscanf(FILTER,"%f",&filtin) == EOF) die("filter incomplete","");
        filter[ii] = filtin;
        fsum += filtin; 
    }
    
    for(ii=0;ii<zz;ii++) {
        for(jj=0;jj<spl;jj++) {
            conv2d(real,&zz,&spl,filter,&xarr,&yarr,&filtdat,&ii,&jj,&fsum);
            freal[ii*spl+jj]=filtdat;
            conv2d(imag,&zz,&spl,filter,&xarr,&yarr,&filtdat,&ii,&jj,&fsum);
            fimag[ii*spl+jj]=filtdat;
            conv2d(amp1,&zz,&spl,filter,&xarr,&yarr,&filtdat,&ii,&jj,&fsum);
            famp1 = filtdat;
            
            conv2d(amp2,&zz,&spl,filter,&xarr,&yarr,&filtdat,&ii,&jj,&fsum);
            famp2 = filtdat;
            
            corr[ii*spl+jj] = sqrt((freal[ii*spl+jj]*freal[ii*spl+jj]+fimag[ii*spl+jj]*fimag[ii*spl+jj])/(famp1*famp2));

            if (corr[ii*spl+jj] > 0.6) {
                rsum+=freal[ii*spl+jj];
            }
            if (corr[ii*spl+jj] > 0.6) {
                isum+=fimag[ii*spl+jj];
            }
        }
    }
    printf("Image analyzed %dx%d...\n",spl,zz);
    
    for(ii=0;ii<zz;ii+=10) {
        for(jj=0;jj<spl;jj++) {
            if(corr[ii*spl+jj]>0.3) {
                //fprintf(OUTP,"%d\t%d\t%.9f\n",jj,zz_r[ii],(float)(sqrt((jj-10000)*(jj-10000)+(zz_r[ii]-5000)*(zz_r[ii]-5000)))/5.0e5);
                phase = atan2(fimag[ii*spl+jj],freal[ii*spl+jj]);
                fprintf(OUTP,"%d\t%d\t%.9f\t%.9f\n",jj,zz_r[ii]+bshift,phase,corr[ii*spl+jj]);
                //fwrite(&phase,1,sizeof(double),OUTP);
            }
        }
    }

    phase = atan2(isum,rsum);
    printf("residual_phase =  %.6f\n   isum = %.2g   rsum = %.2g\n",phase,isum,rsum);
    printf("spectral_spectrationXdta = %.6f\n",spec_sep*dta);
    printf("residual_shift = %.12f\n",phase/(2*M_PI*spec_sep*dta));

    // free memory and close corresponding files
    free(mf);
    free(mb);
    free(sf);
    free(sb);
    free(zz_r);
    free(amp1);
    free(amp2);
    free(corr);
    free(freal);
    free(fimag);
    free(real);
    free(imag);
    fclose(MF);
    fclose(MB);
    fclose(SF);
    fclose(SB);
    fclose(BBM);
    fclose(BBS);
    fclose(FILTER);
    fclose(OUTP);

    return(1);

}


