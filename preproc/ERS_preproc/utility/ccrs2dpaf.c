/* program to transform CCRS format to DPAF format of ERS data
 method: CCRS has 12060 bytes per line. DPAF has 11644. Read 12060 bytes from CCRS file and print 11644 bytes to new file */

#include "../include/soi.h"
#include<string.h>

main(argc,argv)
int argc;
char *argv[];
{
        char *ifile, *ofile;    /* input file and output file*/
        FILE *fopen(), *fpi, *fpo;
        unsigned char *indata;
        int ccrs_length = 12060, dpaf_length = 11644;
        int n;

        if(argc<3){
           fprintf(stderr,"Usage: %s ccrs_raw out_file \n",argv[0]);
           exit(-1);
        }

        if((fpi=fopen(argv[1],"r")) == NULL){
           fprintf(stderr,"Can't open file %s \n",argv[1]);
           exit(-1);
        }

	if((n = fseek(fpi,ccrs_length,0)) != 0){
          perror(argv[0]);
          exit(-1);
        }

        if((fpo=fopen(argv[2],"w")) == NULL){
           fprintf(stderr,"Can't open output file %s \n",argv[2]);
           exit(-1);
        }

/* allocate memory */
        if((indata = (unsigned char *) malloc(ccrs_length*sizeof(unsigned char)))==NULL){
           fprintf(stderr,"Sorry, couldn't allocate memory for input indata.\n");
           exit(-1);
        }

/* read data */
        while(n=(fread((void *)indata,sizeof(unsigned char),ccrs_length,fpi))==ccrs_length){
		memcpy(indata+210,indata+200,8);
                if(n=(fwrite((void *)indata,sizeof(unsigned char),dpaf_length,fpo))!=dpaf_length){
                        fprintf(stderr,"Problem writing data \n");
                }
        }

        fclose(fpi);
        fclose(fpo);
}
