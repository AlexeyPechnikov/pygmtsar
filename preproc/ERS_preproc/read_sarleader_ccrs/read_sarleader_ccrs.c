/* This program is supposed to read raw SAR data tapes 	*/
/* 							*/
/* 		R. Mellors 	July 1997		*/
/*				IGPP_SIO		*/
/*							*/
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include "../include/SARtape.h"

int main(argc,argv)
int	argc;
char	**argv;
{
char	*filename,string[5];
int	nitems,i,num_data_points;
FILE	*file1,*outfile;
struct	SAR_info sar;
struct	sarleader_binary slfsb,slvsb,slplb;

if (argc < 2) {
	fprintf(stderr,"Usage: read_SAR_tape SARLEADER \n");
	exit(1);
	}

filename = argv[1];

file1 = fopen(filename,"r");
outfile = fopen("ldrfile.log","w");

fprintf(outfile,"listing header values to %s\n",filename);

/* allocate memory */
sar.fixseg = (struct sarleader_fdr_fixseg *) malloc(sizeof(struct sarleader_fdr_fixseg));
sar.varseg = (struct sarleader_fdr_varseg *) malloc(sizeof(struct sarleader_fdr_varseg));
sar.dss = (struct sarleader_dss *) malloc(sizeof(struct sarleader_dss));
sar.platform = (struct platform *) malloc(sizeof(struct platform));

/* read the file */
nitems = fread(&slfsb,sizeof(struct sarleader_binary),1,file1);
fprintf(outfile," read %d bytes %d items from file %s : position %ld\n",(sizeof(struct sarleader_binary)),nitems,filename,ftell(file1));
fprintf(outfile,SARLEADER_FDR_BINARY_WCS,SARLEADER_FDR_BINARY_RVL(&slfsb));

fscanf(file1,SARLEADER_FDR_FIXSEG_RCS,SARLEADER_FDR_FIXSEG_RVL(sar.fixseg));
fprintf(outfile,SARLEADER_FDR_FIXSEG_WCS,SARLEADER_FDR_FIXSEG_RVL(sar.fixseg));

fscanf(file1,SARLEADER_FDR_VARSEG_RCS,SARLEADER_FDR_VARSEG_RVL(sar.varseg));
fprintf(outfile,SARLEADER_FDR_VARSEG_WCS,SARLEADER_FDR_VARSEG_RVL(sar.varseg));

nitems = fread(&slvsb,sizeof(struct sarleader_binary),1,file1);
fprintf(outfile," read %d bytes %d items from file %s : position %ld\n",(sizeof(struct sarleader_binary)),nitems,filename,ftell(file1));
fprintf(outfile,SARLEADER_FDR_BINARY_WCS,SARLEADER_FDR_BINARY_RVL(&slfsb));

fscanf(file1,SARLEADER_DSS_RCS,SARLEADER_DSS_RVL(sar.dss));
fprintf(outfile,SARLEADER_DSS_WCS,SARLEADER_DSS_RVL(sar.dss));

/* why do I need to skip 1664 bytes ahead; what am I missing ? */
fseek(file1,1664,SEEK_CUR);

nitems = fread(&slplb,sizeof(struct sarleader_binary),1,file1);
fprintf(outfile," read %d bytes %d items from file %s : position %ld\n",(sizeof(struct sarleader_binary)),nitems,filename,ftell(file1));
fprintf(outfile,SARLEADER_FDR_BINARY_WCS,SARLEADER_FDR_BINARY_RVL(&slplb));

fscanf(file1,PLATFORM_RCS,PLATFORM_RVL(sar.platform));
fprintf(outfile,PLATFORM_WCS,PLATFORM_RVL(sar.platform));

sar.position = (struct position_vector *) malloc(sizeof(struct position_vector));

sscanf(sar.platform->num_data_points," %1c",string);
for (i=0;i<(atoi(string));i++)
	{
	fscanf(file1,POSITION_VECTOR_RCS,POSITION_VECTOR_RVL(sar.position));
	fprintf(outfile,POSITION_VECTOR_WCS,POSITION_VECTOR_RVL(sar.position));
	}

/* now make the PRM file */
make_prm(sar);
}
