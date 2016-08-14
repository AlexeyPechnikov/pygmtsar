/* program to read sample window start time,
pulse repetition interval, and count lines
and image format counter
in CEOS SAR data file 

	rjm 	UCSD/SDSU Sept 97
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "SARtape.h"
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include "data_param.h"

#define SWAP_2(x) ( (((x) & 0xff) << 8) | ((unsigned short)(x) >> 8) )
#define SWAP_4(x) ( ((x) << 24) | \
                (((x) << 8) & 0x00ff0000) | \
                (((x) >> 8) & 0x0000ff00) | \
                ((x) >> 24) )
#define FIX_SHORT(x) (*(unsigned short *)&(x) = SWAP_2(*(unsigned short *)&(x)))
#define FIX_INT(x)   (*(unsigned int *)&(x)   = SWAP_4(*(unsigned int *)&(x)))
#define FIX_FLOAT(x) FIX_INT(x)

#define MIN_PRF	1640.0
#define MAX_PRF 1720.0
#define SEC_PER_PRI_COUNT 210.94e-09
#define LINELENGTH 12060
#define PREFIX 412
#define SUFFIX 416
#define IFC_OFF 200
#define SOL 299792456.0

int is_big_endian_(void);
int is_big_endian__(void);

int main(argc,argv)
int	argc;
char	*argv[];
{

char	*data,logfilename[255];
//int	file_size,year;
int	prior_pri_dn;
FILE	*indata,*logfile;
struct   sarleader_binary slb;
struct   sardata_rec sdr;
struct   SAR_info sar;
struct	 lineparam info;

int	endian;
int	logflag;
int	nlines,linelength,prefix,suffix;
int	num_patches,good_bytes_per_line,first_sample;
int	iwrite = 0;

//char *iptr;
int  ncnt, print_start;
unsigned short *icu_time1, *icu_time2;
double icu_time, icu_time_old;

double SC_clock_start,SC_clock_stop;
double clock_start,clock_stop;

struct lineparam {
int	ifc;
unsigned short	swst_dn;
unsigned short	pri_dn;
double 		pri;
double		swst;
};

double calc_pri();
double calc_swst();

logflag = 0;

/*------------------------*/
/* check endian           */
/*------------------------*/

        endian = is_big_endian_();

if (argc < 2) 
	{
	fprintf(stderr,"Usage: %s datafile\n",argv[0]);
	exit(1);
	}

indata = fopen(argv[1],"r");

if ((argc > 3) && ((strcmp(argv[1],"-log")==0))) {
	strcpy(logfilename,"dataheader.log");
	logfile = fopen(logfilename,"w");
	logflag = 1;
	}

if (indata == NULL) {
	fprintf(stderr,"error opening data file\n");
	exit(1);
	}

/* read top of file */
(void)fread(&slb,sizeof(struct sarleader_binary),1,indata);

sar.fixseg = (struct sarleader_fdr_fixseg *) malloc(sizeof(struct sarleader_fdr_fixseg));

/* read nominal parameters */
fscanf(indata,SARLEADER_FDR_FIXSEG_RCS,SARLEADER_FDR_FIXSEG_RVL(sar.fixseg));

sar.dataheader = (struct sardata_header *)malloc(sizeof(struct sardata_header));
fscanf(indata, SARDATA_HEADER_RCS, SARDATA_HEADER_RVL(sar.dataheader));

/* reset file pointer and read first line of data */

fseek(indata,LINELENGTH,0);
(void)fread(&slb,sizeof(struct sarleader_binary),1,indata);
(void)fread(&sdr,sizeof(struct sardata_rec),1,indata);

if (logflag) {
	fprintf(logfile,SARLEADER_FDR_FIXSEG_WCS,SARLEADER_FDR_FIXSEG_RVL(sar.fixseg));
	fprintf(logfile,SARDATA_HEADER_WCS, SARDATA_HEADER_RVL(sar.dataheader));
	/*fprintf(logfile,SARDATA_REC_WCS, SARDATA_REC_RVL(&sdr));*/
	}

/* used in SAR processor - each patch is 2800 long */
nlines = atoi(sar.dataheader->n_records);
num_patches = nlines/2800;

/* double-check length of line */
linelength = atoi(sar.dataheader->record_length);

/* check data preifx and suffic; add 12 to prefix for top of file */
prefix = atoi(sar.dataheader->n_bytes_prefix_data) + 12;
suffix = atoi(sar.dataheader->n_bytes_suffix_data);

/* SAR processor needs to know where the first sample is
so skip the dataheader (specified in equivalent samples - 2 bytes per sample)
and subtract the data suffix (in bytes) */ 

first_sample = prefix/2;
good_bytes_per_line = linelength - suffix;

/* set the spacecraft time */

if(endian == -1) {
   FIX_INT(sdr.year);
   FIX_INT(sdr.day_of_year);
   FIX_INT(sdr.msecs_of_day);
}

SC_clock_start = sdr.year*1000+(sdr.day_of_year)+(sdr.msecs_of_day)/(1000.0*3600.0*24.0);
clock_start = sdr.day_of_year+(sdr.msecs_of_day)/(1000.0*3600.0*24.0);
/* use PRI times the number of lines */ 

fprintf(stdout, "num_lines	 	= %d\n",nlines-1);
fprintf(stdout, "good_bytes_per_line 	= %d\n",good_bytes_per_line);
fprintf(stdout, "bytes_per_line 		= %d\n",linelength);
fprintf(stdout, "first_sample 		= %d\n",first_sample);
fprintf(stdout, "num_patches		= %d\n",num_patches);

if (linelength != LINELENGTH) {
	fprintf(stderr, "**** header_linelength (%d) != default linelength (LINELENGTH)\n",linelength);
	}
if (prefix != PREFIX) {
	fprintf(stderr, "**** header_prefix (%d) != default prefix (PREFIX)\n",prefix);
	}
if (suffix != SUFFIX) {
	fprintf(stderr, "**** header_suffix (%d) != default suffix (SUFFIX)\n",suffix);
	}

fseek(indata,LINELENGTH,0);

prior_pri_dn = 0;
data = (char *) malloc(LINELENGTH);
ncnt=0;
print_start = 0;
icu_time_old = 0.;
while(fread(data,LINELENGTH,1,indata)!=0){
        ncnt = ncnt+1;
	memcpy((char *) &info,(data+IFC_OFF),sizeof(struct lineparam));
        icu_time1=(unsigned short *)(data+IFC_OFF-6);
        icu_time2=(unsigned short *)(data+IFC_OFF-4);

/* swap bytes if necessary */

        if(endian == -1) {
           FIX_INT(info.ifc);
           FIX_SHORT(info.swst_dn);
           FIX_SHORT(info.pri_dn);
           FIX_SHORT(*icu_time1);
           FIX_SHORT(*icu_time2);
        }

	info.pri= calc_pri(info);
	info.swst = calc_swst(info);

        icu_time=(double)(*icu_time1)*65536.0+(double)(*icu_time2);
        if(((icu_time-icu_time_old) == 1) & (print_start == 0)) {
                fprintf(stdout,"icu_start               = %.3lf\n",(icu_time-(ncnt-2)*(info.pri*256)));
                print_start = 1;
        }
        icu_time_old = icu_time;

	if (info.pri_dn != prior_pri_dn && iwrite == 0) {

/* don't need to print near_range here because it is obtained from ers_line_fixer */
		/*fprintf(stdout, "near_range		= %lf \n",(info.swst*SOL/2.0)); */
		SC_clock_stop = SC_clock_start+(nlines*info.pri)/(24.0*3600.0);
		clock_stop = clock_start+(nlines*info.pri)/(24.0*3600.0);
		fprintf(stdout, "SC_clock_start		= %16.10lf\n",SC_clock_start);
		fprintf(stdout, "SC_clock_stop		= %16.10lf\n",SC_clock_stop);
		fprintf(stdout, "clock_start		= %16.12lf\n",clock_start);
		fprintf(stdout, "clock_stop		= %16.12lf\n",clock_stop);
		fprintf(stdout, "PRF 			= %lf\n",(1.0/info.pri));
		prior_pri_dn = info.pri_dn;
		iwrite = 1;
		}

	}
}
/* these are taken verbatim from fix_line */
double calc_pri(info)
struct lineparam info;
{
	return(((double)info.pri_dn + 2.0)* SEC_PER_PRI_COUNT);
}
double calc_swst(info)
struct lineparam info;
{
	return((double) info.swst_dn*SEC_PER_PRI_COUNT+9.0*info.pri-6.6E-6);
/* based on Johan Jacob Mohr and Søren Nørvang Madsen, Member, IEEE,Geometric Calibration of ERS
   Satellite SAR Images, 2001, calibriation should be 6.6*/
}
/*_______________________________*/
/* check endian of machine      */
/* 1 if big; -1 if little       */
/*_______________________________*/
int is_big_endian_()
{
        union
        {
        long l;
        char c[sizeof (long) ];
        } u;
        u.l = 1;
        return( u.c[sizeof(long) - 1] ==  1 ? 1 : -1);
}
int is_big_endian__()
{
        return is_big_endian_();
}
