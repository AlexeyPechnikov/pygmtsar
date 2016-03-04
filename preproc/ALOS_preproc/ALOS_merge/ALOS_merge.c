/********************************************************************************
 * Program to merge adjacent frames of ALOS palsar data                         * 
 *******************************************************************************/
/********************************************************************************
 * Creator:  Rob Mellors                                                        *
 *           (San Diego State University)                                       *
 * Date   :  10/03/2007                                                         *
 ********************************************************************************/
/********************************************************************************
 * Modification history:                                                        *
 * Date:                                                                        *
 *******************************************************************************/

#include "image_sio.h"
#include "lib_functions.h"

char    *USAGE = "Usage: \n"
"ALOS_merge PRMfile1 PRMfile2 [-output outfile] [-V] [-nodopp]\n\n"
" appends file2 onto file1 (omitting duplicate lines)\n"
" updates num_lines and doppler (average) in PRM file\n"
" output named merge.raw by default\n"
" if outfile specified, creates outfile.raw and outfile.PRM\n"
" -V verbose\n"
" -nodopp do not calculate doppler (default is calculate doppler)\n"
"example: ALOS_merge IMG-HH-ALPSRP057130830-H1.0__A.PRM IMG-HH-ALPSRP057130840-H1.0__A.PRM -outout foo -V\n";

void parse_ALOS_merge(int, int, char **, char *, char *);
void assign_output(char *,char *, char *);
void fill_PRM(struct PRM *, FILE *);
long read_ALOS_raw_data(FILE *, FILE *, struct PRM *, int *);
int 
main (int argc, char **argv)
{
char	RawOutName[128], tmp[128];
char	PRMOutName[128];
int	start_frame, narg;
long	nlines1, nlines2, nlines;
double	overlap;
FILE	*PRMfile1, *RAWfile1;
FILE	*PRMfile2, *RAWfile2;
FILE	*PRMfile3, *RAWfile3;

struct PRM p1;			
struct PRM p2;			
struct PRM p3;			

	narg = 3;
	if (argc < narg) die(USAGE,"");

	verbose = 0;

	/* default name */
	strcpy(tmp, "merge");
	if (argc > narg) parse_ALOS_merge(argc, narg, argv, tmp, USAGE);
	assign_output(RawOutName, PRMOutName, tmp);

	/* open files */
	if ((PRMfile1 = fopen(argv[1],"r")) == NULL) die("can't open prfile",argv[1]);
	if ((PRMfile2 = fopen(argv[2],"r")) == NULL) die("can't open prfile",argv[2]);
	if ((PRMfile3 = fopen(PRMOutName,"w")) == NULL) die("can't open prfile",PRMOutName);

	/* read prm values; put file1 prm into file 3 */
	fill_PRM(&p1, PRMfile1);
	fill_PRM(&p2, PRMfile2);
	rewind(PRMfile1);
	fill_PRM(&p3, PRMfile1);

	strcpy(p3.input_file,RawOutName);

	/* open raw files */
	if ((RAWfile1 = fopen(p1.input_file,"r")) == NULL) die("cannot open ",p1.input_file);
	if ((RAWfile2 = fopen(p2.input_file,"r")) == NULL) die("cannot open ",p2.input_file);

	/* output */
	if ((RAWfile3 = fopen(RawOutName,"w")) == NULL) die("cannot open output file ",RawOutName);

	/* check whether scenes overlap in time */
	overlap = p1.SC_clock_stop - p2.SC_clock_start;
	if (overlap < 0) die("frames do not overlap in time","");
	if (verbose) fprintf(stderr," frame 1 %lf %lf\n frame 2 %lf %lf overlap %lf\n"
		,p1.SC_clock_start, p1.SC_clock_stop
		,p2.SC_clock_start, p2.SC_clock_stop, overlap);

	/* read files and write */
	/* write out all lines with frame numner > last_frame */
	start_frame = -1;
	if (verbose) fprintf(stderr," reading file 1 %s\n",p1.input_file);
	nlines1 = read_ALOS_raw_data(RAWfile1, RAWfile3, &p1, &start_frame);

	if (verbose) fprintf(stderr," reading file 2 %s\n",p2.input_file);
	nlines2 = read_ALOS_raw_data(RAWfile2, RAWfile3, &p2, &start_frame);
	fclose(RAWfile3);
	/* end of data writing */

	nlines = nlines1 + nlines2;
	if (verbose) fprintf(stderr," nlines1 %ld nlines2 %ld nlines %ld \n",nlines1,nlines2,nlines);

	/* write out all non-null parameters */
	p3.SC_clock_stop = p2.SC_clock_stop;
	p3.ht_end = p2.ht_end;
	p3.num_lines = nlines1 + nlines2;
	p3.num_patches = (int)((1.0*nlines)/(1.0*p3.num_valid_az));

	/* now recalculate doppler for the whole file */
	if (dopp == 1) calc_dop(&p3); 

	/* write out PRM file */
	put_sio_struct(p3, PRMfile3);

	return(EXIT_SUCCESS);
}
/*-----------------------------------------------*/
void fill_PRM(struct PRM *p, FILE *Pfile)
{
	/* set all prm parameters in structure to NULL values */
	null_sio_struct(p);

	/* now rad new ones 			*/
	get_sio_struct(Pfile, p);
}
/*-----------------------------------------------*/
void assign_output(char *rawname, char *prmname, char *base)
{
	strcpy(rawname, base);
	strcat(rawname,".raw");

	strcpy(prmname, base);
	strcat(prmname,".PRM");
}
/*-----------------------------------------------*/
/* reads options 				*/
/* start with third argument			*/

void parse_ALOS_merge(int na, int nstart, char **a, char *filename, char *USAGE)
{
int	n;

verbose = 0;
debug = 0;
dopp = 1;		/* default is to calculate doppler */

for (n = nstart; n < na; n++) {
	if (!strcmp(a[n], "-output")) {
		n++;
		if (n > na) die (" no option after -output!\n",""); 
		strcpy(filename, a[n]);
		fprintf(stderr," setting output file %s \n",filename);
		n++;
	} else if (!strcmp(a[n], "-V")) {
		verbose = 1;
		fprintf(stderr," verbose output \n");
		n++;
	} else if (!strcmp(a[n], "-v")) {
		verbose = 1;
		fprintf(stderr," verbose output \n");
		n++;
	} else if (!strcmp(a[n], "-nodopp")) {
		dopp = 0;
		fprintf(stderr," not calculating doppler centroid \n");
		n++;
	} else if (!strcmp(a[n], "-debug")) {
		debug = 1;
		fprintf(stderr," debugging output \n");
		n++;
	} else {
		fprintf(stderr," %s *** option not recognized ***\n\n",a[n]);
		fprintf(stderr,"%s",USAGE);
		exit(1);
		}
	}
}
/*-----------------------------------------------*/
/* read raw files								*/

long read_ALOS_raw_data (FILE *imagefile, FILE *outfile, struct PRM *prm, int *start_frame) 
{

	struct	sardata_info sdr;

	char	*data;
	int	line_prefix_size, line_data_size;
	int	first_frame, frame=0;			/* frame counter */
	long	n, nlines;

	line_prefix_size = sizeof(struct sardata_info);
	line_data_size = prm->bytes_per_line - line_prefix_size;

	data = malloc(line_data_size * sizeof(char));

	if (debug) fprintf(stderr,".... reading data (%d bytes)\n", line_prefix_size+line_data_size);

	n = 0;
	nlines = 0;
	while ( (fread((void *) &sdr, line_prefix_size, 1, imagefile)) == 1 ) {

		if (debug) fprintf(stderr, SARDATA__WCS, SARDATA_RVL(sdr));

		frame = sdr.frame_counter;

		if (n == 0) {
			first_frame = frame;
			if (verbose) fprintf(stderr, "first frame %d (start_frame %d)\n", first_frame, *start_frame);
			}


		/* write line if not a duplicate */
		if (frame > *start_frame) {
			if (verbose) if (nlines == 0) fprintf(stderr," first line written %ld frame %d \n", n, frame);

			/* read the data */
			fread(data, sizeof(char), line_data_size, imagefile);

			/* write header */
			fwrite(&sdr, sizeof(char), line_prefix_size, outfile);
			/* write data */
			fwrite(data, sizeof(char), line_data_size, outfile);

			nlines++;
			} else {
			fseek(imagefile, line_data_size, SEEK_CUR);
			}
        	n++;
 		}
      
	/* this is the first frame of the next file */
	*start_frame = frame;

	if (verbose) fprintf(stderr, "last frame %d (nlines written %ld)\n",*start_frame,nlines);

	free(data);

	return(nlines);
}
/***************************************************************************/
double get_clock(struct sardata_info sdr, double tbias)
{
double time;

	//nsd = 24.0*60.0*60.0;	/* seconds in a day */

	time = ((double) sdr.sensor_acquisition_year)*1000 +
		(double) sdr.sensor_acquisition_DOY +
		(double) sdr.sensor_acquisition_msecs_day/1000.0/86400.0 +
		tbias/86400.0;

	return(time);
}
/***************************************************************************/
