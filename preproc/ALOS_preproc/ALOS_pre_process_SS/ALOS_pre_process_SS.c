/********************************************************************************
 * Program to prepare ALOS L1.0 data for InSAR processing.  Read in data and    *
 * leader file and write out prm and raw file for siosar                        *
********************************************************************************
 * Creator:  Rob Mellors and David T. Sandwell                                  *
 *           (San Diego State University, Scripps Institution of Oceanography)  *
 * Date   :  10/03/2007                                                         *
 ********************************************************************************/
/********************************************************************************
 * Modification history:                                                        *
 * Date:                                                                        *
 * 4/23/07   added code to check endianess RJM                                  *
 * 07/15/08  added a command-line option to force earth radius DTS              *
 * 10/02/08  added capability to pre_process ScanSAR WB1 data DTS               *
 * *****************************************************************************/

#include"image_sio.h"
#include"lib_functions.h"

char    *USAGE = "\n\nUsage: ALOS_pre_process_SS imagefile LEDfile [-near near_range] [-radius RE] [-swath swath#] [-burst_skip] [-num_burst] [-swap] [-V] [-debug] [-quiet] \n"
"\ncreates data.raw and writes out parameters (PRM format) to stdout\n"
"\nimagefile 	ALOS Level 1.0 complex file (CEOS format):\n"
"LEDfile 	ALOS Level 1.0 LED file (CEOS leaderfile format):\n"
"\n options: \n"
"-near near_range       specify the near_range (m) \n"
"-radius RE             specify the local earth radius (m) \n"
"-swath                 specify swath number 1-5 [default 4]  \n"
"-burst_skip            number of burst to skip before starting output (1559 lines/burst) \n"
"-num_burst             number of burst to process [default all] \n"
"                       there are 72 bursts in a WB1 frame \n"
"-swap                  do byte-swap (should be automatic) \n"
"-fd1 [DOPP]            sets doppler centroid [fd1] to DOPP\n"
"-V                     verbose write information) \n"
"-debug                 write even more information \n"
"-quiet                 don't write any information \n"
"Example:\n"
"ALOS_pre_process_SS IMG-HH-ALPSRS049842950-W1.0__D LED-ALPSRS049842950-W1.0__D -near 847916 -radius 6371668.872945 -burst_skip 5 -num_burst 36 \n\n"
" burst #   look_angle  #lines_burst \n"
"   1         20.1        247     \n"
"   2         26.1        356     \n"
"   3         30.6        274     \n"
"   4         34.1        355     \n"
"   5         36.5        327   \n\n";

long read_ALOS_data_SS (FILE *, FILE *, struct PRM *, long *, int *, int *, int *);
void parse_ALOS_commands(int, char **, char *, struct PRM *, int *, int *, int *);
void set_ALOS_defaults(struct PRM *);
void print_ALOS_defaults(struct PRM *);
void swap_ALOS_data_info(struct sardata_info *);
void get_files(struct PRM *, FILE **, FILE **, char *, char *, int, int);

int main (int argc, char **argv) 
{
FILE	*imagefile, *ldrfile;
FILE	*rawfile[11], *prmfile[11];
char	prmfilename[128];
int	nPRF;
long	byte_offset;
int     nsub = 3, burst_skip = 0, num_burst = 1000 ;
struct PRM prm;
struct ALOS_ORB orb;

	if (argc < 3) die (USAGE,"");

	/* set flags  */
	dopp = 0;
	debug = verbose = swap = quiet_flag = 0;

	nPRF = 0;

	null_sio_struct(&prm);
	set_ALOS_defaults(&prm);

	/* read command line */
	parse_ALOS_commands(argc, argv, USAGE, &prm, &nsub, &burst_skip, &num_burst);

	if (verbose) print_ALOS_defaults(&prm);
	if (is_big_endian_() == -1) {swap = 1;fprintf(stderr,".... swapping bytes\n");} else {swap = 0;} 

	/* IMG and LED files should exist already */
	if ((imagefile = fopen(argv[1], "r")) == NULL) die ("couldn't open Level 1.0 IMG file \n",argv[1]);
	if ((ldrfile = fopen(argv[2], "r")) == NULL) die ("couldn't open LED file \n",argv[2]); 

	/* if it exists, copy to prm structure */
	strcpy(prm.led_file,argv[2]);

	/* name and open output files and header files for raw data (but input for later processing) */
	get_files(&prm, &rawfile[nPRF], &prmfile[nPRF], prmfilename, argv[1], nPRF, nsub+1);

	/* read sarleader; put info into prm; write log file if specified 		*/
	read_ALOS_sarleader(ldrfile, &prm, &orb);

	/* read Level 1.0 file;  put info into prm; convert to *.raw format 		*/
	/* if PRF changes halfway through, create new set of header and data files      */
	/* byte_offset is non-zero only if the prf changes				*/
	/* byte_offset gets set to point in file at prf change				*/

	byte_offset = -1;

	/* set the chirp extension to 500 if FBD  or WB1 fs = 16000000 */
       	if (prm.fs < 17000000.) prm.chirp_ext = 500;

	/* check to be sure the burst_skip*1559 lines does not exceed the file size */
        if (burst_skip*1559 > 100000) fprintf(stderr," warning skip could exceed the length of the file \n");

	/* read_ALOS_data returns 0 if all data file is read;
	returns byte offset if the PRF changes  */
	byte_offset = read_ALOS_data_SS(imagefile, rawfile[nPRF], &prm, &byte_offset, &nsub, &burst_skip, &num_burst);

	/* calculate parameters from orbit */
	ALOS_ldr_orbit(&orb, &prm);

	/* calculate doppler from raw file */
	if (dopp == 1) calc_dop(&prm);

	/* write ascii output, SIO format */
	put_sio_struct(prm, prmfile[nPRF]);

	return(EXIT_SUCCESS);
}
/*------------------------------------------------------*/
void get_files(struct PRM *prm, FILE **rawfile, FILE **prmfile, char *prmfilename, char *name, int n, int nsub)
{
	/* name and open output file for raw data (but input for later processing)      */
	/* if more than 1 set of output files, append an integer (beginning with 2)     */

	if (n == 0) {
		sprintf(prm->input_file,"%s_SW%d.raw", name, nsub);
		sprintf(prmfilename,"%s_SW%d.PRM", name, nsub);
	} else {
		sprintf(prm->input_file,"%s_SW%d.raw.%d",name, nsub,n+1);
		sprintf(prmfilename,"%s_SW%d.PRM.%d", name, nsub, n+1);
	}

	/* now open the files */
	if ((*rawfile = fopen(prm->input_file,"w")) == NULL) die("can't open ",prm->input_file);

	if ((*prmfile = fopen(prmfilename, "w")) == NULL) die ("couldn't open output PRM file \n",prmfilename);

}
/*------------------------------------------------------*/
