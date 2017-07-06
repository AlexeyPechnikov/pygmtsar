/********************************************************************************
 * Program to prepare ALOS L1.1 data for InSAR processing.  Read in data and    *
 * leader file and write out prm and integer SLC  file for siosar               *
********************************************************************************
 * Creator:  Rob Mellors and David T. Sandwell                                  *
 *           (San Diego State University, Scripps Institution of Oceanography)  *
 * Date   :  10/03/2007                                                         *
 ********************************************************************************/
/********************************************************************************
 * Modification history:                                                        *
 * Date:                                                                        *
 * 4/23/07 added code to check endianess RJM                                    *
 * 07/15/08   added a command-line option to force earth radius DTS             *
 * 10/10/10   added command line options for SLC_factor and rbias               *
 * *****************************************************************************/

#include"image_sio.h"
#include"lib_functions.h"

char    *USAGE = "\n\nUsage: ALOS_pre_process_SLC imagefile LEDfile [-radius RE] [-swap] [-V] [-debug] [-quiet] \n"
"\ncreates data.SLC and writes out parameters (PRM format) to stdout\n"
"\nimagefile 	ALOS Level 1.1 complex file (CEOS format):\n"
"LEDfile 	ALOS Level 1.1 LED file (CEOS leaderfile format):\n"
"\n options: \n"
"-near near_range       specify the near_range (m) \n"
"-radius RE             specify the local earth radius (m) \n"
"-SLC_factor fact       scale factor to convert float to int SLC [1.0] \n"
"-swap                  do byte-swap (should be automatic) \n"
"-ALOS1                 ALOS2 L1.1 data format \n"
"-ALOS2                 ALOS2 L1.1 data format (default)\n"
"-LED                   write generic LED file\n"
"-noLED                 oldstyle use ldr file for orbits directlyn"
"-V                     verbose write information \n"
"-debug                 write even more information \n"
"-quiet                 don't write any information \n"
"-tbias tbias           correct the clock bias (positive value means plus)\n"
"-rbias tbias           correct the range bias (positive value means increase near range)\n"
"Example:\n"
"ALOS_pre_process_SLC IMG-HH-ALOS2011986990-140813-HBQR1.1__A LED-ALOS2011986990-140813-HBQR1.1__A -SLC_factor 1. -rbias -70.0000 -tbias 0.068759 \n";

long read_ALOS_data_SLC (FILE *, FILE *, struct PRM *, long *);
void parse_ALOS_commands(int, char **, char *, struct PRM *);
void set_ALOS_defaults(struct PRM *);
void print_ALOS_defaults(struct PRM *);
void swap_ALOS_data_info(struct sardata_info *);
void get_files(struct PRM *, FILE **, FILE **, char *, char *, int);

int ledflag;
int main (int argc, char **argv) 
{
FILE	*imagefile, *ldrfile;
FILE	*rawfile[11], *prmfile[11];
char	prmfilename[128];
int	nPRF;
long	byte_offset;
struct 	PRM prm;
struct 	ALOS_ORB orb;
/*char   	date[8];*/

	if (argc < 3) die (USAGE,"");

	/* set flags  */
	roi = debug = verbose = swap = quiet_flag = 0;
        slc_fact = 0.01;
        tbias = 0.0;
        rbias = 0.0;
        prefix_off = 132;

        /* default is to use the new LED orbit */
        ledflag = 1;

	nPRF = 0;

	null_sio_struct(&prm);
	set_ALOS_defaults(&prm);

	/* read command line */
	parse_ALOS_commands(argc, argv, USAGE, &prm);

        /* shift the start time if this is ALOS1 */
        if(prefix_off == 0) tbias = tbias - 0.0020835;

	if (verbose) print_ALOS_defaults(&prm);
	if (is_big_endian_() == -1) {swap = 1;fprintf(stderr,".... swapping bytes\n");} else {swap = 0;} 

	/* IMG and LED files should exist already */
	if ((imagefile = fopen(argv[1], "r")) == NULL) die ("couldn't open Level 1.1 IMG file \n",argv[1]);
	if ((ldrfile = fopen(argv[2], "r")) == NULL) die ("couldn't open LED file \n",argv[2]); 

	/* if it exists, copy to prm structure */
	strcpy(prm.led_file,argv[2]);

	/* name and open output files and header files for raw data (but input for later processing) */
	get_files(&prm, &rawfile[nPRF], &prmfile[nPRF], prmfilename, argv[1], nPRF);

	/* read sarleader; put info into prm; write log file if specified 		*/
	read_ALOS_sarleader(ldrfile, &prm, &orb);

// AUGUST 2016
       /* write out orbit params in generic LED format */
       if (ledflag) write_ALOS_LED(&orb, &prm, argv[1]);
// AUGUST 2016

	/* read Level 1.1 file;  put info into prm; convert to *.SLC format 		*/
	/* if PRF changes halfway through, create new set of header and data files      */
	/* byte_offset is non-zero only if the prf changes				*/
	/* byte_offset gets set to point in file at prf change				*/

	byte_offset = -1;
	while (byte_offset != 0){

		/* if prf changes, create new prm and data files			*/
		if (nPRF > 0 ) {
			if (verbose) fprintf(stderr,"creating multiple files due to PRF change (*.%d) \n",nPRF+1);
			get_files(&prm, &rawfile[nPRF], &prmfile[nPRF], prmfilename, argv[1], nPRF);
			}

		/* set the chirp extension to 500 if FBD fs = 16000000 */
        	if (prm.fs < 17000000.) {
			prm.chirp_ext = 500;
			prm.chirp_slope =  -5.18519e+11;
 	                prm.SLC_scale = I2SCALE;  /* set dfact */

		} else {
			prm.chirp_slope = -1.03704e+12;
 	                prm.SLC_scale = I2SCALE * 2;
		}

		/* read_ALOS_data returns 0 if all data file is read;
		returns byte offset if the PRF changes  */
		/* calculate parameters from orbit */
		byte_offset = read_ALOS_data_SLC(imagefile, rawfile[nPRF], &prm, &byte_offset);

		ALOS_ldr_orbit(&orb, &prm);

		/* force chirp slope if asked to */
		if (force_slope == 1) prm.chirp_slope = forced_slope;

		/* set parameters for integer SLC */

		/* write ascii output, SIO format */
                prm.near_range = prm.near_range + rbias;

                /* if this is ALOS-2 data after 2014 set the Dopplers to zero */
                if(prm.SC_clock_start/1000. > 2014.) {
                   prm.fd1 = 0.;
                   prm.fdd1 = 0.;
                   prm.fddd1 = 0.;
                }
		put_sio_struct(prm, prmfile[nPRF]);

		/* write roi_pac output 
		if (roi) {
			// first part of rsc file
			write_roi(argv[1], ldrfile, prm, orb, date);
			// orbit file 
			write_roi_orbit(orb, date);
			}*/

		nPRF++;
		}

	return(EXIT_SUCCESS);
}
/*------------------------------------------------------*/
void get_files(struct PRM *prm, FILE **rawfile, FILE **prmfile, char *prmfilename, char *name, int n)
{
	/* name and open output file for raw data (but input for later processing)      */
	/* if more than 1 set of output files, append an integer (beginning with 2)     */

	if (n == 0) {
		sprintf(prm->input_file,"%s.SLC", name);
		sprintf(prmfilename,"%s.PRM", name);
	} else {
		sprintf(prm->input_file,"%s.SLC.%d",name,n+1);
		sprintf(prmfilename,"%s.PRM.%d", name, n+1);
	}
	sprintf(prm->SLC_file,"%s.SLC", name);
	strcpy(prm->dtype, "a");

	/* now open the files */
	if ((*rawfile = fopen(prm->input_file,"w")) == NULL) die("can't open ",prm->input_file);

	if ((*prmfile = fopen(prmfilename, "w")) == NULL) die ("couldn't open output PRM file \n",prmfilename);

}
/*------------------------------------------------------*/
