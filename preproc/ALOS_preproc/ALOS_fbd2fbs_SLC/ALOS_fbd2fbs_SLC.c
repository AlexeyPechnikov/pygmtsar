/***************************************************************************
 * ALOS_fbd2fbs_SLC reads a FBD-HH SLC file processed with				   *
 * ALOS_pre_proc_SLC and upsamples it to FBS-HH bandwidth. The algorithm   *
 * takes the fft each complex range line and zero-pads in the frequency    *
 * domain. This is possible because the FBS bandwidth is exactly two times * 
 * the FBD. Note that the interpolated data may exceed the original data   * 
 * span of 0-31 so the numbers are rescaled to lie between 0 and 127 which *
 * still only one byte of storage.                                         *
 * The code uses cfft1d which seems to be the standard interface in the    *
 * InSAR community.                                                        *
 ***************************************************************************/

/***************************************************************************
 * Creator:  David T. Sandwell and Rob Mellors                             *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  08/04/2007                                                    *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * DATE  10/08/2013                                                        *
 * Name: Xiaohua Xu														   *
 * Modification: Modified form ALOS_fbd2fbs, changing the input file to    *
 * SLC. clip127 is switched to clip signed short int.                      *
 *                                                                         *
 ***************************************************************************/

#include "image_sio.h"
#include "siocomplex.h"
#include "lib_functions.h"

#define clip(A) ( ((A) > 32767) ? 32767 : (((A) < -32768) ? -32768 : A) )

char    *USAGE = "ALOS_FBD2FBS_SLC FBD.PRM FBS.PRM \n"
"  FBD.PRM   PRM file for input  image in fine beam dual   polarization (FBD 14 MHz) (input) \n"
"  FBS.PRM   PRM file for outout image in fine beam single polarization (FBS 28 MHz) (output) \n";

int main (int argc, char **argv)
{
FILE	*prmfile, *datafile, *prmout, *dataout;
short int	*indata, *outdata, *extr;
fcomplex *cin, *cout;
float rtest, itest;
int	j, k, nffti, nffto, N=128;
int 	ibufsize, obufsize, fbdsamp, fbssamp;
struct 	PRM r;

	if (argc < 3) die (USAGE,"");
	
	/* flags defined in global_flags.h */
	verbose = 0;	debug = 0;

	/* fill the struct with default parameters */
	null_sio_struct(&r);

	/* open input PRM file and read the parameters */
	if ((prmfile = fopen(argv[1],"r")) == NULL) die("Can't open ",argv[1]);
	get_sio_struct(prmfile, &r);
	
	/* open input raw data file */
	if ((datafile = fopen(r.input_file,"rb")) == NULL) die("Can't open ",r.input_file);
	/* open output PRM file	*/
	if ((prmout  = fopen(argv[2],"w")) == NULL) die("Can't open ",argv[2]);

	/* assemble the output filename and open for writing */
	sscanf(argv[2],"%s",r.input_file);
	r.input_file[strlen(argv[2])-4]=0;
	strcat(r.input_file,".SLC");
	
	/* open output file for single look complex	 image */
	if ((dataout = fopen(r.input_file,"wb")) == NULL) die("Can't open ",r.input_file);

	/* compute the sizes for the input and output buffers and allocate the memory */
	ibufsize = r.num_rng_bins*2;
	if((indata = (short int *) malloc(ibufsize*sizeof(short int))) == NULL){
		fprintf(stderr,"Sorry, couldn't allocate memory for input indata.\n");
		exit(-1);
        }

	fbdsamp = r.num_rng_bins;
	fbssamp = 2 * fbdsamp; 
	obufsize = 2 * fbssamp;

	if((outdata = (short int *) malloc(obufsize*sizeof(short int))) == NULL){
		fprintf(stderr,"Sorry, couldn't allocate memory for output outdata.\n");
		exit(-1);
        }
	
	if((extr = (short int *) malloc(N * sizeof(short int))) == NULL){
		fprintf(stderr,"Sorry, couldn't allocate memory for supplimental outdata.\n");
		exit(-1);
	}
	for (k=0; k<N; k++){
		extr[k] = (short int)(0);
	}
	
	/* find best length of fft (use power of two) for both input and output	*/
	nffti = find_fft_length(fbdsamp);
	nffto = find_fft_length(fbssamp);
	if (debug) fprintf(stderr," nffti %d nffto %d \n",nffti,nffto);

	/* allocate the memory for the complex arrays */
	if((cin = (fcomplex *) malloc(nffti*sizeof(fcomplex))) == NULL){
		fprintf(stderr,"Sorry, couldn't allocate memory for fbd \n");
		exit(-1);
	}
	if((cout = (fcomplex *) malloc(nffto*sizeof(fcomplex))) == NULL){
		fprintf(stderr,"Sorry, couldn't allocate memory for fbs \n");
		exit(-1);
	}
		
	/* read the input and output SLC files */
	for (k=0; k< r.num_lines; k++) {
		fread((void *)indata,sizeof(short int),ibufsize,datafile);

	/* fill the complex array with complex indata */
		for (j=0; j< nffti; j++) {
			if((j < fbdsamp) && (((short int) indata[2*j]) != NULL_DATA) && (((short int) indata[2*j+1]) != NULL_DATA)){
				cin[j].r = (float)(indata[2*j]);
				cin[j].i = (float)(indata[2*j+1]);
            }
            else{
				cin[j].r = 0.0;
				cin[j].i = 0.0;
            }
	  }
	/* interpolate from fbd to fbs */
		rng_expand(cin,nffti,cout,nffto);
	/* convert the complex back to bytes  */

		for (j=0; j< obufsize/2; j++) {

            rtest = rintf(cout[j].r);
			itest = rintf(cout[j].i);			

	/* sometimes the range can exceed -32768~32767 so  
	   clip the numbers to be in the correct range */
			outdata[2*j] = (short int) clip(rtest);
			outdata[2*j+1] = (short int) clip(itest);
		}
		fwrite((void *)outdata,sizeof(short int),obufsize,dataout);
		fwrite((void *)extr,sizeof(short int),N,dataout);
	}

	
	/* compute the changes to the output PRM-file */
	r.xmi = 0;
	r.xmq = 0;
	r.chirp_ext = r.chirp_ext * 2;
	r.good_bytes = 4*fbssamp;
	r.bytes_per_line = 4*fbssamp+N;
	r.num_rng_bins = fbssamp+N/2;
	r.fs = r.fs * 2.;
        if (debug) fprintf(stderr," %d %d %d %d %f %f \n",r.chirp_ext,r.good_bytes,r.bytes_per_line,r.num_rng_bins,r.fs,r.chirp_slope);
	/*  write the output PRM file */
	if ((prmout = fopen(argv[2],"w")) == NULL) die("can't open prfile",argv[2]);
	put_sio_struct(r,prmout);

	free(indata);
	free(outdata);
	free(cin);
	free (cout);
	fclose(prmfile);
  	fclose(prmout);
	fclose(datafile);
	fclose(dataout);
}
