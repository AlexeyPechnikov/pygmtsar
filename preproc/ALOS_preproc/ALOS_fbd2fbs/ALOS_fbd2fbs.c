/***************************************************************************
 * ALOS_fbd_fbs reads a raw FBD-HH file processed with ALOS_pre_proc and   * 
 * upsamples it to FBS-HH bandwidth.  The algorithm takes the fft each     *
 * complex range line and zero-pads in the frequency domain.  This is      *
 * possible because the FBS bandwidth is exactly two times the FBD.        *
 * Note that the interpolated data may exceed the original data span of    *
 * 0-31 so the numbers are rescaled to lie between 0 and 127 which still   *
 * only one byte of storage.                                               *
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
 * DATE                                                                    *
 *                                                                         *
 ***************************************************************************/

#include"image_sio.h"
#include"siocomplex.h"
#include "lib_functions.h"

#define clip127(A) ( ((A) > 127) ? 127 : (((A) < 0) ? 0 : A) )

char    *USAGE = "ALOS_FBD2FBS FBD.PRM FBS.PRM \n"
"  FBD.PRM   PRM file for input  image in fine beam dual   polarization (FBD 14 MHz) (input) \n"
"  FBS.PRM   PRM file for outout image in fine beam single polarization (FBS 28 MHz) (output) \n";

int main (int argc, char **argv)
{
FILE	*prmfile, *datafile, *prmout, *dataout;
unsigned char	*indata, *outdata;
fcomplex *cin, *cout;
float rtest, itest;
int	i, j, k, nffti, nffto;
int 	ibufsize, obufsize, fbdsamp, fbssamp, headsize;
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
	if ((datafile = fopen(r.input_file,"r")) == NULL) die("Can't open ",r.input_file);

	/* open output PRM file	*/
	if ((prmout  = fopen(argv[2],"w")) == NULL) die("Can't open ",argv[2]);

	/* assemble the output filename and open for writing */
        sscanf(argv[2],"%s",r.input_file);
        r.input_file[strlen(argv[2])-4]=0;
        strcat(r.input_file,".raw");
	/* open output file for single look complex	 image */
	if ((dataout = fopen(r.input_file,"w")) == NULL) die("Can't open ",r.input_file);

	/* compute the sizes for the input and output buffers and allocate the memory */
	ibufsize = r.bytes_per_line;
	if((indata = (unsigned char *) malloc(2*ibufsize*sizeof(unsigned char))) == NULL){
          fprintf(stderr,"Sorry, couldn't allocate memory for input indata.\n");
          exit(-1);
        }
        fbdsamp = r.good_bytes/2-r.first_sample;
        fbssamp = 2 * fbdsamp; 
	headsize = 2 * r.first_sample;
        obufsize = 2*(fbssamp+r.first_sample);
	if((outdata = (unsigned char *) malloc(2*obufsize*sizeof(unsigned char))) == NULL){
          fprintf(stderr,"Sorry, couldn't allocate memory for output outdata.\n");
          exit(-1);
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

	/* read and write the input and output raw files */
          for (k=0; k< r.num_lines; k++) {
	  fread((void *)indata,sizeof(unsigned char),ibufsize,datafile);
	  fwrite((void *)indata,sizeof(unsigned char),headsize,dataout); 

	/* fill the complex array with complex indata */
	  for (j=0; j< nffti; j++) {
            i = j + r.first_sample;
	    if((j < fbdsamp) && (((int) indata[2*i]) != NULL_DATA) && (((int) indata[2*i+1]) != NULL_DATA)){
	      cin[j].r = (float)(indata[2*i]-r.xmi);
	      cin[j].i = (float)(indata[2*i+1]-r.xmq);
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
            i = j + r.first_sample;

	/* increase dynamic range by 2 and set the mean value to 63.5 */
            rtest = rintf(2.*cout[j].r+63.5);
	    itest = rintf(2.*cout[j].i+63.5);

	/* sometimes the range can exceed 0-127 so  
	   clip the numbers to be in the correct range */
            outdata[2*i] = (unsigned char) clip127(rtest);
            outdata[2*i+1] = (unsigned char) clip127(itest);
	  }
  	  fwrite((void *)outdata+headsize,sizeof(unsigned char),obufsize-headsize,dataout); 
	}

	/* compute the changes to the output PRM-file */
	r.xmi = 63.5;
	r.xmq = 63.5;
	r.chirp_ext = r.chirp_ext * 2;
        r.good_bytes = 2*(fbssamp+r.first_sample);
	r.bytes_per_line = obufsize;
	r.num_rng_bins = fbssamp+r.chirp_ext;
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
