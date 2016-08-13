/***************************************************************************
 * ALOS_filter_SLC applies a low-pass filter to an SLC image     	   *
 ***************************************************************************/

/***************************************************************************
 * Creator:  David T. Sandwell and Rob Mellors                             *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  10/09/2014                                                    *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 ***************************************************************************/

#include"image_sio.h"
#include"siocomplex.h"

#define clip(A) ( ((A) > 32767) ? 32767 : (((A) < -32768) ? -32768 : A) )

char    *USAGE = "ALOS_filter_SLC file.PRM  \n"
                 "  file.PRM  PRM file for input  image \n";

main (int argc, char **argv)
{
FILE	*prmfile, *datafile, *dataout;
short int *indata, *outdata;
fcomplex *cin, *cout;
float rtest, itest;
int	i, j, k, np, nffti, i0;
int 	ibufsize, obufsize, nsamp, headsize;
size_t	n;
struct 	PRM r;

	if (argc < 2) die (USAGE,"");
	
	/* flags defined in global_flags.h */
	verbose = 1;	debug = 1;

	/* fill the struct with default parameters */
	null_sio_struct(&r);

	/* open input PRM file and read the parameters */
	if ((prmfile = fopen(argv[1],"r")) == NULL) die("Can't open ",argv[1]);
	get_sio_struct(prmfile, &r);
	
	/* open input raw data file */
	if ((datafile = fopen(r.input_file,"rb")) == NULL) die("Can't open ",r.input_file);

	/* assemble the output filename and open for writing */
	r.input_file[strlen(argv[1])-4]=0;
	strcat(r.input_file,".SLCF");
	
	/* open output file for single look complex	 image */
	if ((dataout = fopen(r.input_file,"wb")) == NULL) die("Can't open ",r.input_file);

	/* compute the sizes for the input and output buffers and allocate the memory */
	ibufsize = r.num_rng_bins*2;
	if((indata = (short int *) malloc(ibufsize*sizeof(short int))) == NULL){
		fprintf(stderr,"Sorry, couldn't allocate memory for input indata.\n");
		exit(-1);
        }

	nsamp = r.num_rng_bins;
	obufsize = r.num_rng_bins*2;

	if((outdata = (short int *) malloc(obufsize*sizeof(short int))) == NULL){
		fprintf(stderr,"Sorry, couldn't allocate memory for output outdata.\n");
		exit(-1);
        }
	
	/* find best length of fft (use power of two) for both input and output	*/
	nffti = find_fft_length(nsamp);
	if (debug) fprintf(stderr," nffti %d \n",nffti);

	/* allocate the memory for the complex arrays */
	if((cin = (fcomplex *) malloc(nffti*sizeof(fcomplex))) == NULL){
		fprintf(stderr,"Sorry, couldn't allocate memory for input \n");
		exit(-1);
	}
	if((cout = (fcomplex *) malloc(nffti*sizeof(fcomplex))) == NULL){
		fprintf(stderr,"Sorry, couldn't allocate memory for output \n");
		exit(-1);
	}
		
	/* read the input and output SLC files */
	for (k=0; k< r.num_lines; k++) {
		fread((void *)indata,sizeof(short int),ibufsize,datafile);

	/* fill the complex array with complex indata */
		for (j=0; j< nffti; j++) {
			if((j < nsamp) && (((short int) indata[2*j]) != NULL_DATA) && (((short int) indata[2*j+1]) != NULL_DATA)){
				cin[j].r = (float)(indata[2*j]);
				cin[j].i = (float)(indata[2*j+1]);
            		}
            		else{
				cin[j].r = 0.0;
				cin[j].i = 0.0;
            		}
	  	}
	/* filter the data */
	rng_filter(cin,nffti,cout);

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
	}

	free(indata);
	free(outdata);
	free(cin);
	free (cout);
	fclose(prmfile);
	fclose(datafile);
	fclose(dataout);
}
