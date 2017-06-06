/***************************************************************************
 * Creator:  Anders Hogrelius                                              *
 *           Earth Consultants International, Inc                          *
 * Date   :  04/17/2017                                                    *
 *                                                                         *
 * Based on make_slc_s1a by Xiaohua(Eric) XU                               *
 *           (Scripps Institution of Oceanography)                         *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * Date   : 05/17/2010, Anders Hogrelius                                   *
 *          Added range and azimuth bias for ERS-1 and ERS-2               *
 *                                                                         *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <time.h>
#include "tiffio.h"
#include "PRM.h"
#include "lib_functions.h"
#include "stateV.h"
#include "xmlC.h"
#include "lib_defs.h"
#include "epr_api-2.3/src/epr_api.h"
#if defined(WIN32) && defined(_DEBUG)
#include <crtdbg.h>
#endif /* if defined(WIN32) && defined(_DEBUG) */
//#include <typeinfo>

int write_orb(state_vector *sv, FILE *fp, int);
int dump_data(EPR_ELogLevel log_level, const char *infile, FILE* outstream, int pixflag, unsigned int l0, unsigned int lN, unsigned int p0, unsigned int pN);
char *remove_ext (char* mystr, char dot, char sep);
int read_header(EPR_ELogLevel log_level, const char *infile, struct PRM * prm, state_vector *sv, int * n_state_vectors);
int monthtonum(char * szMonth);

char *USAGE = "\n\nUsage: envi_slc_decode name_of_data_file (*.E1 *.E2 or *.N1) \n"
"\nExample: envi_slc_decode SAR_IMS_1PNESA19920914_183429_00000018C087_00442_06098_0000.E1 \n"
"\nOutput: SAR_IMS_1PNESA19920914_183429_00000018C087_00442_06098_0000.SLC SAR_IMS_1PNESA19920914_183429_00000018C087_00442_06098_0000.PRM SAR_IMS_1PNESA19920914_183429_00000018C087_00442_06098_0000.LED\n";

int main(int argc, char **argv){
    
    FILE *OUTPUT_PRM,*OUTPUT_SLC,*OUTPUT_LED;
    char *filename=0;
    char tmp_str[200];
    struct PRM prm;
    state_vector sv[200];
    int n=0;
    
    /* EPR_LogLevel can be set to e_log_debug, e_log_info, e_log_warning or e_log_error */
    enum EPR_LogLevel eloglevel = e_log_error;

    if (argc != 2) die (USAGE,"");
      
    // initiate the prm
    null_sio_struct(&prm);

    if (read_header(eloglevel, argv[1], &prm, sv, &n) != 0) 
	{ 
	die ("Couldn't read envisat header in input file: \n",argv[1]);
	}

    // generate the PRM file, test to make sure we can write to it by opening the file
    strcpy(tmp_str,argv[1]);
    filename=remove_ext (tmp_str, '.', '/');
    strcpy(tmp_str,filename);
    free(filename);
    strcat(tmp_str,".PRM");



    if ((OUTPUT_PRM = fopen(tmp_str,"wb")) == NULL)
	{ 
	die ("Couldn't open prm file: \n",tmp_str);
	}
    else
	{
	put_sio_struct(prm, OUTPUT_PRM);
	fclose(OUTPUT_PRM);
    	}
    
    // generate the LED file, test to make sure we can write to it by opening the file
    strcpy(tmp_str,argv[1]);
    filename=remove_ext (tmp_str, '.', '/');
    strcpy(tmp_str,filename);
    free(filename);
    strcat(tmp_str,".LED");

    if ((OUTPUT_LED = fopen(tmp_str,"wb")) == NULL)
	{ 
	die ("Couldn't open slc file: \n",tmp_str);
	}
    else
	{
	write_orb(sv,OUTPUT_LED,n);
	fclose(OUTPUT_LED);
    	}
  
    // generate the SLC file, test to make sure we can write to it by opening the file and closing it for writing with dump data
    strcpy(tmp_str,argv[1]);
    filename=remove_ext (tmp_str, '.', '/');
    strcpy(tmp_str,filename);
    free(filename);
    strcat(tmp_str,".SLC");

    if ((OUTPUT_SLC = fopen(tmp_str,"wb")) == NULL)
	{ 
	die ("Couldn't open slc file: \n",tmp_str);
	}
    else
	{
	dump_data(eloglevel, argv[1], OUTPUT_SLC, 1, 0, 0, 0, 0);
	fclose(OUTPUT_SLC);
    	}

    
    return 0;

}

int write_orb(state_vector *sv, FILE *fp, int n){
    int i;
    double dt;
    
    dt = trunc((sv[1].sec)*100.0)/100.0-trunc((sv[0].sec)*100.0)/100.0;
    if(n<=1) return(-1);
    fprintf(fp,"%d %d %d %.6lf %lf \n",n,sv[0].yr,sv[0].jd,sv[0].sec,dt);
    for(i=0;i<n;i++){
        fprintf(fp,"%d %d %.6lf %.6lf %.6lf %.6lf %.8lf %.8lf %.8lf \n",sv[i].yr,sv[i].jd,sv[i].sec,sv[i].x,sv[i].y,sv[i].z,sv[i].vx,sv[i].vy,sv[i].vz);
    }
    return(1);
}


int dump_data(EPR_ELogLevel log_level, const char *infile, FILE* outstream, int pixflag, unsigned int l0, unsigned int lN, unsigned int p0, unsigned int pN) 
  {

  EPR_SProductId* product_id;
  EPR_SDatasetId* MAIN_PROC_PM_ID;
  EPR_SDatasetId* MDS1;
  EPR_SRecord*    rec1 = NULL;
  EPR_SRecord*    rec5 = NULL;
  EPR_SField      numlines_field;
  EPR_SField      numpixels_field;
  EPR_SField      line_field;
  ulong           numlines;
  ulong           numberoflines;
  ulong           numpixels;
  int             cnt;
  int             x,y;  /* loop counter go upto 25.000 or so */
  short           realpart_short;/* written to output file */
  short           imagpart_short;/* written to output file */
  int             status;
 

  /* Initialize the API. Set log-level to DEBUG and use default log-output (stdout) */
  epr_init_api(log_level, epr_log_message, NULL);

  /* Open the product; an argument is a path to product data file */
  /* PRODUCT dataset record field element */
  product_id      = epr_open_product(infile);
  /* product DATASET record field element */
  MAIN_PROC_PM_ID = epr_get_dataset_id(product_id, "MAIN_PROCESSING_PARAMS_ADS");
  MDS1            = epr_get_dataset_id(product_id, "MDS1");
  /* product dataset RECORD field element */
  rec1 = epr_read_record(MAIN_PROC_PM_ID,         0, NULL);
  /* product dataset record FIELD element */
  numlines_field  = *(epr_get_field(rec1, "num_output_lines"));
  numpixels_field = *(epr_get_field(rec1, "num_samples_per_line"));
  /*
  epr_free_record(rec1);
  */
  epr_print_field(&numlines_field, stdout);
  epr_print_field(&numpixels_field, stdout);
  numlines        = epr_get_field_elem_as_uint(&numlines_field, 0);
  numpixels       = epr_get_field_elem_as_uint(&numpixels_field, 0);


  if (pixflag)
    {
    l0 = 1;
    lN = numlines;
    p0 = 1;
    pN = numpixels;
    }

  /* loop over data record to get data and dump it to file */
  numberoflines = epr_get_num_records(MDS1);
  if(log_level == e_log_debug)
    {  
  printf("numberoflines: %f\n", (float)numberoflines);
  printf("numlines: %f\n", (float)numlines);
  printf("numpixels: %f\n", (float)numpixels);
  printf("l0: %f\n", (float)l0);
  printf("lN: %f\n", (float)lN);
  printf("p0: %f\n", (float)p0);
  printf("pN: %f\n", (float)pN);
    }

  /* check if number of records indeed equals the number of lines */
  if (numberoflines != numlines)
    {
    printf("numlines not equal in check, ASAR format error?.");
    return 1;
    }


  /* --- Check if input as acceptable ---------------------------- */
  if (l0 < 1)  
    {
    printf("l0<1 not allowed.  first line is 1 not 0.\n");
    return 1;
    }
  if (p0 < 1)  
    {
    printf("p0<1 not allowed.  first line is 1 not 0.\n");
    return 1;
    }
  if (lN > numlines)  
    {
    printf("lN>numlines not allowed.\n");
    return 1;
    }
  if (pN > numpixels)  
    {
    printf("pN>numpixels not allowed.\n");
    return 1;
    }


  /* --- read in whole line of cpx data -------------------------- */
  
  for (y=l0;y<=lN;y++)
    {
    rec5       = epr_read_record(MDS1, y-1, NULL);
    line_field = *(epr_get_field(rec5, "proc_data"));
    cnt        = 2*(p0-1);/* p0 starts at 1 for first element (BUGFIX!) */
    /* write out selected pixels */
    for (x=p0;x<=pN;x++)
      {
      realpart_short = epr_get_field_elem_as_short(&line_field,cnt);
      cnt++;
      imagpart_short = epr_get_field_elem_as_short(&line_field,cnt);
      cnt++;
      status = fwrite(&realpart_short,2,1,outstream);
      if (status != 1 && log_level == e_log_debug) fprintf(stderr,"fwrite could not write to disk?");
      status = fwrite(&imagpart_short,2,1,outstream);
      if (status != 1 && log_level == e_log_debug) fprintf(stderr,"fwrite could not write to disk?");
      }
    /* this program seemed to fill all memory for some reason?  try to free it */
    epr_free_record(rec5);
    }
  // below it's closing the file twice, may cause bug on Linux
  // fclose(outstream);
  epr_close_product(product_id);
  /* Closes product reader API, release all allocated resources */
  epr_close_api();
  return 0;
  }

// remove_ext: removes the "extension" from a file spec.
//   mystr is the string to process.
//   dot is the extension separator.
//   sep is the path separator (0 means to ignore).
// Returns an allocated string identical to the original but
//   with the extension removed. It must be freed when you're
//   finished with it.
// If you pass in NULL or the new string can't be allocated,
//   it returns NULL.

char *remove_ext (char* mystr, char dot, char sep) 
  {
    char *retstr, *lastdot, *lastsep;

    // Error checks and allocate string.

    if (mystr == NULL)
        return NULL;
    if ((retstr = malloc (strlen (mystr) + 1)) == NULL)
        return NULL;

    // Make a copy and find the relevant characters.

    strcpy (retstr, mystr);
    lastdot = strrchr (retstr, dot);
    lastsep = (sep == 0) ? NULL : strrchr (retstr, sep);

    // If it has an extension separator.

    if (lastdot != NULL) {
        // and it's before the extenstion separator.

        if (lastsep != NULL) {
            if (lastsep < lastdot) {
                // then remove it.

                *lastdot = '\0';
            }
        } else {
            // Has extension separator with no path separator.

            *lastdot = '\0';
        }
    }

    // Return the modified string.

    return retstr;
  }


int read_header(EPR_ELogLevel log_level, const char *infile, struct PRM * prm, state_vector *sv, int * n_state_vectors)
  {
  EPR_SProductId* 		product_id;
  EPR_SDatasetId* 		MDS1_SQ_ADS;
  EPR_SDatasetId* 		MAIN_PROC_PM_ID;
  EPR_SDatasetId* 		DOP_CENTROID_COEFFS_ADS;
  EPR_SDatasetId* 		CHIRP_PARAMS_ADS;
  EPR_SDatasetId* 		GEOLOCATION_GRID_ADS;
  EPR_SRecord*    		mph;
  EPR_SRecord*    		sph;
  EPR_SRecord*    		rec0 = NULL;
  EPR_SRecord*    		rec1 = NULL;
  EPR_SRecord*    		rec2 = NULL;
  EPR_SRecord*    		rec3 = NULL;
  EPR_SRecord*    		rec4 = NULL;
//  EPR_SField*     		azitime0_field;
  EPR_EErrCode    		err_code=e_err_none;
  int             		status;
  char 				tmp_c[200];
  char 				filename[200];
  double 			rbias = 0.;
  double 			tbias = 0.;
  int 				tmp_i;
  double 			c_speed = 299792458.0;
  char*				tmp_string;  

  EPR_SField      		num_looks_range_field;
  uint				num_looks_range=0;
  EPR_SField      		range_samp_rate_field;
  double 			range_samp_rate;
  EPR_SField      		product_field;
  const char * 			product_name;
  int				SC_identity=0;
  const char 			dot = '.';
  char *			Str_SC_identity;
  EPR_SField      		radar_freq_field;
  double 			radar_freq;
  EPR_SField      		pulse_length_field;
  double			pulse_length;
  EPR_SField      		look_bandwidth_field;
  double			look_bandwidth;
  EPR_SField      		output_mean_field;
  double 			output_mean_i;
  double 			output_mean_q;
  EPR_SField      		line_time_interval_field;
  double 			line_time_interval;
  EPR_SField      		slant_range_time_field;
  double 			slant_range_time_ns;
  double			slant_range_time_s;
  double			near_range;
  EPR_SField                    geogrid_slant_range_time_field1;
  EPR_SField                    geogrid_slant_range_time_field2;
  double                        geogrid_slant_range_time1_ns;
  double                        geogrid_slant_range_time2_ns;
  EPR_SField      		pass_field;
  const char *			pass;
  char				orbdir[1];
  EPR_SField      		start_time_field;
  const char *			start_time;
  char 				s_name[200];
  char 				s_out[200];
  char 				year[5];
  char 				month[4];
  char 				day[3];
  EPR_SField      		stop_time_field;
  const char *			stop_time;
  EPR_SField      		num_samples_per_line_field;
  int 				num_samples_per_line;
  EPR_SField      		num_output_lines_field;
  int 				num_output_lines;
  EPR_SField      		orbit_state_vector_time_field;
  EPR_SField      		orbit_state_vector_value_field;
  char				orbit_state_vector_time_field_name[50];
  char				orbit_state_vector_value_field_name[50];
  const EPR_STime *		orbit_state_vector_time_value;
  int				year_for_state_vectors;
  int				day_for_state_vectors;
  EPR_SField                    zero_doppler_time_field;
  const EPR_STime *             zero_doppler_time_mjd;
  EPR_SField                    attach_flag_field;
  int                           attach_flag_value;
  EPR_SField                    dop_coef_field;
  double                        dop_coef_value_D0;
  double                        dop_coef_value_D1;
  double                        dop_coef_value_D2;
  double                        dop_coef_value_D3;
  double                        dop_coef_value_D4;
  EPR_SField                    dop_conf_field;
  double                        dop_conf_value;
  EPR_SField                    dop_thresh_flag_field;
  int                           dop_thresh_flag_value;
  int 				q;
  double			dr;


  // define some of the variables
  prm->first_line = 1;
  prm->st_rng_bin = 1;
  prm->rshift = 0;
  prm->ashift = 0;
  prm->sub_int_r = 0.0;
  prm->sub_int_a = 0.0;
  prm->stretch_r = 0.0;
  prm->stretch_a = 0.0;
  prm->a_stretch_r = 0.0;
  prm->a_stretch_a = 0.0;
  prm->first_sample = 1;
  strasign(prm->dtype,"a",0,0);


  //Routine to save the filename without extension in the variable filename
  strcpy(filename,infile);
  tmp_string=remove_ext(filename, '.', '/');
  strcpy(filename,tmp_string);
  free(tmp_string);


  /* Initialize the API. Set log-level and use default log-output (stdout) */
  /* How to set error logging without exiting? */
  /*status = epr_init_api(e_log_debug, epr_log_message, epr_log_message);*/
  status = epr_init_api(log_level, epr_log_message, NULL);
  if (status != 0)
    {
    printf("read_header: fatal error in epr_init_api\n");
    printf("exiting.\n");
    return 1;
    };

  /* Open the product; an argument is a path to product data file */
  /* PRODUCT dataset record field element */

  if (log_level == e_log_debug)
  {
  printf("\nOpening product.\n");
  printf("-------------------------------------------------\n");
  }

  product_id = epr_open_product(infile);
  err_code   = epr_get_last_err_code();
  if (err_code != e_err_none) 
    {
    printf("read_header: fatal error in epr_open_product\n");
    printf("exiting.\n");
    return 1;
    }

  /* Read fields that I want */
  /* product DATASET record field element */
  /* --- MAIN PRODUCT HEADER ---------------------------------------- */
  if (log_level == e_log_debug)
  {
  printf("\nTrying to read record: MPH\n");
  printf("-------------------------------------------------\n");
  }
  mph      = epr_get_mph(product_id);
  err_code = epr_get_last_err_code();
  if (err_code == e_err_none) 
    {
    if (log_level == e_log_debug)
      epr_print_record(mph,  stdout);
    }
  else
    {
    printf("read_header: likely fatal error in epr_get_mph\n");
    epr_clear_err();
    printf("exiting.\n");
    return 1;
    }

  /* --- SECOND PRODUCT HEADER -------------------------------------- */
  if (log_level == e_log_debug)
  {
  printf("\nTrying to read record: SPH\n");
  printf("-------------------------------------------------\n");
  }
  sph      = epr_get_sph(product_id);
  err_code = epr_get_last_err_code();
  if (err_code == e_err_none) 
    {
    if (log_level == e_log_debug)
      epr_print_record(sph,  stdout);
    }
  else
    {
    printf("read_header: likely fatal error in epr_get_sph\n");
    epr_clear_err();
    printf("exiting.\n");
    return 1;
    }

  /* dsd  = epr_get_dsd(product_id); */
  /* epr_print_record(dsd,  stdout);*/

  /* --- MDS1_SQ_ADS ------------------------------------------------ */
  if (log_level == e_log_debug)
  {
  printf("\nTrying to read record: MDS1_SQ_ADS\n");
  printf("-------------------------------------------------\n");
  }

  MDS1_SQ_ADS = epr_get_dataset_id(product_id, "MDS1_SQ_ADS");
  err_code    = epr_get_last_err_code();
  if (err_code == e_err_none)
    {
    rec0 = epr_read_record(MDS1_SQ_ADS, 0, NULL);
    err_code = epr_get_last_err_code();
    if (err_code == e_err_none)
      {
      if (log_level == e_log_debug)
        epr_print_record(rec0, stdout);
      }
    else
      {
      printf("read_header: non-fatal error in epr_read_record MDS1_SQ_ADS\n");
      printf("could not read this record.  may not be a problem.\n");
      epr_clear_err();
      }
    }
  else
    {
    printf("read_header: error in epr_get_dataset_id MDS1_SQ_ADS\n");
    printf("could not read this record.  may not be a problem.\n");
    epr_clear_err();
    }

  /* --- MAIN_PROC_PM_ID -------------------------------------------- */
    if (log_level == e_log_debug)
  {
  printf("\nTrying to read record: MAIN_PROCESSING_PARAMS_ADS\n");
  printf("-------------------------------------------------\n");
  }

  MAIN_PROC_PM_ID = epr_get_dataset_id(product_id, "MAIN_PROCESSING_PARAMS_ADS");
  err_code    = epr_get_last_err_code();
  if (err_code == e_err_none)
    {
    rec1 = epr_read_record(MAIN_PROC_PM_ID,         0, NULL);
    err_code = epr_get_last_err_code();
    if (err_code == e_err_none)
      {
      if (log_level == e_log_debug)
        epr_print_record(rec1, stdout);
      }
    else
      {
      printf("read_header: error in epr_read_record MAIN_PROCESSING_PARAMS_ADS\n");
      printf("could not read this record.  may not be a problem.\n");
      epr_clear_err();
      }
    }
  else
    {
    printf("read_header: error in epr_get_dataset_id MAIN_PROCESSING_PARAMS_ADS\n");
    printf("could not read this record.  may not be a problem.\n");
    epr_clear_err();
    }

  /* --- DOP_CENTROID_COEFFS_ADS ------------------------------------ */
  if (log_level == e_log_debug)
  {
  printf("\nTrying to read record: DOP_CENTROID_COEFFS_ADS\n");
  printf("-------------------------------------------------\n");
  }

  DOP_CENTROID_COEFFS_ADS = epr_get_dataset_id(product_id, "DOP_CENTROID_COEFFS_ADS");
  err_code    = epr_get_last_err_code();
  if (err_code == e_err_none)
    {
    rec2 = epr_read_record(DOP_CENTROID_COEFFS_ADS, 0, NULL);
    err_code = epr_get_last_err_code();
    if (err_code == e_err_none)
      {
      if (log_level == e_log_debug)
        epr_print_record(rec2, stdout);
      }
    else
      {
      printf("read_header: error in epr_read_record DOP_CENTROID_COEFFS_ADS\n");
      printf("could not read this record.  may not be a problem.\n");
      epr_clear_err();
      }
    }
  else
    {
    printf("read_header: error in epr_get_dataset_id DOP_CENTROID_COEFFS_ADS\n");
    printf("could not read this record.  may not be a problem.\n");
    epr_clear_err();
    }

 /* --- CHIRP_PARAMS_ADS ------------------------------------------- */
  if (log_level == e_log_debug)
  {
  printf("\nTrying to read record: CHIRP_PARAMS_ADS\n");
  printf("-------------------------------------------------\n");
  }

  /* it seems ESRIN does not include this.  For now do not read it at all
  Problems with BAM data, report Z.Perksi
  #%// Bert Kampes, 12-May-2004 */
  CHIRP_PARAMS_ADS = epr_get_dataset_id(product_id, "CHIRP_PARAMS_ADS");
  err_code         = epr_get_last_err_code();
  if (err_code == e_err_none)
    {
    rec3 = epr_read_record(CHIRP_PARAMS_ADS,        0, NULL);
    err_code = epr_get_last_err_code();
    if (err_code == e_err_none)
      {
      if (log_level == e_log_debug)
        epr_print_record(rec3, stdout);
      }
    else
      {
      printf("read_header: error in epr_read_record CHIRP_PARAMS_ADS\n");
      printf("could not read this record.  may not be a problem.\n");
      epr_clear_err();
      }
    }
  else
    {
    printf("read_header: error in epr_get_dataset_id CHIRP_PARAMS_ADS\n");
    printf("could not read this record.  may not be a problem.\n");
    epr_clear_err();
    }



  /* --- GEOLOCATION_GRID_ADS --------------------------------------- */
  if (log_level == e_log_debug)
  {
  printf("\nTrying to read record: GEOLOCATION_GRID_ADS\n");
  printf("-------------------------------------------------\n");
  }

  GEOLOCATION_GRID_ADS = epr_get_dataset_id(product_id, "GEOLOCATION_GRID_ADS");
  err_code         = epr_get_last_err_code();
  if (err_code == e_err_none)
    {
    rec4 = epr_read_record(GEOLOCATION_GRID_ADS,    0, NULL);
    err_code = epr_get_last_err_code();
    if (err_code == e_err_none)
      {
      if (log_level == e_log_debug)
        epr_print_record(rec4, stdout);
      }
    else
      {
      printf("read_header: error in epr_read_record GEOLOCATION_GRID_ADS\n");
      printf("could not read this record.  may not be a problem.\n");
      epr_clear_err();
      }
    }
  else
    {
    printf("read_header: error in epr_get_dataset_id GEOLOCATION_GRID_ADS\n");
    printf("could not read this record.  may not be a problem.\n");
    epr_clear_err();
    }


  /* product dataset record FIELD element */
  /*
  azitime0_field = epr_get_field(rec1, "first_zero_doppler_time");
  epr_print_field(azitime0_field, stdout);
  */

  /*Set all values in the prm struct */


  num_looks_range_field	= *(epr_get_field(rec1, "num_looks_range"));
  num_looks_range		= epr_get_field_elem_as_uint(&num_looks_range_field, 0);
  prm->nlooks = num_looks_range;  	//nlooks


  range_samp_rate		= 0.0;
  range_samp_rate_field	= *(epr_get_field(rec1, "range_samp_rate"));
  range_samp_rate 		= epr_get_field_elem_as_double(&range_samp_rate_field, 0);
  prm->fs = range_samp_rate;		//rng_samp_rate


  product_field			= *(epr_get_field(mph, "PRODUCT"));
  product_name			= epr_get_field_elem_as_str(&product_field);
  Str_SC_identity		= strrchr(product_name, dot); //Returns the last part of the original file name for identification of the platform

  //printf("Str_SC_identity %s\n",Str_SC_identity);

  if(strcasecmp(Str_SC_identity, ".N1") == 0)
	{
	SC_identity=6; 
	}
  else if(strcasecmp(Str_SC_identity, ".E1") == 0)
	{
	SC_identity=1;
	}
  else if(strcasecmp(Str_SC_identity, ".E2") == 0)
	{
	SC_identity=2;
	}
  else {
	SC_identity=0;
	}

  prm->SC_identity = SC_identity; /* (1)-ERS1 (2)-ERS2 (3)-Radarsat (4)-Envisat (5)-ALOS (6)-Envisat_SLC  (7)-TSX (8)-CSK (9)-RS2 (10)-Sentinel-1a*/



  radar_freq			= 0.0;
  radar_freq_field		= *(epr_get_field(rec1, "radar_freq"));
  radar_freq 			= epr_get_field_elem_as_double(&radar_freq_field, 0);
  prm->lambda 			= c_speed/radar_freq; //radar_wavelength



  pulse_length			= 0.0;
  pulse_length_field		= *(epr_get_field(rec1, "image_parameters.tx_pulse_len_value"));
  pulse_length			= epr_get_field_elem_as_double(&pulse_length_field, 0);
  //printf("pulse_length %f\n",pulse_length);


  look_bandwidth		= 0.0;
  look_bandwidth_field		= *(epr_get_field(rec1, "bandwidth.look_bw_range"));
  look_bandwidth		= epr_get_field_elem_as_double(&look_bandwidth_field, 0);
  //printf("look_bandwidth %f\n",look_bandwidth);
  prm->chirp_slope 		= look_bandwidth / pulse_length;
  prm->pulsedur 		= pulse_length;
 
  

  output_mean_i			= 0.0;
  output_mean_q			= 0.0;
  output_mean_field		= *(epr_get_field(rec0, "output_mean"));
  output_mean_i			= epr_get_field_elem_as_double(&output_mean_field, 0);
  output_mean_q			= epr_get_field_elem_as_double(&output_mean_field, 1);
  //printf("output_mean_i %f\noutput_mean_q %f\n",output_mean_i,output_mean_q);
  prm->xmi 			= output_mean_i; //I_mean
  prm->xmq 			= output_mean_q; //Q_mean

  line_time_interval		= 0.0;
  line_time_interval_field	= *(epr_get_field(rec1, "line_time_interval"));
  line_time_interval 		= epr_get_field_elem_as_double(&line_time_interval_field, 0);   
  prm->prf 			= 1/line_time_interval; //Pulse Repetition Frequency


  slant_range_time_ns		= 0.0;
  slant_range_time_s		= 0.0;
  near_range			= 0.0;
//  slant_range_time_field	= *(epr_get_field(rec2, "slant_range_time"));
//  slant_range_time_ns 		= epr_get_field_elem_as_double(&slant_range_time_field, 0);   



  geogrid_slant_range_time_field1        = *(epr_get_field(rec4, "first_line_tie_points.slant_range_times"));
  geogrid_slant_range_time_field2        = *(epr_get_field(rec4, "last_line_tie_points.slant_range_times"));
  geogrid_slant_range_time1_ns           = epr_get_field_elem_as_double(&geogrid_slant_range_time_field1, 0);
  geogrid_slant_range_time2_ns           = epr_get_field_elem_as_double(&geogrid_slant_range_time_field2, 0);

  slant_range_time_ns		= (geogrid_slant_range_time1_ns + geogrid_slant_range_time2_ns) / 2;


  slant_range_time_s            = slant_range_time_ns * 0.000000001;
  near_range                    = slant_range_time_s * c_speed / 2;

  /* convert to pixel space */
  dr				= 0.5 * c_speed / prm->fs;  
  

  /* make a bias correction the range for ENVI and ERS*/
  if(SC_identity == 1) {
     rbias = -10.52 * dr;
  }
  else if(SC_identity == 2) {
     rbias = -5.3 * dr;
  }
  else if(SC_identity == 6) {
     rbias = -0.8 * dr;
  }
  else {
    printf(" SC_identity out of range ");
  }

  prm->near_range 		= near_range + rbias; //near_range


  /* Read Doppler Processing Parameters */

  /* Zero Doppler azimuth time at which estimate applies */
  zero_doppler_time_field       = *(epr_get_field(rec2, "zero_doppler_time"));
  zero_doppler_time_mjd         = epr_get_field_elem_as_mjd(&zero_doppler_time_field);

  /* Attachment Flag (always set to zero for this ADSR */
  attach_flag_field		= *(epr_get_field(rec2, "attach_flag"));  
  attach_flag_value		= epr_get_field_elem_as_double(&attach_flag_field, 0);

  /* 2-way slant range time origin (t0) */
  slant_range_time_ns		= 0;
  slant_range_time_field        = *(epr_get_field(rec2, "slant_range_time"));
  slant_range_time_ns           = epr_get_field_elem_as_double(&slant_range_time_field, 0);

  /* Doppler centroid coefficients as a function of slant range time: D0, D1, D2, D3, and D4. */
  dop_coef_field		= *(epr_get_field(rec2, "dop_coef"));
  dop_coef_value_D0		= epr_get_field_elem_as_double(&dop_coef_field, 0);
  dop_coef_value_D1             = epr_get_field_elem_as_double(&dop_coef_field, 1);
  dop_coef_value_D2             = epr_get_field_elem_as_double(&dop_coef_field, 2);
  dop_coef_value_D3             = epr_get_field_elem_as_double(&dop_coef_field, 3);
  dop_coef_value_D4             = epr_get_field_elem_as_double(&dop_coef_field, 4);

  /* Doppler Centroid Confidence Measure */
  dop_conf_field		= *(epr_get_field(rec2, "dop_conf"));
  dop_conf_value		= epr_get_field_elem_as_double(&dop_conf_field, 0);

  /* Doppler Confidence Below Threshold Flag */
  dop_thresh_flag_field		= *(epr_get_field(rec2, "dop_thresh_flag"));
  dop_thresh_flag_value		= epr_get_field_elem_as_double(&dop_thresh_flag_field, 0);

  /*
  printf("\nINFO:\nDoppler Centroid Coefficients ADSR\n");
  printf("Zero Doppler azimuth time at which estimate applies: d=%d (days), j=%d (seconds), m=%d (microseconds)\n",zero_doppler_time_mjd->days,zero_doppler_time_mjd->seconds,zero_doppler_time_mjd->microseconds);
  printf("Attachment Flag (always set to zero for this ADSR): %d\n",attach_flag_value);
  printf("2-way slant range time origin (t0): %f (ns)\n",slant_range_time_ns);
  printf("Doppler centroid coefficients as a function of slant range time, D0, D1, D2, D3, and D4: D0=%f (Hz), D1=%f (Hz/s), D2=%f (Hz/s2), D3=%f (Hz/s3), D4=%f (Hz/s4)\n",dop_coef_value_D0,dop_coef_value_D1,dop_coef_value_D2,dop_coef_value_D3,dop_coef_value_D4);
  printf("Doppler Centroid Confidence Measure: %f\n",dop_conf_value);

  printf("Value between 0 and 1, 0 = poorest confidence, 1= highest confidence\n");
  printf("If multiple Doppler Centroid estimates were performed, this value is the lowest confidence value attained.\n");

  printf("Doppler Confidence Below Threshold Flag: %d\n",dop_thresh_flag_value);

  printf("0 = confidence above threshold, Doppler Centroid calculated from data\n");
  printf("1 = confidence below threshold, Doppler Centroid calculated from orbit parameters\n");

  printf("\n");
  */


  /* Set fd1 to the value of D0 for now and make the assumption that it is close enough, we'll make the precise calculations later */
  /* David noted that this value is about twice of what it is supposed to be, divide by 2 for now */
  //prm->fd1 = dop_coef_value_D0/2;
  prm->fd1 = 0;



  /* End Read Doppler Processing Parameters */

  prm->ra = 6378137.00; //equatorial_radius
  prm->rc = 6356752.31; //polar_radius  

  pass_field			= *(epr_get_field(sph, "PASS"));
  pass				= epr_get_field_elem_as_str(&pass_field);
  strncpy(orbdir, pass, 1);
  strasign(prm->orbdir,orbdir,0,0); //orbdir

  strasign(prm->lookdir,"R",0,0);

  strcpy(tmp_c,filename);
  strcat(tmp_c,".SLC");
  strcpy(prm->input_file,tmp_c);
    
  strcpy(tmp_c,filename);
  strcat(tmp_c,".LED");
  strcpy(prm->led_file,tmp_c);
    
  strcpy(tmp_c,filename);
  strcat(tmp_c,".SLC");
  strcpy(prm->SLC_file,tmp_c);

  prm->SLC_scale = 1.0;
 

  start_time_field		= *(epr_get_field(sph, "FIRST_LINE_TIME"));
  start_time			= epr_get_field_elem_as_str(&start_time_field);
  strcpy(tmp_c,start_time);
  for (q=0; q<strlen(tmp_c); q++)
	{
	if(tmp_c[q] == ' ')
		tmp_c[q] = 'T';
	}

  /* Ugly routine to convert the date format to be the same as in Sentinel-1 data */ 
  year[4] = '\0';  
  month[3] = '\0';
  day[2] = '\0';

  strncpy(year,&tmp_c[7],4);  
  strncpy(month,&tmp_c[3],3);  
  strncpy(day,&tmp_c[0],2);  
  //Convert month to number
  if(monthtonum(month) < 10)
  	sprintf(tmp_c,"%s-0%d-%s%s",year,monthtonum(month),day,&tmp_c[11]);
  else
  	sprintf(tmp_c,"%s-%d-%s%s",year,monthtonum(month),day,&tmp_c[11]);

  cat_nums(s_name,tmp_c);
  str_date2JD(s_out, s_name);

  /* make a time bias correction for ENVI and ERS*/
  if(SC_identity == 1) {
     tbias = 19.2501/prm->prf/86400.;
  }
  else if(SC_identity == 2) {
     tbias = 0.;
  }
  else if(SC_identity == 6) {
     tbias = -4.25/prm->prf/86400.;
  }
  else {
    printf(" SC_identity out of range ");
  }


  prm->clock_start = atof(s_out)+1+tbias;
  tmp_c[4] = '\0';
  prm->SC_clock_start = prm->clock_start + 1000.*atof(tmp_c);
		
  stop_time_field		= *(epr_get_field(sph, "LAST_LINE_TIME"));
  stop_time			= epr_get_field_elem_as_str(&stop_time_field);
  strcpy(tmp_c,stop_time); 
  for (int q=0; q<strlen(tmp_c); q++)
	{
	if(tmp_c[q] == ' ')
		tmp_c[q] = 'T';
	}

  /* Ugly routine to convert the date format to be the same as in Sentinel-1 data */ 
  year[4] = '\0';  
  month[3] = '\0';
  day[2] = '\0';

  strncpy(year,&tmp_c[7],4);  
  strncpy(month,&tmp_c[3],3);  
  strncpy(day,&tmp_c[0],2);  
  //Save year and day for use in state vectors
  year_for_state_vectors = str2double(year);
  

  //Convert month to number
  if(monthtonum(month) < 10)
  	sprintf(tmp_c,"%s-0%d-%s%s",year,monthtonum(month),day,&tmp_c[11]);
  else
  	sprintf(tmp_c,"%s-%d-%s%s",year,monthtonum(month),day,&tmp_c[11]);
  
  cat_nums(s_name,tmp_c);
  str_date2JD(s_out, s_name);
  //Save day for use in state vectors
  day_for_state_vectors = atof(s_out)+1;


  //Why doesn't str2double work here?
  prm->clock_stop = atof(s_out)+1;
  tmp_c[4] = '\0';
  prm->SC_clock_stop = prm->clock_stop + 1000.*atof(tmp_c);


  strasign(prm->iqflip,"n",0,0); //Flip_iq
  strasign(prm->deskew,"n",0,0); //deskew
  strasign(prm->offset_video,"n",0,0);

  num_samples_per_line		= 0.0;
  tmp_i				= 0;
  num_samples_per_line_field	= *(epr_get_field(rec1, "num_samples_per_line"));
  num_samples_per_line 	= (int)epr_get_field_elem_as_double(&num_samples_per_line_field, 0);
  /* This does not work with ERS/ENVISAT data
  tmp_i 			= num_samples_per_line - num_samples_per_line%4; */
  tmp_i 			= num_samples_per_line;
  prm->bytes_per_line		= tmp_i*4;
 

  prm->good_bytes = prm->bytes_per_line;
  prm->caltone = 0.0;
  prm->pctbwaz = 0.0; //rm_az_band
  prm->pctbw = 0.2; //rm_rng_band
  prm->rhww = 1.0; //rng_spec_wgt
  strasign(prm->srm,"0",0,0); //scnd_rng_mig
  prm->az_res = 0.0;
  //prm.antenna_side = -1;
  prm->fdd1 = 0.0;
  prm->fddd1 = 0.0;

  num_output_lines		= 0.0;
  num_output_lines_field	= *(epr_get_field(rec1, "num_output_lines"));
  num_output_lines 		= (int)epr_get_field_elem_as_double(&num_output_lines_field, 0);
  tmp_i				= 0;
  tmp_i				= num_output_lines%4;
  prm->num_lines 		= num_output_lines - tmp_i;
 
  prm->nrows = prm->num_lines;
  prm->num_valid_az = prm->num_lines;
  prm->num_patches = 1;
  prm->num_rng_bins = prm->bytes_per_line/4;
  prm->chirp_ext	= 0;
    
  printf("PRM set for Image File...\n");


 //Read orbit state vectors, there should always be five of them
  for(int n_orbit_state_vector=1;n_orbit_state_vector<6;n_orbit_state_vector++)
{
	// Time values
	sprintf(orbit_state_vector_time_field_name,"orbit_state_vectors.%d.state_vect_time_1",n_orbit_state_vector);

	orbit_state_vector_time_field           = *(epr_get_field(rec1, orbit_state_vector_time_field_name));
	orbit_state_vector_time_value           = epr_get_field_elem_as_mjd(&orbit_state_vector_time_field);

	sv[n_orbit_state_vector-1].yr           = year_for_state_vectors;
	sv[n_orbit_state_vector-1].jd           = day_for_state_vectors;
	sv[n_orbit_state_vector-1].sec  = orbit_state_vector_time_value->seconds;

	// x pos
	sprintf(orbit_state_vector_value_field_name,"orbit_state_vectors.%d.x_pos_1",n_orbit_state_vector);
	orbit_state_vector_value_field  = *(epr_get_field(rec1, orbit_state_vector_value_field_name));
	sv[n_orbit_state_vector-1].x            = (int)epr_get_field_elem_as_double(&orbit_state_vector_value_field, 0)/100;
	// y pos
	sprintf(orbit_state_vector_value_field_name,"orbit_state_vectors.%d.y_pos_1",n_orbit_state_vector);
	orbit_state_vector_value_field  = *(epr_get_field(rec1, orbit_state_vector_value_field_name));
	sv[n_orbit_state_vector-1].y            = (int)epr_get_field_elem_as_double(&orbit_state_vector_value_field, 0)/100;
	// z pos
	sprintf(orbit_state_vector_value_field_name,"orbit_state_vectors.%d.z_pos_1",n_orbit_state_vector);
	orbit_state_vector_value_field  = *(epr_get_field(rec1, orbit_state_vector_value_field_name));
	sv[n_orbit_state_vector-1].z            = (int)epr_get_field_elem_as_double(&orbit_state_vector_value_field, 0)/100;

	// x vel
	sprintf(orbit_state_vector_value_field_name,"orbit_state_vectors.%d.x_vel_1",n_orbit_state_vector);
	orbit_state_vector_value_field  = *(epr_get_field(rec1, orbit_state_vector_value_field_name));
	sv[n_orbit_state_vector-1].vx           = (int)epr_get_field_elem_as_double(&orbit_state_vector_value_field, 0)/100000;
	// y vel
	sprintf(orbit_state_vector_value_field_name,"orbit_state_vectors.%d.y_vel_1",n_orbit_state_vector);
	orbit_state_vector_value_field  = *(epr_get_field(rec1, orbit_state_vector_value_field_name));
	sv[n_orbit_state_vector-1].vy           = (int)epr_get_field_elem_as_double(&orbit_state_vector_value_field, 0)/100000;
	// z vel
	sprintf(orbit_state_vector_value_field_name,"orbit_state_vectors.%d.z_vel_1",n_orbit_state_vector);
	orbit_state_vector_value_field  = *(epr_get_field(rec1, orbit_state_vector_value_field_name));
	sv[n_orbit_state_vector-1].vz           = (int)epr_get_field_elem_as_double(&orbit_state_vector_value_field, 0)/100000;

	*n_state_vectors = n_orbit_state_vector;
	}

  printf("LED set for Image File...\n");


  /* Close product_id and release rest of the allocated memory */
  printf("\n");
  epr_close_product(product_id);
  /* Closes product reader API, release all allocated resources */
  epr_close_api();


	return 0;
  }


// Converts month to number 1 - 12
int monthtonum(char * szMonth)
  {
  const char * months[12] = {"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};
  int nMonth = 0;
  int i;

  for (i = 0; i < 12; ++i)
   {
    if(strncasecmp(szMonth, months[i], 3) == 0)
      nMonth = i+1;
   }
  	return nMonth;
  }
