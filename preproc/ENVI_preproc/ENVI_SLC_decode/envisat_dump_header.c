/* BK: simple program to get header info for ENVISAT
*   13-Jun-2003
*   simply dumps all header info.  a unix script will be used to 
*   convert this to doris input
    $Id: envisat_dump_header.c,v 1.4 2004/05/13 18:13:39 kampes Exp $
*/
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "epr_api-2.3/src/epr_api.h"
#if defined(WIN32) && defined(_DEBUG)
#include <crtdbg.h>
#endif /* if defined(WIN32) && defined(_DEBUG) */


int main(int argc, char** argv) 
  {
  EPR_SProductId* product_id;
  EPR_SDatasetId* MDS1_SQ_ADS;
  EPR_SDatasetId* MAIN_PROC_PM_ID;
  EPR_SDatasetId* DOP_CENTROID_COEFFS_ADS;
  EPR_SDatasetId* CHIRP_PARAMS_ADS;
  EPR_SDatasetId* GEOLOCATION_GRID_ADS;
  EPR_SRecord*    mph;
  EPR_SRecord*    sph;
  EPR_SRecord*    rec0;
  EPR_SRecord*    rec1;
  EPR_SRecord*    rec2;
  EPR_SRecord*    rec3;
  EPR_SRecord*    rec4;
//  EPR_SField*     azitime0_field;
  EPR_EErrCode    err_code=e_err_none;
  int             status;
  const char*     product_file_path;
  /* const char* datasetname; */



  /* --- Handle input ----------------------------------------- */
  printf("+-----------------------------------------------+\n");
  printf("|  This is envisat_dump_header by Bert Kampes   |\n");
  printf("+-----------------------------------------------+\n");
  if (argc <= 1)
    {
    printf("Usage: envisat_dump_header envisat-product\n");
    printf("  where envisat-product is the input filename\n");
    printf("Example:\n");
    printf("  envisat_dump_header ASA_IMS_1PNDPA20021025_175208_000000162010_00356_03416_0005.N1\n\n");
    exit(1);
    }
  product_file_path = argv[1];
  printf("\nInitializing api.\n");
  printf("-------------------------------------------------\n");

  /* Initialize the API. Set log-level to DEBUG and use default log-output (stdout) */
  /* How to set error logging without exiting? */
  /*status = epr_init_api(e_log_debug, epr_log_message, epr_log_message);*/
  status = epr_init_api(e_log_debug, epr_log_message, NULL);
  if (status != 0)
    {
    printf("envisat_dump_header: fatal error in epr_init_api\n");
    printf("exiting.\n");
    exit(1);
    };

  /* Open the product; an argument is a path to product data file */
  /* PRODUCT dataset record field element */
  printf("\nOpening product.\n");
  printf("-------------------------------------------------\n");
  product_id = epr_open_product(product_file_path);
  err_code   = epr_get_last_err_code();
  if (err_code != e_err_none) 
    {
    printf("envisat_dump_header: fatal error in epr_open_product\n");
    printf("exiting.\n");
    exit(1);
    }



  /* Read fields that I want for Doris */
  /* product DATASET record field element */
  /* --- MAIN PRODUCT HEADER ---------------------------------------- */
  printf("\nTrying to read record: MPH\n");
  printf("-------------------------------------------------\n");
  mph      = epr_get_mph(product_id);
  err_code = epr_get_last_err_code();
  if (err_code == e_err_none) 
    {
    epr_print_record(mph,  stdout);
    }
  else
    {
    printf("envisat_dump_header: likely fatal error in epr_get_mph\n");
    printf("continuing.\n");
    epr_clear_err();
    }



  /* --- SECOND PRODUCT HEADER -------------------------------------- */
  printf("\nTrying to read record: SPH\n");
  printf("-------------------------------------------------\n");
  sph      = epr_get_sph(product_id);
  err_code = epr_get_last_err_code();
  if (err_code == e_err_none) 
    {
    epr_print_record(sph,  stdout);
    }
  else
    {
    printf("envisat_dump_header: likely fatal error in epr_get_sph\n");
    epr_clear_err();
    }


  /* dsd  = epr_get_dsd(product_id); */
  /* epr_print_record(dsd,  stdout);*/

  /* --- MDS1_SQ_ADS ------------------------------------------------ */
  printf("\nTrying to read record: MDS1_SQ_ADS\n");
  printf("-------------------------------------------------\n");
  MDS1_SQ_ADS = epr_get_dataset_id(product_id, "MDS1_SQ_ADS");
  err_code    = epr_get_last_err_code();
  if (err_code == e_err_none)
    {
    rec0 = epr_read_record(MDS1_SQ_ADS, 0, NULL);
    err_code = epr_get_last_err_code();
    if (err_code == e_err_none)
      {
      epr_print_record(rec0, stdout);
      }
    else
      {
      printf("envisat_dump_header: non-fatal error in epr_read_record MDS1_SQ_ADS\n");
      printf("could not read this record.  may not be a problem.\n");
      epr_clear_err();
      }
    }
  else
    {
    printf("envisat_dump_header: error in epr_get_dataset_id MDS1_SQ_ADS\n");
    printf("could not read this record.  may not be a problem.\n");
    epr_clear_err();
    }



  /* --- MAIN_PROC_PM_ID -------------------------------------------- */
  printf("\nTrying to read record: MAIN_PROCESSING_PARAMS_ADS\n");
  printf("-------------------------------------------------\n");
  MAIN_PROC_PM_ID = epr_get_dataset_id(product_id, "MAIN_PROCESSING_PARAMS_ADS");
  err_code    = epr_get_last_err_code();
  if (err_code == e_err_none)
    {
    rec1 = epr_read_record(MAIN_PROC_PM_ID,         0, NULL);
    err_code = epr_get_last_err_code();
    if (err_code == e_err_none)
      {
      epr_print_record(rec1, stdout);
      }
    else
      {
      printf("envisat_dump_header: error in epr_read_record MAIN_PROCESSING_PARAMS_ADS\n");
      printf("could not read this record.  may not be a problem.\n");
      epr_clear_err();
      }
    }
  else
    {
    printf("envisat_dump_header: error in epr_get_dataset_id MAIN_PROCESSING_PARAMS_ADS\n");
    printf("could not read this record.  may not be a problem.\n");
    epr_clear_err();
    }



  /* --- DOP_CENTROID_COEFFS_ADS ------------------------------------ */
  printf("\nTrying to read record: DOP_CENTROID_COEFFS_ADS\n");
  printf("-------------------------------------------------\n");
  DOP_CENTROID_COEFFS_ADS = epr_get_dataset_id(product_id, "DOP_CENTROID_COEFFS_ADS");
  err_code    = epr_get_last_err_code();
  if (err_code == e_err_none)
    {
    rec2 = epr_read_record(DOP_CENTROID_COEFFS_ADS, 0, NULL);
    err_code = epr_get_last_err_code();
    if (err_code == e_err_none)
      {
      epr_print_record(rec2, stdout);
      }
    else
      {
      printf("envisat_dump_header: error in epr_read_record DOP_CENTROID_COEFFS_ADS\n");
      printf("could not read this record.  may not be a problem.\n");
      epr_clear_err();
      }
    }
  else
    {
    printf("envisat_dump_header: error in epr_get_dataset_id DOP_CENTROID_COEFFS_ADS\n");
    printf("could not read this record.  may not be a problem.\n");
    epr_clear_err();
    }



  /* --- CHIRP_PARAMS_ADS ------------------------------------------- */
  printf("\nTrying to read record: CHIRP_PARAMS_ADS\n");
  printf("-------------------------------------------------\n");
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
      epr_print_record(rec3, stdout);
      }
    else
      {
      printf("envisat_dump_header: error in epr_read_record CHIRP_PARAMS_ADS\n");
      printf("could not read this record.  may not be a problem.\n");
      epr_clear_err();
      }
    }
  else
    {
    printf("envisat_dump_header: error in epr_get_dataset_id CHIRP_PARAMS_ADS\n");
    printf("could not read this record.  may not be a problem.\n");
    epr_clear_err();
    }



  /* --- GEOLOCATION_GRID_ADS --------------------------------------- */
  printf("\nTrying to read record: GEOLOCATION_GRID_ADS\n");
  printf("-------------------------------------------------\n");
  GEOLOCATION_GRID_ADS = epr_get_dataset_id(product_id, "GEOLOCATION_GRID_ADS");
  err_code         = epr_get_last_err_code();
  if (err_code == e_err_none)
    {
    rec4 = epr_read_record(GEOLOCATION_GRID_ADS,    0, NULL);
    err_code = epr_get_last_err_code();
    if (err_code == e_err_none)
      {
      epr_print_record(rec4, stdout);
      }
    else
      {
      printf("envisat_dump_header: error in epr_read_record GEOLOCATION_GRID_ADS\n");
      printf("could not read this record.  may not be a problem.\n");
      epr_clear_err();
      }
    }
  else
    {
    printf("envisat_dump_header: error in epr_get_dataset_id GEOLOCATION_GRID_ADS\n");
    printf("could not read this record.  may not be a problem.\n");
    epr_clear_err();
    }



  /* product dataset record FIELD element */
  /*
  azitime0_field = epr_get_field(rec1, "first_zero_doppler_time");
  epr_print_field(azitime0_field, stdout);
  */

  /* Close product_id and release rest of the allocated memory */
  printf("\n");
  epr_close_product(product_id);
  /* Closes product reader API, release all allocated resources */
  epr_close_api();
  printf("\n+-----------------------------------------------+\n");
  printf("|  Thank you for using envisat_dump_header      |\n");
  printf("+-----------------------------------------------+\n\n");

  return 0;
  }


