/*
 * $Id: epr_product.c,v 1.1.1.1 2004-10-28 19:22:23 norman Exp $
 *
 * Copyright (C) 2002 by Brockmann Consult (info@brockmann-consult.de)
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation. This program is distributed in the hope it will
 * be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "epr_api.h"
#include "epr_core.h"
#include "epr_string.h"
#include "epr_ptrarray.h"
#include "epr_swap.h"
#include "epr_field.h"
#include "epr_record.h"
#include "epr_dataset.h"
#include "epr_param.h"
#include "epr_dsd.h"
#include "epr_msph.h"
#include "epr_band.h"
#include "epr_bitmask.h"

#include "epr_dddb.h"

uint epr_compute_scene_width(const EPR_SProductId* product_id);
uint epr_compute_scene_height(const EPR_SProductId* product_id);

/*********************************** PRODUCT ***********************************/

/*
   Function:    epr_open_product
   Access:      public API
   Changelog:   2002/01/05  nf initial version
 */
/**
 * Opens the ENVISAT product file with the given file path
 */
EPR_SProductId* epr_open_product(const char* product_file_path) {
    EPR_SProductId* product_id = NULL;
    char message_buffer[80];
    int s_par;
    uint compare_ok = 0;

    epr_clear_err();
    if (!epr_check_api_init_flag()) {
        return NULL;
    }

    if (product_file_path == NULL) {
        epr_set_err(e_err_null_pointer,
                    "epr_open_product: product file path must not be NULL");
        return NULL;
    }

    product_id = (EPR_SProductId*) calloc(1, sizeof (EPR_SProductId));
    if (product_id == NULL) {
        epr_set_err(e_err_out_of_memory,
                    "epr_open_product: out of memory");
        return NULL;
    }
    product_id->magic = EPR_MAGIC_PRODUCT_ID;

    epr_assign_string(&product_id->file_path, product_file_path);

    if (product_id->file_path == NULL) {
        epr_set_err(e_err_out_of_memory,
                    "epr_open_product: out of memory");
        return NULL;
    }

    /* Convert to OS compatible path */
    epr_make_os_compatible_path(product_id->file_path);

    product_id->istream = fopen(epr_trim_string(product_id->file_path), "rb");
    if (product_id->istream == NULL) {
        if (errno == ENOENT) {
            epr_set_err(e_err_file_not_found,
                        "epr_open_product: file not found");
        } else {
            epr_set_err(e_err_file_access_denied,
                        "epr_open_product: file open failed");
        }
        return NULL;
    }

    epr_log(e_log_debug, "product opened:");
    epr_log(e_log_debug, epr_trim_string(product_id->file_path));

    /* Set file pointer to start of product identifier */
    if (fseek(product_id->istream, EPR_PRODUCT_ID_OFFSET, SEEK_SET) != 0) {
        epr_set_err(e_err_file_access_denied,
                    "epr_open_product: file seek failed");
        epr_close_product(product_id);
        return NULL;
    }

    if (fread(product_id->id_string,
              1,
              EPR_PRODUCT_ID_STRLEN,
              product_id->istream) != (uint) EPR_PRODUCT_ID_STRLEN) {
        epr_set_err(e_err_file_access_denied,
                    "epr_open_product: file read failed");
        epr_close_product(product_id);
        return NULL;
    }

    /* Product identifier filter*/
    if ((strncmp(EPR_ENVISAT_PRODUCT_MERIS, product_id->id_string, 3) != 0) &&
            (strncmp(EPR_ENVISAT_PRODUCT_ASAR,  product_id->id_string, 3) != 0) &&
            (strncmp(EPR_ENVISAT_PRODUCT_SAR,  product_id->id_string, 3) != 0) &&
            (strncmp(EPR_ENVISAT_PRODUCT_AATSR, product_id->id_string, 3) != 0)) {
        epr_set_err(e_err_invalid_product_id,
                    "epr_open_product: invalid product identifier");
        epr_close_product(product_id);
        return NULL;
    }

    if (product_id->id_string[9] != 'P') {
        char* ch = product_id->id_string + 9;
        if (*ch == 'C') {
            epr_log(e_log_info, "child product detected, mapping to 'P'");
        } else {
            sprintf(message_buffer, "unknown product sub-type '%c', mapping to 'P'", *ch);
            epr_log(e_log_warning, message_buffer);
        }
        *ch = 'P';
    }

    /* Set file to end of file in order to determine file size */
    if (fseek(product_id->istream, 0, SEEK_END) != 0) {
        epr_set_err(e_err_file_access_denied,
                    "epr_open_product: file seek failed");
        epr_close_product(product_id);
        return NULL;
    }

    /* Get file size */
    product_id->tot_size = (uint) ftell(product_id->istream);
    if (product_id->tot_size == (uint) -1) {
        epr_set_err(e_err_file_access_denied,
                    "epr_open_product: failed to determine file size");
        epr_close_product(product_id);
        return NULL;
    }
    sprintf(message_buffer, "product size: %u", product_id->tot_size);
    epr_log(e_log_debug, message_buffer);

    /* Set file pointer back to start */
    if (fseek(product_id->istream, 0, SEEK_SET) != 0) {
        epr_set_err(e_err_file_access_denied,
                    "epr_open_product: file seek failed");
        epr_close_product(product_id);
        return NULL;
    }

    product_id->record_info_cache = epr_create_ptr_array(32);
    product_id->param_table = epr_create_param_table();

    epr_log(e_log_info, "reading MPH");
    product_id->mph_record = epr_read_mph(product_id);

    epr_log(e_log_info, "reading SPH");
    product_id->sph_record = epr_read_sph(product_id);
    s_par = epr_set_dyn_dddb_params(product_id);

    epr_log(e_log_info, "reading all DSDs");
    product_id->dsd_array = epr_read_all_dsds(product_id);
    compare_ok = epr_compare_param(product_id);
    if (compare_ok == 0) {
        epr_set_err(e_err_invalid_value,
                    "epr_open_product: MPH_SIZE+SPH_SIZE must be equal to DSD[0].DS_OFFSET");
        epr_close_product(product_id);
        return NULL;
    }

    if (strncmp(EPR_ENVISAT_PRODUCT_MERIS, product_id->id_string, 3) == 0) {
        product_id->meris_iodd_version = epr_detect_meris_iodd_version(product_id);
        sprintf(message_buffer, "product type %s (MERIS IODD%d)", product_id->id_string, product_id->meris_iodd_version);
        epr_log(e_log_info, message_buffer);
    }

    epr_log(e_log_info, "creating dataset identifiers");
    product_id->dataset_ids = epr_create_dataset_ids(product_id);
    if (product_id->dataset_ids == NULL) {
        epr_close_product(product_id);
        return NULL;
    }

    epr_log(e_log_info, "creating band identifiers");
    product_id->band_ids = epr_create_band_ids(product_id);
    if (product_id->band_ids == NULL) {
        epr_close_product(product_id);
        return NULL;
    }

    /* Get scene size */
    product_id->scene_width = epr_compute_scene_width(product_id);
    product_id->scene_height = epr_compute_scene_height(product_id);
    sprintf(message_buffer, "product scene raster size: %u x %u", product_id->scene_width, product_id->scene_height);
    epr_log(e_log_debug, message_buffer);

    return product_id;
}



/*
   Function:    epr_close_product
   Access:      public API
   Changelog:   2002/01/05  nf initial version
 */
/**
 * Closes the ENVISAT product file determined by the gibven file identifier.
 */
int epr_close_product(EPR_SProductId* product_id) {
    /* Nothing to close, return */
    if (product_id == NULL) {
        return e_err_none;
    }

    epr_clear_err();
    if (!epr_check_api_init_flag()) {
        return epr_get_last_err_code();
    }

    assert(product_id->istream != NULL);
    if (fclose(product_id->istream) != 0) {
        epr_set_err(e_err_file_close_failed,
                    "epr_close_product: product file close failed");
        return epr_get_last_err_code();
    }
    product_id->istream = NULL;

    epr_log(e_log_info, "product closed: file path: ");
    epr_log(e_log_info, product_id->file_path);

    epr_free_product_id(product_id);

    return e_err_none;
}


/*
   Function:    epr_free_product_id
   Access:      private API implementation helper
   Changelog:   2002/01/05  nf initial version
 */
/**
 * Destructor for an <code>EPR_SProductId</code> structure instance.
 *
 * @param product_id the product file identifier to be destructed
 */
void epr_free_product_id(EPR_SProductId* product_id) {
    if (product_id == NULL)
        return;

    product_id->istream = NULL;

    epr_free_string(product_id->file_path);
    product_id->file_path = NULL;

    product_id->id_string[0] = '\0';

    epr_free_param_table(product_id->param_table);
    product_id->param_table = NULL;

    if (product_id->record_info_cache != NULL) {
        EPR_SRecordInfo* record_info = NULL;
        uint record_info_index = 0;

        for (record_info_index = 0; record_info_index < product_id->record_info_cache->length; record_info_index++) {
            record_info = (EPR_SRecordInfo*)epr_get_ptr_array_elem_at(product_id->record_info_cache, record_info_index);
            epr_free_record_info(record_info);
        }

        epr_free_ptr_array(product_id->record_info_cache);
        product_id->record_info_cache = NULL;
    }

    if (product_id->dsd_array != NULL) {
        EPR_SDSD* dsd = NULL;
        uint dsd_index = 0;

        for (dsd_index = 0; dsd_index < product_id->dsd_array->length; dsd_index++) {
            dsd = (EPR_SDSD*)epr_get_ptr_array_elem_at(product_id->dsd_array, dsd_index);
            epr_free_dsd(dsd);
        }

        epr_free_ptr_array(product_id->dsd_array);
        product_id->dsd_array = NULL;
    }

    if (product_id->mph_record != NULL) {
        epr_free_record(product_id->mph_record);
        product_id->mph_record = NULL;
    }

    if (product_id->sph_record != NULL) {
        epr_free_record(product_id->sph_record);
        product_id->sph_record = NULL;
    }

    if (product_id->dataset_ids != NULL) {
        EPR_SDatasetId* dataset_id = NULL;
        uint d_index = 0;

        for (d_index = 0; d_index < product_id->dataset_ids->length; d_index++) {
            dataset_id = (EPR_SDatasetId*)epr_get_ptr_array_elem_at(product_id->dataset_ids, d_index);
            epr_free_dataset_id(dataset_id);
        }

        epr_free_ptr_array(product_id->dataset_ids);
        product_id->dataset_ids = NULL;
    }

    if (product_id->band_ids != NULL) {
        EPR_SBandId* band_id = NULL;
        uint b_index = 0;

        for (b_index = 0; b_index < product_id->band_ids->length; b_index++) {
            band_id = (EPR_SBandId*)epr_get_ptr_array_elem_at(product_id->band_ids, b_index);
            epr_free_band_id(band_id);
        }

        epr_free_ptr_array(product_id->band_ids);
        product_id->band_ids = NULL;
    }

    product_id->tot_size = 0;

    free(product_id);
}


/**
 * Gets the scene width in pixel.
 *
 * @param product_id the product identifier, must not be <code>NULL</code>
 * @return width pixel number, or <code>0</code> if an error occured.
 */
uint epr_get_scene_width(const EPR_SProductId* product_id) {
    epr_clear_err();

    if (product_id == NULL) {
        epr_set_err(e_err_null_pointer,
                    "epr_get_scene_width: product_id must not be NULL");
        return (uint)0;
    }
    return product_id->scene_width;
}

/**
 * Gets the scene height in pixel.
 *
 * @param product_id the product identifier, must not be <code>NULL</code>
 * @return height pixel number, or <code>0</code> if an error occured.
 */
uint epr_get_scene_height(const EPR_SProductId* product_id) {
    epr_clear_err();

    if (product_id == NULL) {
        epr_set_err(e_err_null_pointer,
                    "epr_get_scene_height: product_id must not be NULL");
        return (uint)0;
    }
    return product_id->scene_height;
}


/*********************************** RECORD ***********************************/

EPR_SRecord* epr_get_sph(const EPR_SProductId* product_id) {
    epr_clear_err();

    if (product_id == NULL) {
        epr_set_err(e_err_null_pointer,
                    "epr_get_sph: product-identifier must not be NULL");
        return NULL;
    }
    return product_id->sph_record;
}

EPR_SRecord* epr_get_mph(const EPR_SProductId* product_id) {
    epr_clear_err();

    if (product_id == NULL) {
        epr_set_err(e_err_null_pointer,
                    "epr_get_mph: product-identifier must not be NULL");
        return NULL;
    }
    return product_id->mph_record;
}


/**
 * Gets the scene width in pixel.
 *
 * @param product_id the product identifier, must not be <code>NULL</code>
 *
 * @return width pixel number, or <code>0</code> if an error occured.
 */
uint epr_compute_scene_width(const EPR_SProductId* product_id) {
    EPR_SRecord* sph_record = NULL;
    uint scan_line_length;

    if (product_id == NULL) {
        epr_set_err(e_err_null_pointer,
                    "epr_compute_scene_width: product ID must not be NULL");
        return (uint)0;
    }

    sph_record = product_id->sph_record;

    if (strncmp(EPR_ENVISAT_PRODUCT_MERIS, product_id->id_string, 3) == 0) {
        const EPR_SField* field = field = epr_get_field(sph_record, "LINE_LENGTH");
        scan_line_length = epr_get_field_elem_as_uint(field, 0);
    } else if (strncmp(EPR_ENVISAT_PRODUCT_AATSR, product_id->id_string, 3) == 0) {
        scan_line_length = EPR_ATS_LINE_LENGTH;
    } else if (strncmp(EPR_ENVISAT_PRODUCT_ASAR, product_id->id_string, 3) == 0) {
        const EPR_SField* field = field = epr_get_field(sph_record, "LINE_LENGTH");
        scan_line_length = epr_get_field_elem_as_uint(field, 0);
    } else if (strncmp(EPR_ENVISAT_PRODUCT_SAR, product_id->id_string, 3) == 0) {
        const EPR_SField* field = field = epr_get_field(sph_record, "LINE_LENGTH");
        scan_line_length = epr_get_field_elem_as_uint(field, 0);
    } else {
        epr_set_err(e_err_illegal_arg,
                    "epr_compute_scene_width: unknown product type");
        scan_line_length = (uint)0;
    }

    return scan_line_length;
}

/**
 * Computes the scene height in pixel of a product. The scene height is
 * the minimum number of records in all measurement datasets.
 *
 * @param product_id the product identifier, must not be <code>NULL</code>
 * @return height pixel number, or <code>0</code> if an error occured.
 */
uint epr_compute_scene_height(const EPR_SProductId* product_id) {
    EPR_SDSD* dsd = NULL;
    uint min_num_mds_recs = 0;
    uint dsd_index;

    if (product_id == NULL) {
        epr_set_err(e_err_null_pointer,
                    "epr_compute_scene_height: product ID must not be NULL");
        return (uint)0;
    }

    for (dsd_index = 0; dsd_index < product_id->dsd_array->length; dsd_index++) {
        dsd = (EPR_SDSD*)epr_get_ptr_array_elem_at(product_id->dsd_array, dsd_index);
        if (epr_equal_names(dsd->ds_type, "M")) {
            if (dsd->num_dsr > min_num_mds_recs) {
                min_num_mds_recs = dsd->num_dsr;
            }
        }
    }

    if (min_num_mds_recs == 0) {
        epr_set_err(e_err_invalid_data_format,
                    "epr_compute_scene_height: product height was zero");
    }

    return min_num_mds_recs;
}



