/*
 * $Id: epr_band.c,v 1.2 2009-03-27 10:25:54 sabine Exp $
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
#include <math.h>

#include "epr_api.h"
#include "epr_core.h"
#include "epr_string.h"
#include "epr_ptrarray.h"
#include "epr_swap.h"
#include "epr_field.h"
#include "epr_record.h"
#include "epr_param.h"
#include "epr_dsd.h"
#include "epr_msph.h"
#include "epr_band.h"
#include "epr_bitmask.h"

#include "epr_dddb.h"

/**
 * Obtains all bands infos from the dddb.
 */
EPR_SPtrArray* epr_create_band_ids(EPR_SProductId* product_id) {
    EPR_SBandId* band_id = NULL;
    EPR_SPtrArray* band_ids = NULL;
    char test_block[1024];
    int bt_index;
    int i;
    const struct BandDescriptorTable* b_tables;
    int num_descr;

    if (product_id == NULL) {
        epr_set_err(e_err_null_pointer,
                    "epr_create_band_ids: product_id must not be NULL");
        return NULL;
    }


    /* @DDDB */

    b_tables = dddb_band_tables;
    bt_index = -1;
    for (i = 0; i < EPR_NUM_BAND_TABLES; i++) {
        const char* id = b_tables[i].name;
        if (strncmp(product_id->id_string, id, 10) == 0) {
            if (product_id->meris_iodd_version == 5) {
                if (strcmp(id, "MER_RR__1P_IODD5") == 0 ||
                        strcmp(id, "MER_FR__1P_IODD5") == 0) {
                    bt_index = i;
                }
            } else if (product_id->meris_iodd_version == 6) {
                if (strcmp(id, "MER_RR__2P_IODD6") == 0 ||
                        strcmp(id, "MER_FR__2P_IODD6") == 0) {
                    bt_index = i;
                }
            } else {
                bt_index = i;
            }
        }
        if (bt_index != -1) {
            break;
        }
    }
    if (bt_index == -1) {
        epr_set_err(e_err_null_pointer,
                    "epr_create_band_ids: unknown product type");
        return NULL;
    }

    band_ids = epr_create_ptr_array(16);
    num_descr = b_tables[bt_index].num_descriptors;
    for (i = 0; i < num_descr; i++) {

        band_id = (EPR_SBandId*) calloc(1, sizeof (EPR_SBandId));
        if (band_id == NULL) {
            epr_set_err(e_err_out_of_memory,
                        "epr_create_band_ids: out of memory");
            return NULL;
        }
        band_id->magic = EPR_MAGIC_BAND_ID;
        band_id->product_id = product_id;

        /* 1: band_name */
        epr_assign_string(&band_id->band_name, b_tables[bt_index].descriptors[i].id);
        /* 2: dataset_name */
        band_id->dataset_ref = epr_get_ref_struct(product_id, b_tables[bt_index].descriptors[i].rec_name);
        if (band_id->dataset_ref.dataset_id == NULL) {
            epr_set_err(e_err_invalid_dataset_name,
                        "epr_create_band_ids: invalid dataset name in DDDB");
            epr_free_band_id(band_id);
            return NULL;
        }
        /* 3: sample_offset */
        band_id->sample_model = b_tables[bt_index].descriptors[i].sample_offset;
        /* 4: band_datatype */
        band_id->data_type = b_tables[bt_index].descriptors[i].type;
        /* 5: spectr_band_index*/
        band_id->spectr_band_index = b_tables[bt_index].descriptors[i].spectral_index;
        /* 6: scaling_method*/
        if (b_tables[bt_index].descriptors[i].scale_method == e_smid_non) {
            band_id->scaling_method = 0;
            band_id->scaling_offset = 0.0;
            band_id->scaling_factor = 1.0;
        } else {
            band_id->scaling_method = b_tables[bt_index].descriptors[i].scale_method;
            /* 7: scaling_offset*/
            strcpy (test_block, b_tables[bt_index].descriptors[i].scale_offset);
            if (test_block == NULL) {
                band_id->scaling_offset = 0.0;
            } else {
                float scaling_offset = (float)atof(test_block);
                if (epr_numeral_suspicion(test_block) == 1) {
                    band_id->scaling_offset = scaling_offset;
                } else {
                    scaling_offset = epr_get_scaling_params(product_id, test_block);
                    if (scaling_offset == -909.909) { /* @todo what an ugly return value. Eeeek!*/
                        epr_set_err(e_err_invalid_dataset_name,
                                    "epr_create_band_ids: invalid dataset name in dddb");
                        epr_free_band_id(band_id);
                        return NULL;
                    }
                    band_id->scaling_offset = scaling_offset;
                }
            }
            /* 8: scaling_factor*/
            strcpy (test_block, b_tables[bt_index].descriptors[i].scale_factor);
            if (test_block == NULL) {
                band_id->scaling_factor = 0.0;
            } else {
                float scaling_factor = (float)atof(test_block);
                if (epr_numeral_suspicion(test_block) == 1) {
                    band_id->scaling_factor = scaling_factor;
                } else {
                    scaling_factor = epr_get_scaling_params(product_id, test_block);
                    if (scaling_factor == -909.909) { /* @todo what an ugly return value. Eeeek!*/
                        epr_set_err(e_err_invalid_dataset_name,
                                    "epr_create_band_ids: invalid dataset name in dddb");
                        epr_free_band_id(band_id);
                        return NULL;
                    }
                    band_id->scaling_factor = scaling_factor;
                }
            }
        }
        /* 9: bit_expr*/
        epr_assign_string(&band_id->bm_expr, b_tables[bt_index].descriptors[i].bitmask_expr);
        /* 10: flags_definition_file*/
        if (b_tables[bt_index].descriptors[i].flag_coding_name != NULL) {
            band_id->flag_coding = epr_create_flag_coding(product_id, b_tables[bt_index].descriptors[i].flag_coding_name);
            if (band_id->flag_coding == NULL) {
                epr_set_err(e_err_out_of_memory,
                            "epr_create_band_ids: out of memory");
                epr_free_band_id(band_id);
                return NULL;
            }
        } else {
            band_id->flag_coding = NULL;
        }
        /* 11: unit*/
        epr_assign_string(&band_id->unit, b_tables[bt_index].descriptors[i].unit);
        /* 12: description*/
        epr_assign_string(&band_id->description, b_tables[bt_index].descriptors[i].description);

        /* lines_flipped*/
        if (strncmp(product_id->id_string, EPR_ENVISAT_PRODUCT_MERIS, 3) == 0
                || strncmp(product_id->id_string, EPR_ENVISAT_PRODUCT_AATSR, 3) == 0) {
            band_id->lines_mirrored = TRUE;
        } else {
            if (strncmp(product_id->id_string, EPR_ENVISAT_PRODUCT_ASAR, 3) == 0
                    && strncmp(product_id->id_string, "ASA_IMG", 7) != 0
                    && strncmp(product_id->id_string, "ASA_APG", 7) != 0) {
                band_id->lines_mirrored = TRUE;
            } else {
                if (strncmp(product_id->id_string, EPR_ENVISAT_PRODUCT_SAR, 3) == 0
                        && strncmp(product_id->id_string, "SAR_IMG", 7) != 0
                        && strncmp(product_id->id_string, "SAR_APG", 7) != 0) {
                    band_id->lines_mirrored = TRUE;
                } else {
                    band_id->lines_mirrored = FALSE;
                }
            }
        }

        epr_add_ptr_array_elem(band_ids, band_id);
    }

    return band_ids;
}


uint epr_get_num_bands(EPR_SProductId* product_id) {
    epr_clear_err();
    if (!epr_check_api_init_flag()) {
        return 0;
    }

    if (product_id == NULL) {
        epr_set_err(e_err_null_pointer,
                    "epr_get_num_bands: product_id must not be NULL");
        return (uint) -1;
    }
    return product_id->band_ids->length;
}

EPR_SBandId* epr_get_band_id_at(EPR_SProductId* product_id, uint index) {
    EPR_SBandId* band_id = NULL;

    epr_clear_err();

    if (product_id == NULL) {
        epr_set_err(e_err_null_pointer,
                    "epr_get_band_id_at: product_id must not be NULL");
        return NULL;
    }
    if (index >= product_id->band_ids->length) {
        epr_set_err(e_err_index_out_of_range,
                    "epr_get_band_id_at: band index out of range");
        return NULL;
    }

    band_id = (EPR_SBandId*)epr_get_ptr_array_elem_at(product_id->band_ids, index);
    return band_id;
}

EPR_SBandId* epr_get_band_id(EPR_SProductId* product_id, const char* band_name) {
    EPR_SBandId* band_id = NULL;
    int num_bands, i;

    epr_clear_err();

    if (product_id == NULL) {
        epr_set_err(e_err_null_pointer,
                    "epr_get_band_id: product_id must not be NULL");
        return NULL;
    }
    if (band_name == NULL) {
        epr_set_err(e_err_null_pointer,
                    "epr_get_band_id: dataset_name must not be NULL");
        return NULL;
    }

    num_bands = epr_get_num_bands(product_id);
    for (i = 0; i < num_bands; i++) {
        band_id = epr_get_band_id_at(product_id, i);
        if (epr_equal_names(band_name, epr_get_band_name(band_id))) {
            return band_id;
        }
    }
    epr_set_err(e_err_invalid_band_name,
                "epr_get_band_id: band not found");
    return NULL;
}

const char* epr_get_band_name(EPR_SBandId* band_id) {
    epr_clear_err();

    if (band_id == NULL) {
        epr_set_err(e_err_null_pointer,
                    "epr_get_band_name: band_id must not be NULL");
        return NULL;
    }
    return band_id->band_name;
}

/**
 * Release the memory allocated through a band ID.
 *
 * @param band the dataset description identifier, if <code>NULL</code> the function
 *        immediately returns zero.
 * @return zero for success, an error code otherwise
 */
void epr_free_band_id(EPR_SBandId* band_id) {
    if (band_id == NULL)
        return;

    band_id->dataset_ref.elem_index = -1;
    band_id->dataset_ref.field_index = -1;
    band_id->dataset_ref.dataset_id = NULL;

    epr_free_and_null_string(&band_id->band_name);
    epr_free_and_null_string(&band_id->bm_expr);

    epr_free_flag_coding(band_id->flag_coding);
    band_id->flag_coding = NULL;

    band_id->spectr_band_index = 0;
    band_id->scaling_offset  = 0;
    band_id->scaling_factor  = 0;
    band_id->data_type = e_tid_unknown;

    epr_free_and_null_string(&band_id->unit);
    epr_free_and_null_string(&band_id->description);

    band_id->lines_mirrored = FALSE;

    free(band_id);
}


/**
 * Gets the scaling params: factor or offset by the given dataset_id, field_index, elem_index
 *
 * @param product_id the the product file identifier
 * @param str the string with the name, separator ('.') and indexes.
 * @return the dataset_id, field_index and elem_index (-1 if no).
 *    <code>NULL</code> if correspondent dataset name was not found.
 */
float epr_get_scaling_params(EPR_SProductId* product_id,  const char* str) {
    EPR_SDatasetRef scal_fact;
    const EPR_SField* field = NULL;
    EPR_SRecord* record = NULL;
    float ziff;

    scal_fact = epr_get_ref_struct(product_id, str);
    if (scal_fact.dataset_id == NULL) {
        return (float)(-909.909);
    }

    /*'Scaling_Factor_GADS'*/
    record = epr_create_record(scal_fact.dataset_id);
    record = epr_read_record(scal_fact.dataset_id, 0, record);

    field = epr_get_field_at(record, scal_fact.field_index - 1);
    ziff = epr_get_field_elem_as_float(field, (uint)(scal_fact.elem_index - 1));

    epr_free_record(record);

    return ziff;
}


/**
 * Gets the scaling factor by the given dataset_id, field_index, elem_index
 *
 * @param product_id the the product file identifier
 * @param str the string with the name, separator ('.') and indexes.
 * @return the dataset_id, field_index and elem_index (-1 if no).
 *    <code>NULL</code> if correspondent dataset name was not found.
 */
float epr_get_scaling_factor(EPR_SProductId* product_id,  const char* str) {
    EPR_SDatasetRef scal_fact;
    const EPR_SField* field = NULL;
    EPR_SRecord* record = NULL;
    float ziff;

    scal_fact = epr_get_ref_struct(product_id, str);
    if (scal_fact.dataset_id == NULL) {
        return (float)(-909.909);
    }

    /*'Scaling_Factor_GADS'*/
    record = epr_create_record(scal_fact.dataset_id);
    record = epr_read_record(scal_fact.dataset_id, 0, record);

    field = epr_get_field_at(record, scal_fact.field_index - 1);
    ziff = epr_get_field_elem_as_float(field, (uint)(scal_fact.elem_index - 1));

    epr_free_record(record);

    return ziff;
}

/**
 * Gets the dataset_id, field_index and elem_index
 *
 * @param product_id the the product file identifier
 * @param str the string with the name, separator ('.') and indexes.
 * @return the dataset_id, field_index and elem_index (-1 if no).
 *    <code>NULL</code> if correspondent dataset name was not found.
 */
EPR_SDatasetRef epr_get_ref_struct(EPR_SProductId* product_id, const char* str) {
    EPR_SDatasetRef ref_struct;
    int pos = 0;
    char* stopstring;
    char* token;

    ref_struct.dataset_id = NULL;
    ref_struct.field_index = -1;
    ref_struct.elem_index = -1;

    token = epr_str_tok(str, ".", &pos);

    ref_struct.dataset_id = epr_get_dataset_id(product_id, token);
    if (ref_struct.dataset_id == NULL) {
        epr_free_and_null_string(&token);
        return ref_struct;
    }
    epr_free_and_null_string(&token);

    token = epr_str_tok(str, ".", &pos);
    if (token == NULL) {
        ref_struct.field_index = -1;
    } else {
        ref_struct.field_index = strtol(token, &stopstring, 10);
    }
    epr_free_and_null_string(&token);

    token = epr_str_tok(str, ".", &pos);
    if (token == NULL) {
        ref_struct.elem_index = -1;
    } else {
        ref_struct.elem_index = strtol(token, &stopstring, 10);
    }
    epr_free_and_null_string(&token);

    return ref_struct;
}

/**
 * Converts the given string into a scaling method identifier.
 *
 * @param str the string to be converted.
 * @return the scaling method identifier represented by the given string.
 *         If the string is equal of '*' the value
 *         <code>e_non_smid</code> is returned.
 */
EPR_EScalingMethod epr_str_to_scaling_method(const char* str) {
    assert(str != NULL);
    if (epr_equal_names(str, "Linear_Scale"))
        return e_smid_lin;
    else if (epr_equal_names(str, "Log_Scale"))
        return e_smid_log;
    else
        return e_smid_non;
}


/**
 * Converts the given string into a sample offset identifier.
 *
 * @param str the string to be converted.
 * @return the sample offset identifier represented by the given string.
 *         If the string is equal of '*' the value
 *         <code>e_none_samoff</code> is returned.
 */
EPR_ESampleModel epr_str_to_sample_offset(const char* str) {
    assert(str != NULL);
    if (epr_equal_names(str, "1OF2"))
        return e_smod_1OF2;
    else if (epr_equal_names(str, "2OF2"))
        return e_smod_2OF2;
    else if (epr_equal_names(str, "3TOI"))
        return e_smod_3TOI;
    else if (epr_equal_names(str, "2TOF"))
        return e_smod_2TOF;
    else
        return e_smod_1OF1;
}

/**
 * Creates a raster to be used for reading bitmasks. The raster returned is always of type <code>byte</code>.
 *
 * @param source_width the width (across track dimension) of the source to be read into the raster. See description of epr_create_compatible_raster.
 * @param source_height the height (along track dimension) of the source to be read into the raster. See description of epr_create_compatible_raster.
 * @param source_step_x the subsampling step across track of the source when reading into the raster. See description of epr_create_compatible_raster.
 * @param source_step_y the subsampling step along track of the source when reading into the raster. See description of epr_create_compatible_raster.
 * @return the new raster instance
 *         or <code>NULL</code> if an error occured.
 */
EPR_SRaster* epr_create_bitmask_raster(uint source_width,
                                       uint source_height,
                                       uint source_step_x,
                                       uint source_step_y) {
    return epr_create_raster(e_tid_uchar,
                             source_width,
                             source_height,
                             source_step_x,
                             source_step_y);
}


/**
 * Creates a raster for the given datatype and dimension.
 *
 * @param data_type the data type identifier
 * @param source_width the source's width
 * @param source_height the source's height
 * @param source_step_x the sub-sampling in X
 * @param source_step_y the sub-sampling in Y
 * @return the new raster instance
 *         or <code>NULL</code> if an error occured.
 */
EPR_SRaster* epr_create_raster(EPR_EDataTypeId data_type,
                               uint source_width,
                               uint source_height,
                               uint source_step_x,
                               uint source_step_y) {
    EPR_SRaster* raster = NULL;
    uint num_elems;

    epr_clear_err();

    if (data_type == e_tid_string ||
            data_type == e_tid_spare ||
            data_type == e_tid_time) {
        epr_set_err(e_err_illegal_data_type, "epr_create_raster: illegal data type");
        return NULL;
    }

    raster = (EPR_SRaster*) calloc(1, sizeof (EPR_SRaster));
    if (raster == NULL) {
        epr_set_err(e_err_out_of_memory, "epr_create_raster: out of memory");
        return NULL;
    }


    raster->magic         = EPR_MAGIC_RASTER;
    raster->data_type     = data_type;
    raster->elem_size     = epr_get_data_type_size(data_type);
    raster->source_height = source_height;
    raster->source_width  = source_width;
    raster->source_step_x = source_step_x;
    raster->source_step_y = source_step_y;
    raster->raster_width  = (source_width  - 1) / source_step_x + 1;
    raster->raster_height = (source_height - 1) / source_step_y + 1;

    num_elems = raster->raster_width * raster->raster_height;

    raster->buffer = calloc(raster->elem_size, num_elems);
    if (raster->buffer == NULL) {
        epr_free_raster(raster);
        epr_set_err(e_err_out_of_memory, "epr_create_raster: out of memory");
        return NULL;
    }

    return raster;
}


/**
 * Creates a compatible raster for the given band.
 *
 * @param band_id the band identifier, must not be <code>NULL</code>
 * @return the new raster instance
 *         or <code>NULL</code> if an error occured.
 */
EPR_SRaster* epr_create_compatible_raster(EPR_SBandId* band_id,
        uint source_width,
        uint source_height,
        uint source_step_x,
        uint source_step_y) {
    epr_clear_err();

    if (band_id == NULL) {
        epr_set_err(e_err_invalid_band, "epr_create_raster: band_id must not be NULL");
        return NULL;
    }
    return epr_create_raster(band_id->data_type,
                             source_width,
                             source_height,
                             source_step_x,
                             source_step_y);
}


void epr_free_raster(EPR_SRaster* raster) {
    epr_clear_err();

    if (raster == NULL)
        return;

    raster->data_type = e_tid_unknown;
    raster->elem_size     = 0;
    raster->raster_height = 0;
    raster->raster_width  = 0;
    raster->source_height = 0;
    raster->source_width  = 0;
    raster->source_step_x = 0;
    raster->source_step_y = 0;

    if (raster->buffer != NULL) {
        free(raster->buffer);
        raster->buffer = NULL;
    }

    free(raster);
}



int epr_read_band_raster(EPR_SBandId* band_id,
                         int offset_x,
                         int offset_y,
                         EPR_SRaster* raster/*, EPR_SRaster** bitmask_raster*/) {

    EPR_SProductId* product_id = NULL;
    EPR_SDatasetId* dataset_id = NULL;
    char* rec_type;

    epr_clear_err();

    if (band_id == NULL) {
        epr_set_err(e_err_invalid_band,
                    "epr_read_band_raster: band_id must not be NULL");
        return epr_get_last_err_code();
    }
    if (band_id->data_type != raster->data_type) {
        epr_set_err(e_err_illegal_data_type,
                    "epr_read_band_raster: illegal raster data type");
        return epr_get_last_err_code();
    }
    if (raster->buffer == NULL) {
        epr_set_err(e_err_illegal_arg,
                    "epr_read_band_raster: raster->buffer must not be NULL");
        return epr_get_last_err_code();
    }
    if ((offset_x<0) || (offset_y<0) || (raster->raster_width<0) || (raster->raster_height<0) || (raster->source_step_x<0) || (raster->source_step_y<0)) {
        epr_set_err(e_err_invalid_value,
                    "epr_read_band_raster: all digit parameter must be positive");
        return epr_get_last_err_code();
    }
    /*  removed because the source_step_x can truly be greater than raster_width.
    if ((raster->source_step_x>raster->raster_width) || (raster->source_step_y>raster->raster_height)) {
        epr_set_err(e_err_invalid_value,
                "epr_read_band_raster: too small raster sizes or large steps");
        return epr_get_last_err_code();
    }
    */
    product_id = band_id->product_id;
    dataset_id = band_id->dataset_ref.dataset_id;
    rec_type = dataset_id->dsd->ds_type;
    if (strcmp(rec_type, "M") == 0) {
        if (epr_read_band_measurement_data(band_id,
                                           offset_x,
                                           offset_y,
                                           raster) != 0) {
            epr_set_err(e_err_file_read_error,
                        "epr_read_band_raster: unsuccessfully reading band measurement data");
            return epr_get_last_err_code();
        }
        if (band_id->bm_expr != NULL) {
            EPR_SRaster* bm_raster;
            int rd_bm;

            bm_raster = epr_create_raster(e_tid_uchar, /*was char*/
                                          raster->source_width,
                                          raster->source_height,
                                          raster->source_step_x,
                                          raster->source_step_y);


            rd_bm = epr_read_bitmask_raster(product_id,
                                            band_id->bm_expr,
                                            offset_x,
                                            offset_y,
                                            bm_raster);

            epr_zero_invalid_pixels(raster, bm_raster);

            epr_free_raster(bm_raster);

        }
    } else if (strcmp(rec_type, "A") == 0) {
        if (epr_read_band_annotation_data
                (band_id, offset_x, offset_y, raster) == 1) {
            epr_set_err(e_err_file_read_error,
                        "epr_read_band_raster: unsuccessfully reading band annotation data");
            return epr_get_last_err_code();
        }
    } else {
        epr_set_err(e_err_invalid_value,
                    "epr_read_band_raster: illegat DS-TYPE; 'A' or'M' will be accepted");
        return epr_get_last_err_code();
    }
    return e_err_none;
}

/**
 * Reads the measurement data and converts its into physical values.
 *
 * @param band_id the information about properties and quantities of ENVISAT data.
 * @param offset_x X-coordinate in pixel co-ordinates (zero-based) of the upper right corner raster to search
 * @param offset_y Y-coordinate in pixel co-ordinates (zero-based) of the upper right corner raster to search
 * @param raster the instance to the buffer information was used
 *
 * @return zero for success, an error code otherwise
 */
int epr_read_band_measurement_data(EPR_SBandId* band_id,
                                   int offset_x,
                                   int offset_y,
                                   EPR_SRaster* raster) {
    EPR_SProductId* product_id = NULL;
    const EPR_SField* field = NULL;
    EPR_SFieldInfo* field_info = NULL;
    EPR_SDatasetId* dataset_id = NULL;
    EPR_SRecord* record = NULL;
    EPR_SRecord* sph_record = NULL;
    EPR_EDataTypeId band_datatype, datatype_id;
    EPR_ESampleModel band_smod;
    uint rec_size;
    uint rec_numb;
    int iY, raster_pos, delta_raster_pos;
    int offset_x_mirrored = 0;
    uint scan_line_length;
    EPR_FLineDecoder decode_func;
    uint scene_width;

    product_id = band_id->product_id;

    if (strncmp(EPR_ENVISAT_PRODUCT_MERIS, product_id->id_string, 3) == 0) {
        sph_record = product_id->sph_record;
        field = epr_get_field(sph_record, "LINE_LENGTH");
        scan_line_length = epr_get_field_elem_as_uint(field, 0);
    } else if (strncmp(EPR_ENVISAT_PRODUCT_AATSR, product_id->id_string, 3) == 0) {
        scan_line_length = EPR_ATS_LINE_LENGTH;
    } else if (strncmp(EPR_ENVISAT_PRODUCT_ASAR, product_id->id_string, 3) == 0) {
        scan_line_length = epr_get_scene_width(product_id);
    } else if (strncmp(EPR_ENVISAT_PRODUCT_SAR, product_id->id_string, 3) == 0) {
        scan_line_length = epr_get_scene_width(product_id);
    } else {
        epr_set_err(e_err_illegal_arg,
                    "epr_read_band_measurement_data: scan line length unknown");
        return epr_get_last_err_code();
    }

    dataset_id = band_id->dataset_ref.dataset_id;
    /*the length of measurement record size*/
    rec_size = dataset_id->dsd->dsr_size;
    /*the number of measurement records*/
    rec_numb = dataset_id->dsd->num_dsr;
    /*data type in the band*/
    band_datatype = band_id->data_type;
    /*data model in the band*/
    band_smod = band_id->sample_model;
    record = epr_create_record(dataset_id);
    field_info = (EPR_SFieldInfo*)epr_get_ptr_array_elem_at(record->info->field_infos, band_id->dataset_ref.field_index - 1);
    datatype_id = field_info->data_type_id;

    /* if the user raster (or part of) is outside bbox in source coordinates*/
    if (offset_x + raster->raster_width > (int)scan_line_length) {
        epr_set_err(e_err_illegal_arg,
                    "epr_read_band_measurement_data: raster x co-ordinates out of bounds");
        epr_free_record(record);
        return epr_get_last_err_code();
    }
    if (offset_y + raster->raster_height > (int)(rec_numb)) {
        epr_set_err(e_err_illegal_arg,
                    "epr_read_band_measurement_data: raster y co-ordinates out of bounds");
        epr_free_record(record);
        return epr_get_last_err_code();
    }
    raster_pos = 0;
    delta_raster_pos = (int)floor((raster->source_width - 1) / raster->source_step_x) + 1;

    /*select the correspondent function to scaling and transform data type*/
    decode_func = select_line_decode_function(band_datatype, band_smod, datatype_id);
    if (decode_func == NULL) {
        epr_set_err(e_err_illegal_data_type,
                    "epr_read_band_measurement_data: internal error: unknown data type");
        epr_free_record(record);
        return epr_get_last_err_code();
    }

    scene_width = band_id->product_id->scene_width;
    if (band_id->lines_mirrored) {
        offset_x_mirrored = (scene_width - 1) - (offset_x + raster->source_width - 1);
    } else {
        offset_x_mirrored = offset_x;
    }

    for (iY = offset_y; (uint)iY < offset_y + raster->source_height; iY += raster->source_step_y ) {

        /*get the next record by the given name*/
        record = epr_read_record(dataset_id, iY, record);
        if (record == NULL) {
            return epr_get_last_err_code();
        }
        /*get the field at its number*/
        field = epr_get_field_at(record, band_id->dataset_ref.field_index - 1);
        /*get the scaled "line" of physical values*/
        decode_func(field->elems, band_id, offset_x_mirrored, raster->source_width, raster->source_step_x, raster->buffer, raster_pos);
        /*locate "data point" for the next "line"*/
        raster_pos += delta_raster_pos;
    }

    if (band_id->lines_mirrored) {
        if (band_datatype == e_tid_float) {
            mirror_float_array((float*)raster->buffer, raster->raster_width, raster->raster_height);
        } else if (band_datatype == e_tid_uchar || band_datatype == e_tid_char) {
            mirror_uchar_array((uchar*)raster->buffer, raster->raster_width, raster->raster_height);
        } else if (band_datatype == e_tid_ushort || band_datatype == e_tid_short) {
            mirror_ushort_array((ushort*)raster->buffer, raster->raster_width, raster->raster_height);
        } else if (band_datatype == e_tid_uint || band_datatype == e_tid_int) {
            mirror_uint_array((uint*)raster->buffer, raster->raster_width, raster->raster_height);
        } else {
            epr_set_err(e_err_illegal_data_type,
                        "epr_read_band_measurement_data: internal error: unknown data type");
            epr_free_record(record);
            return epr_get_last_err_code();
        }
    }

    epr_free_record(record);

    return 0;
}


/**
 * Reads the annotation data and converts its into physical values.
 *
 * @param band_id the information about properties and quantities of ENVISAT data.
 * @param offset_x X-coordinate in pixel co-ordinates (zero-based) of the upper right corner raster to search
 * @param offset_y Y-coordinate in pixel co-ordinates (zero-based) of the upper right corner raster to search
 * @param raster the instance to the buffer information was used
 *
 * @return zero for success, an error code otherwise
 */
int epr_read_band_annotation_data(EPR_SBandId* band_id,
                                  int offset_x,
                                  int offset_y,
                                  EPR_SRaster* raster) {
    EPR_SProductId* product_id = NULL;
    const EPR_SField* field = NULL;
    const EPR_SField* field_beg = NULL;
    const EPR_SField* field_end = NULL;
    EPR_SFieldInfo* field_info = NULL;
    EPR_SDatasetId* dataset_id = NULL;
    EPR_SRecord* record = NULL;
    EPR_SRecord* record_beg = NULL;
    EPR_SRecord* record_end = NULL;
    EPR_SRecord* sph_record = NULL;
    EPR_EDataTypeId band_datatype = 0, datatype_id = 0;
    EPR_ESampleModel band_smod = 0;
    uint rec_size = 0;
    uint rec_numb = 0;
    uint lines_per_tie_pt, samples_per_tie_pt, scan_line_length;
    int iY, raster_pos, delta_raster_pos;
    EPR_FArrayTransformer transform_array_func = NULL;
    int y_beg, y_end, y_beg_old, y_end_old;
    int offset_x_mirrored = 0;
    uint num_elems = 0;
    float y_mod = 0;
    float scan_offset_x = 0;
    float scan_offset_y = 0;
    void* line_beg_buffer = NULL;
    void* line_end_buffer = NULL;

    product_id = band_id->product_id;

    dataset_id = band_id->dataset_ref.dataset_id;
    /*the length of annotation record size*/
    rec_size = dataset_id->dsd->dsr_size;
    /*the number of annotation records*/
    rec_numb = dataset_id->dsd->num_dsr;
    /*data type in the band*/
    band_datatype = band_id->data_type;
    /*data model in the band*/
    band_smod = band_id->sample_model;
    record = epr_create_record(dataset_id);
    field_info = (EPR_SFieldInfo*)epr_get_ptr_array_elem_at(record->info->field_infos, band_id->dataset_ref.field_index - 1);
    datatype_id = field_info->data_type_id;


    /*find LINES_PER_TIE_PT & SAMPLES_PER_TIE_PT for different products*/
    if (strncmp(EPR_ENVISAT_PRODUCT_MERIS, product_id->id_string, 3) == 0) {
        /*elements number in the band (e.g.71)*/
        scan_offset_x = 0.0F; /*!! was 0.5F !!*/
        scan_offset_y = 0.0F; /*!! was 0.5F !!*/
        num_elems = field_info->num_elems;
        sph_record = product_id->sph_record;
        field = epr_get_field(sph_record, "LINES_PER_TIE_PT");
        lines_per_tie_pt = epr_get_field_elem_as_uint(field, 0);
        field = epr_get_field(sph_record, "SAMPLES_PER_TIE_PT");
        samples_per_tie_pt = epr_get_field_elem_as_uint(field, 0);
        field = epr_get_field(sph_record, "LINE_LENGTH");
        scan_line_length = epr_get_field_elem_as_uint(field, 0);
    } else if (strncmp(EPR_ENVISAT_PRODUCT_AATSR, product_id->id_string, 3) == 0) {
        scan_offset_y = 0.0F; /*!! EPR-7: was 0.5F !!*/
        scan_line_length = EPR_ATS_LINE_LENGTH;
        lines_per_tie_pt = EPR_AATSR_LINES_PER_TIE_PT;
        num_elems = field_info->num_elems;
        if (num_elems == EPR_ATS_NUM_PER_POINT_ACROSS_LOCAT) {
            scan_offset_x = -19.0F;
            samples_per_tie_pt = 25;
        } else if (num_elems == EPR_ATS_NUM_PER_POINT_ACROSS_SOLAR) {
            scan_offset_x = 6.0F;
            samples_per_tie_pt = 50;
        } else {
            epr_set_err(e_err_invalid_value, "epr_read_band_annotation_data: internal error: illegal value for samples_per_tie_pt");
            epr_free_record(record);
            return epr_get_last_err_code();
        }
    } else if ((strncmp(EPR_ENVISAT_PRODUCT_ASAR, product_id->id_string, 3) == 0) ||
               (strncmp(EPR_ENVISAT_PRODUCT_SAR, product_id->id_string, 3) == 0)) {
        EPR_SDatasetId* dataset_id = NULL;
        uint num_rec;
        scan_offset_x = 0.5F; /* @todo CHECK THIS FOR ASAR! */
        scan_offset_y = 0.5F;
        scan_line_length = epr_get_scene_width(product_id);
        samples_per_tie_pt = scan_line_length / (EPR_ASAR_NUM_PER_POINT_ACROSS_LOCAT - 1);
        dataset_id = epr_get_dataset_id(product_id, "GEOLOCATION_GRID_ADS");
        num_rec = epr_get_num_records(dataset_id);
        lines_per_tie_pt = epr_get_scene_height(product_id) / (num_rec - 1);
        num_elems = field_info->num_elems;
    } else {
        epr_set_err(e_err_illegal_arg,
                    "epr_read_band_annotation_data: unhandled ENVISAT product type");
        epr_free_record(record);
        return epr_get_last_err_code();
    }

    /*memory allocate for the increasingly begin tie point line*/
    line_beg_buffer = calloc(sizeof(float), num_elems);
    if (line_beg_buffer == NULL) {
        epr_set_err(e_err_out_of_memory, "epr_read_band_annotation_data: out of memory");
        epr_free_record(record);
        return epr_get_last_err_code();
    }
    /*memory allocate for the increasingly end tie point line*/
    line_end_buffer = calloc(sizeof(float), num_elems);
    if (line_end_buffer == NULL)  {
        epr_set_err(e_err_out_of_memory, "epr_read_band_annotation_data: out of memory");
        epr_free_record(record);
        free(line_beg_buffer);
        return epr_get_last_err_code();
    }
    /* if the user raster (or its part) is outside of orbit in source coordinates*/
    if (offset_x + raster->raster_width > (int)scan_line_length) {
        epr_set_err(e_err_illegal_arg,
                    "epr_read_band_data: raster x co-ordinates out of bounds");
        epr_free_record(record);
        free(line_beg_buffer);
        free(line_end_buffer);
        return epr_get_last_err_code();
    }
    if (offset_y + raster->raster_height > (int)(rec_numb * lines_per_tie_pt)) {
        epr_set_err(e_err_illegal_arg,
                    "epr_read_band_data: raster y co-ordinates out of bounds");
        epr_free_record(record);
        free(line_beg_buffer);
        free(line_end_buffer);
        return epr_get_last_err_code();
    }
    raster_pos = 0;

    delta_raster_pos = (int)floor((raster->source_width - 1) / raster->source_step_x) + 1;

    /*select the correspondent function to scaling and transform data type*/
    transform_array_func = select_transform_array_function(band_datatype, datatype_id);
    if (transform_array_func == NULL) {
        epr_set_err(e_err_illegal_data_type,
                    "epr_read_band_annotation_data: internal error: illegal data type");
        epr_free_record(record);
        free(line_beg_buffer);
        free(line_end_buffer);
        return epr_get_last_err_code();
    }
    y_beg_old = 9999;
    y_end_old = 9999;

    if (band_id->lines_mirrored) {
        offset_x_mirrored = num_elems - (offset_x + raster->source_width - 1) - 1;
    } else {
        offset_x_mirrored = offset_x;
    }

    for (iY = offset_y; (uint)iY < offset_y + raster->source_height; iY += raster->source_step_y ) {

        /*find the increasing neighbour begin and end tie point lines*/
        y_mod = ((float)iY - scan_offset_y) / lines_per_tie_pt;
        y_beg = (uint)floor(y_mod);

        if (y_beg < 0) {
            y_beg = 0;
        }
        if ((uint)y_beg > dataset_id->dsd->num_dsr - 2) {
            y_beg = dataset_id->dsd->num_dsr - 2;
        }

        y_mod -= y_beg;
        y_end = y_beg + 1;

        /*as long as between increasing neighbour tie point lines, not to change them*/
        if (y_beg_old != y_beg) {
            record_beg = epr_read_record(dataset_id, y_beg, record_beg);
            y_beg_old = y_beg;
        }
        if (y_end_old != y_end) {
            record_end = epr_read_record(dataset_id, y_end, record_end);
            y_end_old = y_end;
        }

        /*get the values for the increasing neighbour tie point lines*/
        field_beg = epr_get_field_at(record_beg, band_id->dataset_ref.field_index - 1);
        field_end = epr_get_field_at(record_end, band_id->dataset_ref.field_index - 1);

        /*transform and scale the values for the increasing neighbour tie point lines*/
        transform_array_func(field_beg->elems, band_id, line_beg_buffer, num_elems);
        transform_array_func(field_end->elems, band_id, line_end_buffer, num_elems);

        /*get the "line" of interpolated physical values from tie point data*/
        decode_tiepoint_band(line_beg_buffer, line_end_buffer,
                             samples_per_tie_pt, num_elems, band_id, offset_x, scan_offset_x, y_mod,
                             raster->source_width, raster->source_step_x, raster->buffer, raster_pos);
        /*locate "data point" for the next "line"*/
        raster_pos += delta_raster_pos;
    }

    if (strncmp(EPR_ENVISAT_PRODUCT_MERIS, product_id->id_string, 3) == 0) {
        mirror_float_array((float*)raster->buffer, raster->raster_width, raster->raster_height);
    } else {
        if (strncmp(EPR_ENVISAT_PRODUCT_AATSR, product_id->id_string, 3) == 0) {
            mirror_float_array((float*)raster->buffer, raster->raster_width, raster->raster_height);
        } else {
            if (strncmp(EPR_ENVISAT_PRODUCT_ASAR, product_id->id_string, 3) == 0
                    && strncmp(product_id->id_string, "ASA_IMG", 7) != 0
                    && strncmp(product_id->id_string, "ASA_APG", 7) != 0) {
                mirror_float_array((float*)raster->buffer, raster->raster_width, raster->raster_height);
            } else {
                if (strncmp(EPR_ENVISAT_PRODUCT_SAR, product_id->id_string, 3) == 0
                        && strncmp(product_id->id_string, "SAR_IMG", 7) != 0
                        && strncmp(product_id->id_string, "SAR_APG", 7) != 0) {
                    mirror_float_array((float*)raster->buffer, raster->raster_width, raster->raster_height);
                }
            }
        }
    }

    epr_free_record(record_beg);
    epr_free_record(record_end);
    epr_free_record(record);
    free(line_beg_buffer);
    free(line_end_buffer);
    return 0;
}


uint epr_get_raster_elem_size(const EPR_SRaster* raster) {
    if (raster == NULL) {
        epr_set_err(e_err_invalid_raster, "epr_get_raster_elem_size: raster must not be NULL");
        return 0;
    }
    return raster->elem_size;
}

void* epr_get_raster_elem_addr(const EPR_SRaster* raster, uint offset) {
    if (raster == NULL) {
        epr_set_err(e_err_invalid_raster, "epr_get_raster_elem_addr: raster must not be NULL");
        return 0;
    }
    return ((uchar*) raster->buffer) + epr_get_raster_elem_size(raster) * offset;
}

void* epr_get_raster_pixel_addr(const EPR_SRaster* raster, uint x, uint y) {
    if (raster == NULL) {
        epr_set_err(e_err_invalid_raster, "epr_get_raster_pixel_addr: raster must not be NULL");
        return 0;
    }
    return epr_get_raster_elem_addr(raster, y * raster->raster_width + x);
}

void* epr_get_raster_line_addr(const EPR_SRaster* raster, uint y) {
    if (raster == NULL) {
        epr_set_err(e_err_invalid_raster, "epr_get_raster_line_addr: raster must not be NULL");
        return 0;
    }
    return epr_get_raster_elem_addr(raster, y * raster->raster_width);
}

uint epr_get_raster_width(EPR_SRaster* raster) {
    if (raster == NULL) {
        epr_set_err(e_err_invalid_raster, "epr_get_raster_width: raster must not be NULL");
        return 0;
    }
    return raster->raster_width;
}

uint epr_get_raster_height(EPR_SRaster* raster) {
    if (raster == NULL) {
        epr_set_err(e_err_invalid_raster, "epr_get_raster_height: raster must not be NULL");
        return 0;
    }
    return raster->raster_height;
}

/******************************************************************/
EPR_FLineDecoder select_line_decode_function(EPR_EDataTypeId band_tid,
        EPR_ESampleModel band_smod,
        EPR_EDataTypeId raw_tid) {
    EPR_FLineDecoder decode_func;
    if ((band_tid == e_tid_char || band_tid == e_tid_uchar)
            && band_smod == e_smod_1OF1
            && (raw_tid == e_tid_char || raw_tid == e_tid_uchar))
        decode_func = decode_line_uchar_1_of_1_to_uchar;
    else if ((band_tid == e_tid_char || band_tid == e_tid_uchar)
             && band_smod == e_smod_1OF2
             && (raw_tid == e_tid_char || raw_tid == e_tid_uchar))
        decode_func = decode_line_uchar_1_of_2_to_uchar;
    else if ((band_tid == e_tid_char || band_tid == e_tid_uchar)
             && band_smod == e_smod_2OF2
             && (raw_tid == e_tid_char || raw_tid == e_tid_uchar))
        decode_func = decode_line_uchar_2_of_2_to_uchar;
    else if ((band_tid == e_tid_short || band_tid == e_tid_ushort)
             && band_smod == e_smod_1OF1
             && (raw_tid == e_tid_short || raw_tid == e_tid_ushort))
        decode_func = decode_line_ushort_1_of_1_to_ushort;
    else if (band_tid == e_tid_float
             && band_smod == e_smod_1OF1
             && raw_tid == e_tid_uchar)
        decode_func = decode_line_uchar_1_of_1_to_float;
    else if (band_tid == e_tid_float
             && band_smod == e_smod_1OF1
             && raw_tid == e_tid_char)
        decode_func = decode_line_char_1_of_1_to_float;
    else if (band_tid == e_tid_float
             && band_smod == e_smod_1OF1
             && raw_tid == e_tid_ushort)
        decode_func = decode_line_ushort_1_of_1_to_float;
    else if (band_tid == e_tid_float
             && band_smod == e_smod_1OF1
             && raw_tid == e_tid_short)
        decode_func = decode_line_short_1_of_1_to_float;
    else if (band_tid == e_tid_float
             && band_smod == e_smod_1OF2
             && raw_tid == e_tid_short)
        decode_func = decode_line_short_1_of_2_to_float;
    else if (band_tid == e_tid_float
             && band_smod == e_smod_2OF2
             && raw_tid == e_tid_short)
        decode_func = decode_line_short_2_of_2_to_float;
    else if (band_tid == e_tid_float
             && band_smod == e_smod_1OF2
             && raw_tid == e_tid_uchar)
        decode_func = decode_line_uchar_1_of_2_to_float;
    else if (band_tid == e_tid_float
             && band_smod == e_smod_2OF2
             && raw_tid == e_tid_uchar)
        decode_func = decode_line_uchar_2_of_2_to_float;
    else if (band_tid == e_tid_float
             && band_smod == e_smod_2TOF
             && raw_tid == e_tid_uchar)
        decode_func = decode_line_uchar_2_to_f_to_float;
    else if (band_tid == e_tid_uint
             && band_smod == e_smod_3TOI
             && raw_tid == e_tid_uchar)
        decode_func = decode_line_uchar_3_to_i_to_uint;
    else {
        return NULL;
    }
    return decode_func;
}


EPR_FArrayTransformer select_transform_array_function(EPR_EDataTypeId band_tid,
        EPR_EDataTypeId raw_tid) {
    EPR_FArrayTransformer transform_array_func;
    if (band_tid == e_tid_float && raw_tid == e_tid_short)
        transform_array_func = transform_array_short_to_float;
    else if (band_tid == e_tid_float && raw_tid == e_tid_ushort)
        transform_array_func = transform_array_ushort_to_float;
    else if (band_tid == e_tid_float && raw_tid == e_tid_int)
        transform_array_func = transform_array_int_to_float;
    else if (band_tid == e_tid_float && raw_tid == e_tid_uint)
        transform_array_func = transform_array_uint_to_float;
    else {
        return NULL;
    }
    return transform_array_func;
}


void decode_line_uchar_1_of_1_to_float(void* source_array,
                                       EPR_SBandId* band_id,
                                       int offset_x,
                                       int raster_width,
                                       int step_x,
                                       void* raster_buffer,
                                       int raster_pos) {
    int x, x1, x2;
    uchar* sa = (uchar*) source_array;
    float* buf = (float*) raster_buffer;

    x1 = offset_x;
    x2 = x1 + raster_width - 1;

    if (band_id->scaling_method == e_smid_log) {
        for (x = x1; x <= x2; x += step_x)  {
            buf[raster_pos++] = (float)pow(10, band_id->scaling_offset + band_id->scaling_factor * sa[x]);
        }
    } else if (band_id->scaling_method == e_smid_lin) {
        for (x = x1; x <= x2; x += step_x)  {
            buf[raster_pos++] = band_id->scaling_offset + band_id->scaling_factor * sa[x];
        }
    } else {
        for (x = x1; x <= x2; x += step_x)  {
            buf[raster_pos++] = sa[x];
        }
    }
}


void decode_line_char_1_of_1_to_float(void* source_array,
                                      EPR_SBandId* band_id,
                                      int offset_x,
                                      int raster_width,
                                      int step_x,
                                      void* raster_buffer,
                                      int raster_pos) {
    int x, x1, x2;
    char* sa = (char*) source_array;
    float* buf = (float*) raster_buffer;

    x1 = offset_x;
    x2 = x1 + raster_width - 1;

    if (band_id->scaling_method == e_smid_log) {
        for (x = x1; x <= x2; x += step_x)  {
            buf[raster_pos++] = (float)pow(10, band_id->scaling_offset + band_id->scaling_factor * sa[x]);
        }
    } else if (band_id->scaling_method == e_smid_lin) {
        for (x = x1; x <= x2; x += step_x)  {
            buf[raster_pos++] = band_id->scaling_offset + band_id->scaling_factor * sa[x];
        }
    } else {
        for (x = x1; x <= x2; x += step_x) {
            buf[raster_pos++] = sa[x];
        }
    }
}


void decode_line_ushort_1_of_1_to_float(void* source_array,
                                        EPR_SBandId* band_id,
                                        int offset_x,
                                        int raster_width,
                                        int step_x,
                                        void* raster_buffer,
                                        int raster_pos) {
    int x, x1, x2;
    ushort* sa = (ushort*) source_array;
    float* buf = (float*) raster_buffer;

    x1 = offset_x;
    x2 = x1 + raster_width - 1;

    if (band_id->scaling_method == e_smid_log) {
        for (x = x1; x <= x2; x += step_x) {
            buf[raster_pos++] = (float)pow(10, band_id->scaling_offset + band_id->scaling_factor * sa[x]);
        }
    } else if (band_id->scaling_method == e_smid_lin) {
        for (x = x1; x <= x2; x += step_x) {
            buf[raster_pos++] = band_id->scaling_offset + band_id->scaling_factor * sa[x];
        }
    } else {
        for (x = x1; x <= x2; x += step_x) {
            buf[raster_pos++] = sa[x];
        }
    }
}


void decode_line_short_1_of_1_to_float(void* source_array,
                                       EPR_SBandId* band_id,
                                       int offset_x,
                                       int raster_width,
                                       int step_x,
                                       void* raster_buffer,
                                       int raster_pos) {
    int x, x1, x2;
    short* sa = (short*) source_array;
    float* buf = (float*) raster_buffer;

    x1 = offset_x;
    x2 = x1 + raster_width - 1;

    if (band_id->scaling_method == e_smid_log) {
        for (x = x1; x <= x2; x += step_x)  {
            buf[raster_pos++] = (float)pow(10, band_id->scaling_offset + band_id->scaling_factor * sa[x]);
        }
    } else if (band_id->scaling_method == e_smid_lin) {
        for (x = x1; x <= x2; x += step_x)  {
            buf[raster_pos++] = band_id->scaling_offset + band_id->scaling_factor * sa[x];
        }
    } else {
        for (x = x1; x <= x2; x += step_x)  {
            buf[raster_pos++] = sa[x];
        }
    }
}

void decode_line_short_1_of_2_to_float(void* source_array,
                                       EPR_SBandId* band_id,
                                       int offset_x,
                                       int raster_width,
                                       int step_x,
                                       void* raster_buffer,
                                       int raster_pos) {
    int x, x1, x2;
    short* sa = (short*) source_array;
    float* buf = (float*) raster_buffer;

    x1 = offset_x;
    x2 = x1 + raster_width - 1;

    if (band_id->scaling_method == e_smid_log) {
        for (x = x1; x <= x2; x += step_x)  {
            buf[raster_pos++] = (float)pow(10, band_id->scaling_offset + band_id->scaling_factor * sa[2 * x]);
        }
    } else if (band_id->scaling_method == e_smid_lin) {
        for (x = x1; x <= x2; x += step_x)  {
            buf[raster_pos++] = band_id->scaling_offset + band_id->scaling_factor * sa[2 * x];
        }
    } else {
        for (x = x1; x <= x2; x += step_x)  {
            buf[raster_pos++] = (float)(sa[2 * x]);
        }
    }
}


void decode_line_short_2_of_2_to_float(void* source_array,
                                       EPR_SBandId* band_id,
                                       int offset_x,
                                       int raster_width,
                                       int step_x,
                                       void* raster_buffer,
                                       int raster_pos) {
    int x, x1, x2;
    short* sa = (short*) source_array;
    float* buf = (float*) raster_buffer;

    x1 = offset_x;
    x2 = x1 + raster_width - 1;

    if (band_id->scaling_method == e_smid_log) {
        for (x = x1; x <= x2; x += step_x) {
            buf[raster_pos++] = (float)pow(10, band_id->scaling_offset + band_id->scaling_factor * sa[2 * x + 1]);
        }
    } else if (band_id->scaling_method == e_smid_lin) {
        for (x = x1; x <= x2; x += step_x)  {
            buf[raster_pos++] = band_id->scaling_offset + band_id->scaling_factor * sa[2 * x + 1];
        }
    } else {
        for (x = x1; x <= x2; x += step_x)  {
            buf[raster_pos++] = (float)(sa[2 * x + 1]);
        }
    }
}


void decode_line_uchar_1_of_1_to_uchar(void* source_array,
                                       EPR_SBandId* band_id,
                                       int offset_x,
                                       int raster_width,
                                       int step_x,
                                       void* raster_buffer,
                                       int raster_pos) {
    int x, x1, x2;
    uchar* sa = (uchar*) source_array;
    uchar* buf = (uchar*) raster_buffer;

    x1 = offset_x;
    x2 = x1 + raster_width - 1;

    for (x = x1; x <= x2; x += step_x)  {
        buf[raster_pos++] = sa[x];
    }
}

void decode_line_uchar_1_of_2_to_uchar(void* source_array,
                                       EPR_SBandId* band_id,
                                       int offset_x,
                                       int raster_width,
                                       int step_x,
                                       void* raster_buffer,
                                       int raster_pos) {
    int x, x1, x2;
    uchar* sa = (uchar*) source_array;
    uchar* buf = (uchar*) raster_buffer;

    x1 = offset_x;
    x2 = x1 + raster_width - 1;

    for (x = x1; x <= x2; x += step_x)  {
        buf[raster_pos++] = sa[2 * x];
    }
}


void decode_line_uchar_2_of_2_to_uchar(void* source_array,
                                       EPR_SBandId* band_id,
                                       int offset_x,
                                       int raster_width,
                                       int step_x,
                                       void* raster_buffer,
                                       int raster_pos) {
    int x, x1, x2;
    uchar* sa = (uchar*) source_array;
    uchar* buf = (uchar*) raster_buffer;

    x1 = offset_x;
    x2 = x1 + raster_width - 1;

    for (x = x1; x <= x2; x += step_x)  {
        buf[raster_pos++] = sa[2 * x + 1];
    }
}

void decode_line_ushort_1_of_1_to_ushort(void* source_array,
        EPR_SBandId* band_id,
        int offset_x,
        int raster_width,
        int step_x,
        void* raster_buffer,
        int raster_pos) {
    int x, x1, x2;
    ushort* sa = (ushort*) source_array;
    ushort* buf = (ushort*) raster_buffer;

    x1 = offset_x;
    x2 = x1 + raster_width - 1;

    for (x = x1; x <= x2; x += step_x)  {
        buf[raster_pos++] = sa[x];
    }
}

void decode_line_uchar_2_to_f_to_float(void* source_array,
                                       EPR_SBandId* band_id,
                                       int offset_x,
                                       int raster_width,
                                       int step_x,
                                       void* raster_buffer,
                                       int raster_pos) {
    int x, x1, x2, shi;
    uchar* sa = (uchar*) source_array;
    float* buf = (float*) raster_buffer;

    x1 = offset_x;
    x2 = x1 + raster_width - 1;

    if (band_id->scaling_method == e_smid_log) {
        for (x = x1; x <= x2; x += step_x)  {
            shi = (((sa[2 * x] & 0xff)) | ((sa[2 * x + 1] & 0xff) << 8)) & 0xffff;
            buf[raster_pos++] = (float)pow(10, band_id->scaling_offset + band_id->scaling_factor * shi);
        }
    } else if (band_id->scaling_method == e_smid_lin) {
        for (x = x1; x <= x2; x += step_x)  {
            shi = (((sa[2 * x] & 0xff)) | ((sa[2 * x + 1] & 0xff) << 8)) & 0xffff;
            buf[raster_pos++] = band_id->scaling_offset + band_id->scaling_factor * shi;
        }
    } else {
        for (x = x1; x <= x2; x += step_x)  {
            shi = (((sa[2 * x] & 0xff)) | ((sa[2 * x + 1] & 0xff) << 8)) & 0xffff;
            buf[raster_pos++] = (float)shi;
        }
    }
}

void decode_line_uchar_1_of_2_to_float(void* source_array,
                                       EPR_SBandId* band_id,
                                       int offset_x,
                                       int raster_width,
                                       int step_x,
                                       void* raster_buffer,
                                       int raster_pos) {
    int x, x1, x2;
    uchar* sa = (uchar*) source_array;
    float* buf = (float*) raster_buffer;

    x1 = offset_x;
    x2 = x1 + raster_width - 1;

    if (band_id->scaling_method == e_smid_log) {
        for (x = x1; x <= x2; x += step_x)  {
            buf[raster_pos++] = (float)pow(10, band_id->scaling_offset + band_id->scaling_factor * (sa[2 * x] & 0xff));
        }
    } else if (band_id->scaling_method == e_smid_lin) {
        for (x = x1; x <= x2; x += step_x)  {
            buf[raster_pos++] = band_id->scaling_offset + band_id->scaling_factor * (sa[2 * x] & 0xff);
        }
    } else {
        for (x = x1; x <= x2; x += step_x)  {
            buf[raster_pos++] = (float)(sa[2 * x] & 0xff);
        }
    }
}



void decode_line_uchar_2_of_2_to_float(void* source_array,
                                       EPR_SBandId* band_id,
                                       int offset_x,
                                       int raster_width,
                                       int step_x,
                                       void* raster_buffer,
                                       int raster_pos) {
    int x, x1, x2;
    uchar* sa = (uchar*) source_array;
    float* buf = (float*) raster_buffer;

    x1 = offset_x;
    x2 = x1 + raster_width - 1;

    if (band_id->scaling_method == e_smid_log) {
        for (x = x1; x <= x2; x += step_x) {
            buf[raster_pos++] = (float)pow(10, band_id->scaling_offset + band_id->scaling_factor * (sa[2 * x + 1] & 0xff));
        }
    } else if (band_id->scaling_method == e_smid_lin) {
        for (x = x1; x <= x2; x += step_x)  {
            buf[raster_pos++] = band_id->scaling_offset + band_id->scaling_factor * (sa[2 * x + 1] & 0xff);
        }
    } else {
        for (x = x1; x <= x2; x += step_x)  {
            buf[raster_pos++] = (float)(sa[2 * x + 1] & 0xff);
        }
    }
}

void decode_line_uchar_3_to_i_to_uint(void* source_array,
                                       EPR_SBandId* band_id,
                                       int offset_x,
                                       int raster_width,
                                       int step_x,
                                       void* raster_buffer,
                                       int raster_pos) {
    int x, x1, x2, n;
    uchar* sa = (uchar*) source_array;
    uint* buf = (uint*) raster_buffer;

    x1 = offset_x;
    x2 = x1 + raster_width - 1;

    for (x = x1; x <= x2; x += step_x)  {
        n = 3 * x;
        buf[raster_pos++] = (sa[n] & 0xff) << 16 | ((sa[n+1] & 0xff) << 8) | ((sa[n+2] & 0xff) );
    }
}

/**
 * Computes the physical values for the annotation data.
 *
 * @param sa_beg the float array of tie points "before" Y-coordinate of the point to search
 * @param sa_end the float array of tie points "after" Y-coordinate of the point to search
 * @param samples_per_tie_pt the "distance" between two neighbour tie points (in scan-line direction)
 * @param num_elems number of elements in one tie point scan-line
 * @param band_id the information about properties and quantities of ENVISAT data.
 * @param offset_x [PIXEL] X-coordinate (0-based) of the upper right raster corner to search
 * @param y_mod [PIXEL] relative location of the point is being searched for (in flight direction)
 * @param raster_width [PIXEL] the width of the raster is being searched for
 * @param s_x [PIXEL] X-step to get the next point (in source coordinates) to search
 * @param raster_buffer the float user array to be filled with physical values
 * @param raster_pos the actual "filled" position in raster_buffer-array
 *
 */
void decode_tiepoint_band(float* sa_beg,
                          float* sa_end,
                          uint samples_per_tie_pt,
                          uint num_elems,
                          EPR_SBandId* band_id,
                          int offset_x,
                          float scan_offset_x,
                          float y_mod,
                          int raster_width,
                          int s_x,
                          float* raster_buffer,
                          int raster_pos) {

    int ix;
    float x_mod;
    uint x_knot;
    int inti_flag = 0;
    int intersection = 0;
    float circle;
    float half_circle;
    float null_point;

    circle = EPR_LONGI_ABS_MAX - EPR_LONGI_ABS_MIN;
    half_circle = 0.5F * circle;
    null_point = 0.5F * (EPR_LONGI_ABS_MAX + EPR_LONGI_ABS_MIN);

    if (strncmp(band_id->band_name, EPR_LONGI_BAND_NAME, strlen(EPR_LONGI_BAND_NAME)) == 0) {
        inti_flag = 1;
    } else {
        inti_flag = 0;
    }

    for (ix = offset_x; ix < offset_x + raster_width; ix += s_x) {
        x_mod = (ix - scan_offset_x) / samples_per_tie_pt;
        if (x_mod >= 0.0F) {
            x_knot = (uint)floor(x_mod);
            if (x_knot >= num_elems - 1) {
                x_knot = num_elems - 2;
            }
        } else {
            x_knot = (uint)0;
        }
        x_mod -= x_knot;

        if (inti_flag == 1) {
            if (fabs((float)(sa_beg[x_knot + 1] - sa_beg[x_knot])) > half_circle ||
                    fabs((float)(sa_beg[x_knot] - sa_end[x_knot])) > half_circle ||
                    fabs((float)(sa_end[x_knot] - sa_end[x_knot + 1])) > half_circle ||
                    fabs((float)(sa_end[x_knot + 1] - sa_beg[x_knot + 1])) > half_circle) {
                intersection = 1;
                if (sa_beg[x_knot] < (float)null_point) {
                    sa_beg[x_knot] += circle;
                }
                if (sa_beg[x_knot + 1] < (float)null_point) {
                    sa_beg[x_knot + 1] += circle;
                }
                if (sa_end[x_knot] < (float)null_point) {
                    sa_end[x_knot] += circle;
                }
                if (sa_end[x_knot + 1] < (float)null_point) {
                    sa_end[x_knot + 1] += circle;
                }
            }
        } else {
            intersection = 0;
        }

        raster_buffer[raster_pos] = epr_interpolate2D(x_mod, y_mod,
                                    sa_beg[x_knot], sa_beg[x_knot + 1],
                                    sa_end[x_knot], sa_end[x_knot + 1]);

        if (inti_flag == 1 &&
                intersection == 1 &&
                raster_buffer[raster_pos] > EPR_LONGI_ABS_MAX) {
            raster_buffer[raster_pos] -= circle;
        }

        raster_pos++;
    }
}


float epr_interpolate2D(float wi, float wj, float x00, float x10, float x01, float x11) {
    return x00 + wi * (x10 - x00) + wj * (x01 - x00) + wi * wj * (x11 + x00 - x01 - x10);
}


void transform_array_short_to_float (void* sourceArray,
                                     EPR_SBandId* band_id,
                                     float* raster_buffer,
                                     uint nel) {
    uint ix;
    short* sa = (short*) sourceArray;

    for (ix = 0; ix < nel; ix ++) {
        raster_buffer[ix] = band_id->scaling_offset + band_id->scaling_factor * sa[ix];
    }
}

void transform_array_ushort_to_float (void* sourceArray,
                                      EPR_SBandId* band_id,
                                      float* raster_buffer,
                                      uint nel) {
    uint ix;
    ushort* sa = (ushort*) sourceArray;

    for (ix = 0; ix < nel; ix ++) {
        raster_buffer[ix] = band_id->scaling_offset + band_id->scaling_factor * sa[ix];
    }
}

void transform_array_int_to_float (void* sourceArray,
                                    EPR_SBandId* band_id,
                                    float* raster_buffer,
                                    uint nel) {
    uint ix;
    int* sa = (int*) sourceArray;

    for (ix = 0; ix < nel; ix ++) {
        raster_buffer[ix] = band_id->scaling_offset + band_id->scaling_factor * sa[ix];
    }
}

void transform_array_uint_to_float (void* sourceArray,
                                     EPR_SBandId* band_id,
                                     float* raster_buffer,
                                     uint nel) {
    uint ix;
    uint* sa = (uint*) sourceArray;

    for (ix = 0; ix < nel; ix ++) {
        raster_buffer[ix] = band_id->scaling_offset + band_id->scaling_factor * sa[ix];
    }
}

void mirror_float_array(float* raster_buffer, uint raster_width, uint raster_height) {
    uint w, h, pol_w, offset;
    float tmp;
    pol_w = raster_width / 2;

    for (h = 0; h < raster_height; h ++) {
        for (w = 0; w < pol_w; w ++) {
            offset = h * raster_width;
            tmp = raster_buffer[w + offset];
            raster_buffer[w + offset] = raster_buffer[raster_width - 1 - w + offset];
            raster_buffer[raster_width - 1 - w + offset] = tmp;
        }
    }
}

void mirror_uchar_array(uchar* raster_buffer, uint raster_width, uint raster_height) {
    uint w, h, pol_w, offset;
    uchar tmp;
    pol_w = raster_width / 2;

    for (h = 0; h < raster_height; h ++) {
        for (w = 0; w < pol_w; w ++) {
            offset = h * raster_width;
            tmp = raster_buffer[w + offset];
            raster_buffer[w + offset] = raster_buffer[raster_width - 1 - w + offset];
            raster_buffer[raster_width - 1 - w + offset] = tmp;
        }
    }
}

void mirror_ushort_array(ushort* raster_buffer, uint raster_width, uint raster_height) {
    uint w, h, pol_w, offset;
    ushort tmp;
    pol_w = raster_width / 2;

    for (h = 0; h < raster_height; h ++) {
        for (w = 0; w < pol_w; w ++) {
            offset = h * raster_width;
            tmp = raster_buffer[w + offset];
            raster_buffer[w + offset] = raster_buffer[raster_width - 1 - w + offset];
            raster_buffer[raster_width - 1 - w + offset] = tmp;
        }
    }
}

void mirror_uint_array(uint* raster_buffer, uint raster_width, uint raster_height) {
    uint w, h, pol_w, offset;
    uint tmp;
    pol_w = raster_width / 2;

    for (h = 0; h < raster_height; h ++) {
        for (w = 0; w < pol_w; w ++) {
            offset = h * raster_width;
            tmp = raster_buffer[w + offset];
            raster_buffer[w + offset] = raster_buffer[raster_width - 1 - w + offset];
            raster_buffer[raster_width - 1 - w + offset] = tmp;
        }
    }
}

void epr_zero_invalid_pixels(EPR_SRaster* raster, EPR_SRaster* bm_raster) {

    uchar* bm_pixels;
    uint bm_pos = 0;
    uint bm_len = raster->raster_width * raster->raster_height;

    assert(bm_raster->data_type == e_tid_char
           || bm_raster->data_type == e_tid_uchar);


    bm_pixels = (uchar*) bm_raster->buffer;
    switch (raster->data_type) {
        case e_tid_char:
        case e_tid_uchar: {
            char* pixels = (char*) raster->buffer;
            for (bm_pos = 0; bm_pos < bm_len; bm_pos++) {
                if (bm_pixels[bm_pos] == 0) {
                    pixels[bm_pos] = 0;
                }
            }
        }
        break;
        case e_tid_short:
        case e_tid_ushort: {
            short* pixels = (short*) raster->buffer;
            for (bm_pos = 0; bm_pos < bm_len; bm_pos++) {
                if (bm_pixels[bm_pos] == 0) {
                    pixels[bm_pos] = 0;
                }
            }
        }
        break;
        case e_tid_int:
        case e_tid_uint: {
            int* pixels = (int*) raster->buffer;
            for (bm_pos = 0; bm_pos < bm_len; bm_pos++) {
                if (bm_pixels[bm_pos] == 0) {
                    pixels[bm_pos] = 0;
                }
            }
        }
        break;
        case e_tid_float: {
            float* pixels = (float*) raster->buffer;
            for (bm_pos = 0; bm_pos < bm_len; bm_pos++) {
                if (bm_pixels[bm_pos] == 0) {
                    pixels[bm_pos] = 0.0F;
                }
            }
        }
        break;
        case e_tid_double: {
            double* pixels = (double*) raster->buffer;
            for (bm_pos = 0; bm_pos < bm_len; bm_pos++) {
                if (bm_pixels[bm_pos] == 0) {
                    pixels[bm_pos] = 0.0;
                }
            }
        }
        break;
        default: {}
        break;
    }
}
