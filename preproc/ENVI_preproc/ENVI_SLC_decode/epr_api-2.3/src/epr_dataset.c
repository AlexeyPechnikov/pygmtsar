/*
 * $Id: epr_dataset.c,v 1.1.1.1 2004-10-28 19:22:22 norman Exp $
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
#include "epr_dataset.h"
#include "epr_record.h"
#include "epr_param.h"
#include "epr_dsd.h"
#include "epr_msph.h"
#include "epr_band.h"
#include "epr_bitmask.h"

#include "epr_dddb.h"


EPR_SDatasetId* epr_create_dataset_id(EPR_SProductId* product_id,
                                      const EPR_SDSD* dsd,
                                      const char* dataset_name,
                                      const struct RecordDescriptor* record_descriptor,
                                      const char* dsd_name,
                                      const char* description)
{
    EPR_SDatasetId* dataset_id = (EPR_SDatasetId*) calloc(1, sizeof (EPR_SDatasetId));
    if (dataset_id == NULL) {
        epr_set_err(e_err_out_of_memory,
                    "epr_create_dataset_id: out of memory");
        return NULL;
    }
    dataset_id->magic = EPR_MAGIC_DATASET_ID;

    dataset_id->product_id = product_id;
    dataset_id->dsd = dsd;
    dataset_id->record_info = NULL;

    epr_assign_string(&dataset_id->dataset_name, dataset_name);
    dataset_id->record_descriptor = record_descriptor;
    epr_assign_string(&dataset_id->dsd_name, dsd_name);
    epr_assign_string(&dataset_id->description, description);

    return dataset_id;
}

/**
 * Release the memory allocated through a dataset ID.
 *
 * @param product_id the file identifier, if <code>NULL</code> the function
 *        immediately returns zero.
 * @return zero for success, an error code otherwise
 */
void epr_free_dataset_id(EPR_SDatasetId* dataset_id)
{
    if (dataset_id == NULL)
        return;

    /* Don't free the product ID, it is NOT owned by a dataset ID */
    dataset_id->product_id = NULL;
    /* Don't free the record information, it is NOT owned by a dataset ID */
    dataset_id->record_info = NULL;
    /* Don't free the DSD, it is NOT owned by a dataset ID */
    dataset_id->dsd = NULL;

    epr_free_and_null_string(&dataset_id->dataset_name);
    epr_free_and_null_string(&dataset_id->dsd_name);
    epr_free_and_null_string(&dataset_id->description);

    free(dataset_id);
}


/**
 * Creates an array of dataset_id for the given ENVISAT product
 *
 * @param product_id the the product file identifier
 * @return the instance of the array
 */
EPR_SPtrArray* epr_create_dataset_ids(EPR_SProductId* product_id)
{
    EPR_SPtrArray* dataset_ids = NULL;
    EPR_SDatasetId* dataset_id = NULL;
    const EPR_SDSD* dsd = NULL;
    uint dsd_index;
    int i;
    const struct DatasetDescriptorTable* p_tables;
    int pt_index;
    int num_descr;

    if (product_id == NULL) {
        epr_set_err(e_err_null_pointer,
                    "create_dataset_ids: product_id must not be NULL");
        return NULL;
    }

    /* @DDDB */

    p_tables = dddb_product_tables;
    pt_index = -1;
    for (i = 0; i < EPR_NUM_PRODUCT_TABLES; i++) {
        const char* id = p_tables[i].name;
        if (strncmp(product_id->id_string, id, 10) == 0) {
            if (product_id->meris_iodd_version == 5) {
                if (strcmp(id, "MER_RR__1P_IODD5") == 0 ||
                    strcmp(id, "MER_FR__1P_IODD5") == 0) {
                    pt_index = i;
                }
            } else if (product_id->meris_iodd_version == 6) {
                if (strcmp(id, "MER_RR__2P_IODD6") == 0 ||
                    strcmp(id, "MER_FR__2P_IODD6") == 0) {
                    pt_index = i;
                }
            } else {
                pt_index = i;
            }
        }
        if (pt_index != -1) {
            break;
        }
    }
    if (pt_index == -1) {
        epr_set_err(e_err_null_pointer,
                    "create_dataset_ids: unknown product type");
        return NULL;
    }

	dataset_ids = epr_create_ptr_array(16);
    num_descr = p_tables[pt_index].num_descriptors;
    for (i = 0; i < num_descr; i++) {
        for (dsd_index = 0; dsd_index < product_id->dsd_array->length; dsd_index++) {
            dsd = (EPR_SDSD*)epr_get_ptr_array_elem_at(product_id->dsd_array, dsd_index);
            if (strncmp(dsd->ds_name, p_tables[pt_index].descriptors[i].ds_name, strlen(epr_strip_string_r(dsd->ds_name))) == 0) {
                dataset_id = epr_create_dataset_id(product_id,
                                                    dsd,
                                                    p_tables[pt_index].descriptors[i].id,
                                                    p_tables[pt_index].descriptors[i].rec_descriptor,
                                                    p_tables[pt_index].descriptors[i].ds_name,
                                                    p_tables[pt_index].descriptors[i].description);
                epr_add_ptr_array_elem(dataset_ids, dataset_id);
                break;
            }
        }
	}

    return dataset_ids;
}


const char* epr_get_dataset_name(EPR_SDatasetId* dataset_id)
{
    epr_clear_err();

    if (dataset_id == NULL) {
        epr_set_err(e_err_null_pointer,
                    "epr_get_dataset_name: band_id must not be NULL");
        return NULL;
    }
    return epr_trim_string(dataset_id->dataset_name);
}

uint epr_get_num_records(const EPR_SDatasetId* dataset_id)
{
    epr_clear_err();

    if (dataset_id == NULL) {
        epr_set_err(e_err_invalid_dataset_name,
                    "epr_get_num_records: invalid dataset name");
        return (uint)-1;
    }
    return dataset_id->dsd->num_dsr;
}

const EPR_SDSD* epr_get_dsd(const EPR_SDatasetId* dataset_id)
{
    epr_clear_err();

    if (dataset_id == NULL) {
        epr_set_err(e_err_invalid_dataset_name,
                    "epr_get_dsd: invalid dataset name");
        return NULL;
    }
    return dataset_id->dsd;
}

const char* epr_get_dsd_name(const EPR_SDatasetId* dataset_id)
{
    epr_clear_err();

    if (dataset_id == NULL) {
        epr_set_err(e_err_invalid_dataset_name,
                    "epr_get_dsd_name: invalid dataset name");
        return NULL;
    }
    return epr_trim_string(dataset_id->dsd_name);
}

uint epr_get_dataset_offset(EPR_SDatasetId* dataset_id)
{
    if (dataset_id == NULL) {
        epr_set_err(e_err_null_pointer,
                    "epr_get_dataset_offset: dataset_id must not be NULL");
        return (uint)0;
    }
    return dataset_id->dsd->ds_offset;
}

/*********************************** RECORD ***********************************/

/*
   Function:    epr_create_record
   Access:      public API
   Changelog:   2002/01/23  mp initial version
 */
/**
 * Creates a new record for the dataset given by the specified
 * dataset identifier.
 */
EPR_SRecord* epr_create_record(EPR_SDatasetId* dataset_id)
{
    EPR_SRecord* record = NULL;

    epr_clear_err();

    if (dataset_id == NULL) {
        epr_set_err(e_err_illegal_arg,
                    "epr_create_record: dataset ID must not be NULL");
        return NULL;
    }
    if (dataset_id->record_info == NULL) {
        dataset_id->record_info = epr_get_record_info(dataset_id);
    }

    record = epr_create_record_from_info(dataset_id->record_info);
    if (record == NULL) {
        epr_set_err(e_err_invalid_record_name,
                    "epr_create_record: invalid record name");
        return NULL;
    }
    return record;
}


/**
 * Reads a full record from ENVISAT product file.
 */
EPR_SRecord* epr_read_record(EPR_SDatasetId* dataset_id,
                             uint record_index,
                             EPR_SRecord* record)
{
    uint field_index;
    uint dsd_offset;
    uint record_size;
    uint data_type_size;
    uint elements_to_read;
    uint elements_read;
    EPR_SField* field = NULL;

    epr_clear_err();

    if (dataset_id == NULL) {
        epr_set_err(e_err_invalid_dataset_name,
                    "epr_read_record: invalid dataset name");
        return NULL;

    }
    if ((record_index < 0) || (record_index >= dataset_id->dsd->num_dsr)) {
        epr_set_err(e_err_invalid_value,
                    "epr_read_record: invalid record_index parameter, must be >=0 and <num_dsr");
        return NULL;
    }

    if (record == NULL) {
        record = epr_create_record(dataset_id);
    } else if (record->info != dataset_id->record_info) {
        epr_set_err(e_err_invalid_record_name,
                    "epr_read_record: invalid record name");
        return NULL;
    }

    /*READ FILE with: dataset_id->product_id->istream after the setting it to begin*/
    rewind(dataset_id->product_id->istream);

    /*GET OFFSET*/
    dsd_offset = dataset_id->dsd->ds_offset;

    record_size = record->info->tot_size;
    if (record_size != dataset_id->dsd->dsr_size) {
        //epr_set_err(e_err_invalid_data_format,
        //    "epr_read_record: wrong record size");
        //return NULL;
    }

    /* Set file pointer to begin of demanded record */
    if (fseek(dataset_id->product_id->istream, dsd_offset + record_size * record_index, SEEK_SET) != 0) {
        epr_set_err(e_err_file_access_denied,
            "epr_read_record: file seek failed");
        return NULL;
    }

    for (field_index = 0; field_index < record->num_fields; field_index++) {
        field = record->fields[field_index];
        elements_to_read = field->info->num_elems ;
        data_type_size = epr_get_data_type_size(field->info->data_type_id);
        assert(data_type_size != 0);
        assert(field->elems != NULL);

        if (elements_to_read * data_type_size != field->info->tot_size) {
            /*epr_log(e_log_info, "Spare");*/
            data_type_size = field->info->tot_size / elements_to_read;
        }

        elements_read = fread(field->elems, data_type_size, elements_to_read, dataset_id->product_id->istream);
        if (elements_read != elements_to_read) {
            epr_set_err(e_err_file_read_error,
                "epr_read_record: file read failed");
            return NULL;
        }

        /*
         * SWAP bytes on little endian (LE) order architectures (I368, Pentium Processors).
         * ENVISAT products are stored in big endian (BE) order.
         */
        if (epr_api.little_endian_order) {
            epr_swap_endian_order(field);
        }
    }
    return record;
}
