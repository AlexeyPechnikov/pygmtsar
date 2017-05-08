/*
 * $Id: epr_record.c,v 1.1.1.1 2004-10-28 19:22:23 norman Exp $
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

/*
 * ===================== Record Info Access ==============================
 */


/*
   Function:    epr_create_record_info
   Access:      public API
   Changelog:   2002/01/16  mp initial version
 */
/**
 * Creates a new record_info for the record by the given
 * dataset name.
 *
 * @param dataset_name the name of the dataset, to which the record
 *        belongs to, must not be <code>NULL</code>
 * @param field_infos the pointer at the strucrure with information
 *        of all fields wich belong to record,
 *          must not be <code>NULL</code>
 * @return the new record instance
 *         or <code>NULL</code> if an error occured.
 */
EPR_SRecordInfo* epr_create_record_info(const char* dataset_name, EPR_SPtrArray* field_infos)
{

    EPR_SRecordInfo* record_info = NULL;
    EPR_SFieldInfo* field_info = NULL;
    int field_infos_index;
    int field_infos_length;
    uint tot_record_size = 0;

    record_info = (EPR_SRecordInfo*) calloc(1, sizeof (EPR_SRecordInfo));
    if (record_info == NULL)
    {
        epr_set_err(e_err_out_of_memory,
                    "epr_create_record_info: out of memory");
        return NULL;
    }
    epr_assign_string(&record_info->dataset_name, dataset_name);
    if (record_info->dataset_name == NULL)
    {
        epr_set_err(e_err_out_of_memory,
                    "epr_create_record_info: out of memory");
        return NULL;
    }

    record_info->field_infos = field_infos;

    field_infos_length = (int)epr_get_ptr_array_length(field_infos);
    for (field_infos_index = 0; field_infos_index < field_infos_length; field_infos_index++)
    {
        field_info = (EPR_SFieldInfo*)epr_get_ptr_array_elem_at(field_infos, field_infos_index);
        tot_record_size += field_info->tot_size;
    }

    record_info->tot_size = tot_record_size;

    return record_info;
}


/**
 * Frees the memory allocated by the given record_info.
 *
 * <p> After calling this function the give record_info pointer gets
 * invalid and should not be used anymore.
 *
 * @param record_info the record to be released, if <code>NULL</code>
 *        the function immediately returns
 */
void epr_free_record_info(EPR_SRecordInfo* record_info)
{
    EPR_SFieldInfo* field_info = NULL;
    uint field_info_index = 0;

    if (record_info == NULL)
        return;

    if (record_info->field_infos != NULL)
    {
        for (field_info_index = 0; field_info_index < record_info->field_infos->length; field_info_index++)
        {
            field_info = (EPR_SFieldInfo*)epr_get_ptr_array_elem_at(record_info->field_infos, field_info_index);
            epr_free_field_info(field_info);
        }
        epr_free_ptr_array(record_info->field_infos);
        record_info->field_infos = NULL;
    }

    epr_free_string(record_info->dataset_name);
    record_info->dataset_name = NULL;

    record_info->field_infos = NULL;
    record_info->tot_size = 0;

    free(record_info);
}


/*
   Function:    epr_get_record_info
   Access:      private API implementation helper
   Changelog:   2002/01/05  mp initial version
 */
/**
 * Returns information about the structure of the records contained in
 * a dataset specified by the given <code>dataset_id</code>.
 *
 * @param dataset_id the the dataset identifier
 *
 * @return the the pointer for the record structure information.
 */
EPR_SRecordInfo* epr_get_record_info(EPR_SDatasetId* dataset_id)
{
    EPR_SProductId* product_id;
    EPR_SRecordInfo* record_info;
    uint num_record_infos;
    uint record_index = 0;

    assert(dataset_id != NULL);

    product_id = dataset_id->product_id;
    assert(product_id != NULL);
    assert(product_id->record_info_cache != NULL);

    /*
     *  checks, if 'record_info_cache' is not empty, then returns this cache,
     *  not making moreover
     */
    num_record_infos = product_id->record_info_cache->length;
    for (record_index = 0; record_index < num_record_infos; record_index++)
    {
        record_info = (EPR_SRecordInfo*) product_id->record_info_cache->elems[record_index];
        if (epr_equal_names(record_info->dataset_name, dataset_id->dataset_name))
            return record_info;
    }

	record_info = epr_read_record_info(product_id, dataset_id);
    if (record_info == NULL) {
        epr_set_err(e_err_file_access_denied,
                   "epr_get_record_info: invalid DDDB resource path: missing any 'ASAR' files");
        return NULL;
    }

    epr_add_ptr_array_elem(product_id->record_info_cache, record_info);

    return record_info;
}



EPR_SRecordInfo* epr_read_record_info(EPR_SProductId* product_id, EPR_SDatasetId* dataset_id)
{
    EPR_SRecordInfo* record_info = NULL;
    EPR_SFieldInfo* field_info = NULL;
    EPR_SPtrArray* field_infos = NULL;
    EPR_EDataTypeId data_type_id = e_tid_unknown;
    uint num_elems = 0;
    uint num_bytes = 0;
    uint more_count = 0;
    char* field_name = NULL;
    char* data_type = NULL;
    char* unit = NULL;
    char* description = NULL;
	int i;
    int rt_index;
    const struct RecordDescriptorTable* r_tables;
    int num_descr;
	int num_r_tables;

    if (product_id == NULL) {
        epr_set_err(e_err_null_pointer,
                    "epr_read_record_info: product_id must not be NULL");
        return NULL;
    }

	/* @DDDB */

    if (strncmp(product_id->id_string, "MER", 3) == 0) {
		r_tables = dddb_meris_rec_tables;
		num_r_tables = EPR_NUM_MERIS_REC_TABLES;
	} else if (strncmp(product_id->id_string, "ATS", 3) == 0) {
		r_tables = dddb_aatsr_rec_tables;
		num_r_tables = EPR_NUM_AATSR_REC_TABLES;
	} else if (strncmp(product_id->id_string, "ASA", 3) == 0) {
		r_tables = dddb_asar_rec_tables;
		num_r_tables = EPR_NUM_ASAR_REC_TABLES;
	} else if (strncmp(product_id->id_string, "SAR", 3) == 0) {
		r_tables = dddb_asar_rec_tables;
		num_r_tables = EPR_NUM_ASAR_REC_TABLES;
	} else {
        epr_set_err(e_err_invalid_product_id,
                    "epr_read_record_info: invalid product identifier");
        return NULL;
	}

    rt_index = -1;
    for (i = 0; i < num_r_tables; i++) {
        if (dataset_id->record_descriptor == r_tables[i].descriptors) {
            rt_index = i;
            break;
        }
    }
    if (rt_index == -1) {
        epr_set_err(e_err_invalid_product_id,
                    "epr_read_record_info: unknown record");
        return NULL;
    }

    field_infos = epr_create_ptr_array(16);
    num_descr = r_tables[rt_index].num_descriptors;
    for (i = 0; i < num_descr; i++) {
        /* 1: field_name */
		field_name = epr_clone_string(r_tables[rt_index].descriptors[i].id);
        /* 2: data_type_id */
		data_type_id = r_tables[rt_index].descriptors[i].type;
        /* 3: unit */
		unit = epr_clone_string(r_tables[rt_index].descriptors[i].unit);
        /* 4: num_elems */
		num_bytes = r_tables[rt_index].descriptors[i].elem_size;
        /* 5: num_elems and more_count*/
        /* @todo: check return value! and check epr_parse_value_count */
		num_elems = epr_parse_value_count(product_id, r_tables[rt_index].descriptors[i].num_elem);
		more_count = 1;
        /* 6: description*/
		description = epr_clone_string(r_tables[rt_index].descriptors[i].description);

		field_info = epr_create_field_info(data_type_id, description, field_name, num_elems, num_bytes, more_count, unit);
        epr_add_ptr_array_elem(field_infos, field_info);

        epr_free_string(field_name);
        epr_free_string(data_type);
        epr_free_string(unit);
        epr_free_string(description);
    }
    record_info = epr_create_record_info(dataset_id->dataset_name, field_infos);

    return record_info;
}


/*
   Function:    epr_read_record_info
   Access:      private API implementation helper
   Changelog:   2002/01/05  mp initial version
 */
/**
 * Reads the record information with the given file path and
 * returns the poiter at it.
 *
 * @param product_id the the product file identifier
 * @param dataset_name the name of the dataset
 * @param db_file_istream the DB-table file identifier
 *
 * @return the the pointer at the record information.
 */
 /*
EPR_SRecordInfo* epr_read_record_info(EPR_SProductId* product_id, EPR_SDatasetId* dataset_id, FILE* db_file_istream)

}
*/

/*
   Function:    epr_create_record
   Access:      public API
   Changelog:   2002/01/23  mp initial version
 */
/**
 * Creates a new record instance from the dataset specified by the given
 * dataset name.
 *
 * @param the pointer at the record information.
 *        must not be <code>NULL</code>
 * @return the new record instance
 *         or <code>NULL</code> if an error occured.
 */
EPR_SRecord* epr_create_record_from_info(EPR_SRecordInfo* record_info)
{
    EPR_SRecord* record = NULL;
    EPR_SFieldInfo* field_info = NULL;
    uint field_infos_index = 0;
    int field_infos_length;

    if (record_info == NULL)
    {
        epr_set_err(e_err_out_of_memory,
                    "epr_create_record_from_info: out of memory");
        return NULL;
    }

    record = (EPR_SRecord*) calloc(1, sizeof(EPR_SRecord));
    if (record == NULL)
    {
        epr_set_err(e_err_out_of_memory,
                    "epr_create_record: out of memory");
        return NULL;
    }

    record->magic = EPR_MAGIC_RECORD;
    record->info = record_info;
    record->num_fields = record_info->field_infos->length;

    record->fields = (EPR_SField**) calloc(record->num_fields, sizeof(EPR_SField*));
    if (record->fields == NULL)
    {
        epr_set_err(e_err_out_of_memory,
                    "epr_create_record: out of memory");
        return NULL;
    }

    field_infos_length = (int)epr_get_ptr_array_length(record_info->field_infos);

    for (field_infos_index = 0; field_infos_index < record_info->field_infos->length; field_infos_index++)
    {
        field_info = (EPR_SFieldInfo*)epr_get_ptr_array_elem_at(record_info->field_infos, field_infos_index);
        record->fields[field_infos_index] = epr_create_field(field_info);
    }
    return record;
}


const EPR_SField* epr_get_field_at(const EPR_SRecord* record, uint field_index){

    EPR_SField* field = NULL;

    epr_clear_err();

    if (record == NULL) {
        epr_set_err(e_err_null_pointer,
                    "epr_get_field_at: record-identifier must not be NULL");
        return NULL;
    }

    if ((field_index < 0) || (field_index >= record->num_fields))
    {
        epr_set_err(e_err_index_out_of_range,
                "epr_get_field_at: field_index out of range");
        return NULL;
    }
    field = record->fields[field_index];
    return field;
}


uint epr_get_num_fields(const EPR_SRecord* record)
{
    epr_clear_err();

    if (record == NULL) {
        epr_set_err(e_err_null_pointer,
                    "epr_get_field_at: record-identifier must not be NULL");
        return (uint)-1;
    }
    return record->num_fields;
}


/**
 * Frees the memory allocated by the given record.
 *
 * <p> After calling this function the give record pointer gets
 * invalid and should not be used anymore.
 *
 * @param record the record to be released, if <code>NULL</code>
 *        the function immediately returns
 */
void epr_free_record(EPR_SRecord* record)
{
    EPR_SField* field = NULL;
    uint fields_index = 0;

    epr_clear_err();

    if (record == NULL)
        return;

    if (record->fields != NULL)
    {
        for (fields_index = 0; fields_index < record->num_fields; fields_index++)
        {
            field = (EPR_SField*)record->fields[fields_index];
            epr_free_field(field);
        }
        free(record->fields);
        record->fields = NULL;
    }

    /* Do NOT free record->info since many records can
       share the same record->info! */
    record->info = NULL;

    record->num_fields = 0;

    free(record);
}
