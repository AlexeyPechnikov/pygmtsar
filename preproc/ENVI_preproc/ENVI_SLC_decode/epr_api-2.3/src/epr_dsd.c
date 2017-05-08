/*
* $Id: epr_dsd.c,v 1.1.1.1 2004-10-28 19:22:23 norman Exp $
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


/**
* Opens dsd for a dataset description,
* obtained from an ENVISAT product file.
*
* @param dsd_index the number of dsd (zero-based), emrty dsd inclusive
*
* @return the the pointer at the dsd information.
*/
EPR_SDSD* epr_create_dsd(int dsd_index)
{
    EPR_SDSD* dsd;
    dsd = (EPR_SDSD*) calloc(1, sizeof (EPR_SDSD));
    if (dsd == NULL) {
        epr_set_err(e_err_out_of_memory,
            "epr_create_dsd: out of memory");
        return NULL;
    }
    dsd->index = dsd_index;
    return dsd;
}



uint epr_get_num_datasets(EPR_SProductId* product_id)
{
    epr_clear_err();

    if (product_id == NULL) {
        epr_set_err(e_err_null_pointer,
            "epr_get_num_datasets: product_id must not be NULL");
        return (uint)-1;
    }
    return product_id->dataset_ids->length;
}

EPR_SDatasetId* epr_get_dataset_id_at(EPR_SProductId* product_id, uint index)
{
    EPR_SDatasetId* dataset_id = NULL;

    epr_clear_err();

    if (product_id == NULL) {
        epr_set_err(e_err_null_pointer,
            "epr_get_dataset_id_at: product_id must not be NULL");
        return NULL;
    }
    if (index >= product_id->dataset_ids->length){
        epr_set_err(e_err_index_out_of_range,
            "epr_get_dataset_id_at: dataset index out of range");
        return NULL;
    }

    dataset_id = (EPR_SDatasetId*)epr_get_ptr_array_elem_at(product_id->dataset_ids, index);
    return dataset_id;
}

EPR_SDatasetId* epr_get_dataset_id(EPR_SProductId* product_id, const char* dataset_name)
{
    EPR_SDatasetId* dataset_id = NULL;
    int datasets_num, i;

    epr_clear_err();

    if (product_id == NULL) {
        epr_set_err(e_err_null_pointer,
            "epr_get_dataset_id: product_id must not be NULL");
        return NULL;
    }
    if (dataset_name == NULL) {
        epr_set_err(e_err_null_pointer,
            "epr_get_dataset_id: dataset_name must not be NULL");
        return NULL;
    }

    datasets_num = epr_get_num_datasets(product_id);
    for (i = 0; i < datasets_num; i ++) {
        dataset_id = epr_get_dataset_id_at(product_id, i);
        if (epr_equal_names(dataset_name, epr_get_dataset_name(dataset_id)))
            return dataset_id;
    }
    epr_set_err(e_err_invalid_band_name,
        "epr_get_dataset_id: dataset_id not found");
    return NULL;
}

int epr_detect_meris_iodd_version(EPR_SProductId* product_id)
{
    EPR_SDSD** elems;
    int size = 0;
    int rec_size = 0;
    int i, ioddFormat = 7;
    char* name;

    /* reflect L1b product format change from IODD5 to IODD6 */
    if (strncmp("MER_RR__1P", product_id->id_string, 10) == 0
        || strncmp("MER_FR__1P", product_id->id_string, 10) == 0) {

        elems = (EPR_SDSD**)product_id->dsd_array->elems;
        size = product_id->dsd_array->length;
        for (i=0; i< size;i++){
            name = elems[i]->ds_name;
            if (strcmp("Flags MDS(16)", name) == 0) {
                rec_size = elems[i]->dsr_size;
                if (rec_size == 2255 || rec_size == 4495) {
                    ioddFormat =  5;
                }
                break;
            }
        }
    }
    /* reflect L2 product format change from IODD6 to IODD7 */
    else if (strncmp("MER_RR__2P", product_id->id_string, 10) == 0
        || strncmp("MER_FR__2P", product_id->id_string, 10) == 0) {

        elems = (EPR_SDSD**)product_id->dsd_array->elems;
        size = product_id->dsd_array->length;
        for (i=0; i<size; i++){
            name = elems[i]->ds_name;
            if (strcmp("Epsilon, OPT   - MDS(19)", name) == 0) {
                ioddFormat = 6;
                break;
            }
        }
    }
    return ioddFormat;
}


/**
* Release the memory allocated through a dataset ID.
*
* @param dsd the dataset description identifier, if <code>NULL</code> the function
*        immediately returns zero.
* @return zero for success, an error code otherwise
*/
void epr_free_dsd(EPR_SDSD* dsd)
{
    if (dsd == NULL)
        return;

    epr_free_string(dsd->ds_name);
    dsd->ds_name   = NULL;

    epr_free_string(dsd->ds_type);
    dsd->ds_type   = NULL;

    epr_free_string(dsd->filename);
    dsd->filename  = NULL;

    dsd->index = 0;
    dsd->ds_offset = 0;
    dsd->ds_size   = 0;
    dsd->num_dsr   = 0;
    dsd->dsr_size  = 0;

    free(dsd);
}


#define EPR_LENGTH_DS_NAME_IDENTIFIER      9
#define EPR_LENGTH_DS_TYPE_IDENTIFIER      8
#define EPR_LENGTH_FILENAME_IDENTIFIER    10
#define EPR_LENGTH_DS_OFFSEN_IDENTIFIER   11
#define EPR_LENGTH_DS_SIZE_IDENTIFIER      9
#define EPR_LENGTH_NUM_DSR_IDENTIFIER      9
#define EPR_LENGTH_DSR_SIZE_IDENTIFIER    10

#define EPR_LENGTH_DS_NAME_FIELD          39
#define EPR_LENGTH_DS_TYPE_FIELD          10
#define EPR_LENGTH_DS_FILENAME_FIELD      74
#define EPR_LENGTH_DS_OFFSEN_FIELD        39
#define EPR_LENGTH_DS_SIZE_FIELD          37
#define EPR_LENGTH_NUM_DSR_FIELD          20
#define EPR_LENGTH_DSR_SIZE_FIELD         28

#define EPR_LENGTH_EMPTY_FIELD            33


/**
* Reads a dataset description from an ENVISAT product file.
*
* @param envisat_source_file the handle of the given ENVISAT product file,
*        must not be <code>NULL</code>
* @param pos number of the dataset description in ENVISAT product file,
* @return a new dataset description or <code>NULL</code> if an error occured.
*/
EPR_SDSD* epr_read_each_dsd(FILE* envisat_source_file, int* pos)
{
    uint l;
    uint l_limit;
    char code_block[EPR_LINE_MAX_LENGTH];
    EPR_SDSD* dsd;
    int  ch = '"';
    char* tmp;

    if (envisat_source_file == NULL)
    {
        epr_set_err(e_err_file_access_denied,
            "epr_read_each_dsd: the product file handle must not be NULL");
        return NULL;
    }

    dsd = (EPR_SDSD*) calloc(1, sizeof (EPR_SDSD));
    if (dsd == NULL)
    {
        epr_set_err(e_err_out_of_memory,
            "epr_read_each_dsd: out of memory");
        return NULL;
    }

    if (* pos == 0)
    {
        l_limit = 9999;
    }
    else l_limit = 0;

    for (l = 0; l <= l_limit; l++)
    {
        fgets(code_block, EPR_LINE_MAX_LENGTH, envisat_source_file);
        if (strncmp (code_block, "DS_NAME=\"", EPR_LENGTH_DS_NAME_IDENTIFIER) == 0)
        {
            /* DS_NAME to be searched for */
            if (((uint)strlen(code_block) != (uint)(EPR_LENGTH_DS_NAME_FIELD)) || ((uint)(strrchr(code_block, ch) - code_block) != (uint)(strlen(code_block) - 2)))
            {
                epr_set_err(e_err_invalid_data_format,
                    "epr_read_each_dsd: invalid dataset name format");
                epr_free_dsd(dsd);
                return NULL;
            }
            dsd->ds_name = epr_sub_string(code_block, EPR_LENGTH_DS_NAME_IDENTIFIER, strlen(code_block) - EPR_LENGTH_DS_NAME_IDENTIFIER - 2);
            if (dsd->ds_name == NULL)
            {
                epr_set_err(e_err_invalid_value,
                    "epr_read_each_dsd: invalid DS_NAME value");
                epr_free_dsd(dsd);
                return NULL;
            }

            /* DS_TYPE to be searched for */
            fgets(code_block, EPR_LINE_MAX_LENGTH, envisat_source_file);
            if (strncmp (code_block, "DS_TYPE=", EPR_LENGTH_DS_TYPE_IDENTIFIER) == 0)
            {
                dsd->ds_type = epr_sub_string(code_block, EPR_LENGTH_DS_TYPE_IDENTIFIER, strlen(code_block) - EPR_LENGTH_DS_TYPE_IDENTIFIER - 1);
                if (dsd->ds_type == NULL)
                {
                    epr_set_err(e_err_invalid_value,
                        "epr_read_each_dsd: invalid DS_TYPE value");
                    epr_free_dsd(dsd);
                    return NULL;
                }
            }

            /* FILENAME to be searched for */
            fgets(code_block, EPR_LINE_MAX_LENGTH, envisat_source_file);
            if (strncmp (code_block, "FILENAME=\"", EPR_LENGTH_FILENAME_IDENTIFIER) == 0)
            {
                if (((uint)strlen(code_block) != (uint)(EPR_LENGTH_DS_FILENAME_FIELD)) || ((uint)(strrchr(code_block, ch) - code_block) != (uint)(strlen(code_block) - 2)))
                {
                    epr_set_err(e_err_invalid_data_format,
                        "epr_read_each_dsd: invalid dataset filename format");
                    epr_free_dsd(dsd);
                    return NULL;
                }
                dsd->filename = epr_sub_string(code_block, EPR_LENGTH_FILENAME_IDENTIFIER, strlen(code_block) - EPR_LENGTH_FILENAME_IDENTIFIER - 1);
                if (dsd->ds_name == NULL)
                {
                    epr_set_err(e_err_invalid_data_format,
                        "epr_read_each_dsd: invalid file name");
                    epr_free_dsd(dsd);
                    return NULL;
                }
            }

            /* DS_OFFSET to be searched for */
            fgets(code_block, EPR_LINE_MAX_LENGTH, envisat_source_file);
            if (strncmp (code_block, "DS_OFFSET=+", EPR_LENGTH_DS_OFFSEN_IDENTIFIER) == 0)
            {
                if (((uint)strlen(code_block) != (uint)(EPR_LENGTH_DS_OFFSEN_FIELD)) || (strncmp(code_block + strlen(code_block) - strlen("<bytes>") - 1, "<bytes>", strlen("<bytes>")) != 0))
                {
                    epr_set_err(e_err_invalid_data_format,
                        "epr_read_each_dsd: invalid dataset filename format");
                    epr_free_dsd(dsd);
                    return NULL;
                }

                tmp = epr_sub_string(code_block, EPR_LENGTH_DS_OFFSEN_IDENTIFIER, strlen(code_block) - strlen("<bytes>") - EPR_LENGTH_DS_OFFSEN_IDENTIFIER - 1);
                dsd->ds_offset = (uint)epr_str_to_number(tmp);
                epr_free_string(tmp);
                if (dsd->ds_offset == -1)
                {
                    epr_set_err(e_err_invalid_data_format,
                        "epr_read_each_dsd: invalid OFFSET value");
                    epr_free_dsd(dsd);
                    return NULL;
                }
            }

            /* DS_SIZE to be searched for */
            fgets(code_block, EPR_LINE_MAX_LENGTH, envisat_source_file);
            if (strncmp (code_block, "DS_SIZE=+", EPR_LENGTH_DS_SIZE_IDENTIFIER) == 0)
            {
                if (((uint)strlen(code_block) != (uint)(EPR_LENGTH_DS_SIZE_FIELD)) || (strncmp(code_block + strlen(code_block) - strlen("<bytes>") - 1, "<bytes>", strlen("<bytes>")) != 0))
                {
                    epr_set_err(e_err_invalid_data_format,
                        "epr_read_each_dsd: invalid dataset filename format");
                    epr_free_dsd(dsd);
                    return NULL;
                }
                tmp = epr_sub_string(code_block, EPR_LENGTH_DS_SIZE_IDENTIFIER, strlen(code_block) - strlen("<bytes>") - EPR_LENGTH_DS_SIZE_IDENTIFIER - 1);
                dsd->ds_size = (uint)epr_str_to_number(tmp);
                epr_free_string(tmp);
                if (dsd->ds_size == -1)
                {
                    epr_set_err(e_err_invalid_data_format,
                        "epr_read_each_dsd: invalid OFFSET value");
                    epr_free_dsd(dsd);
                    return NULL;
                }
            }

            /* NUM_DSR to be searched for */
            fgets(code_block, EPR_LINE_MAX_LENGTH, envisat_source_file);
            if (strncmp (code_block, "NUM_DSR=+", EPR_LENGTH_NUM_DSR_IDENTIFIER) == 0)
            {
                if ((uint)strlen(code_block) != (uint)(EPR_LENGTH_NUM_DSR_FIELD))
                {
                    epr_set_err(e_err_invalid_data_format,
                        "epr_read_each_dsd: invalid dataset record number format");
                    epr_free_dsd(dsd);
                    return NULL;
                }
                tmp = epr_sub_string(code_block, EPR_LENGTH_NUM_DSR_IDENTIFIER, strlen(code_block) - EPR_LENGTH_NUM_DSR_IDENTIFIER - 1);
                dsd->num_dsr = (uint)epr_str_to_number(tmp);
                epr_free_string(tmp);
                if (dsd->num_dsr == -1)
                {
                    epr_set_err(e_err_invalid_data_format,
                        "epr_read_each_dsd: invalid dsr number value");
                    epr_free_dsd(dsd);
                    return NULL;
                }
            }

            /* DSR_SIZE to be searched for */
            fgets(code_block, EPR_LINE_MAX_LENGTH, envisat_source_file);
            if (strncmp (code_block, "DSR_SIZE=+", EPR_LENGTH_DSR_SIZE_IDENTIFIER) == 0)
            {
                if (((uint)strlen(code_block) != (uint)(EPR_LENGTH_DSR_SIZE_FIELD)) || (strncmp(code_block + strlen(code_block) - strlen("<bytes>") - 1, "<bytes>", strlen("<bytes>")) != 0))
                {
                    epr_set_err(e_err_invalid_data_format,
                        "epr_read_each_dsd: invalid dataset record size format");
                    epr_free_dsd(dsd);
                    return NULL;
                }
                tmp = epr_sub_string(code_block, EPR_LENGTH_DSR_SIZE_IDENTIFIER, strlen(code_block) - strlen("<bytes>") - EPR_LENGTH_DSR_SIZE_IDENTIFIER - 1);
                dsd->dsr_size = (uint)epr_str_to_number(tmp);
                epr_free_string(tmp);
                if (dsd->dsr_size == -1)
                {
                    epr_set_err(e_err_invalid_data_format,
                        "epr_read_each_dsd: invalid record size value");
                    epr_free_dsd(dsd);
                    return NULL;
                }
            }

            /* EMPTY LINE BETWEEN DSD's */
            fgets(code_block, EPR_LINE_MAX_LENGTH, envisat_source_file);
            if ((strlen(code_block) > 0) && (code_block[0] != ' '))
            {
                epr_set_err(e_err_invalid_data_format,
                    "epr_read_each_dsd: invalid code_block, must be empty");
                epr_free_dsd(dsd);
                return NULL;
            }
            *pos = *pos + 1;
            return dsd;
        }
    }
    epr_free_dsd(dsd);
    return NULL;
}


/**
* Finds the first dataset description from an ENVISAT product file.
*
* @param envisat_source_file the handle of the given ENVISAT product file,
*        must not be <code>NULL</code>
* @param sph_length [bytes] the length of SPH part from the given ENVISAT product file,
*        must not be <code>NULL</code>
* @return the offset to first searched for dsd or <code>0</code> if not found.
*/
uint epr_find_first_dsd(FILE* envisat_source_file, uint sph_length)
{
    uint l;
    char code_block[EPR_LINE_MAX_LENGTH];

    if (envisat_source_file == NULL)
    {
        epr_set_err(e_err_file_access_denied,
            "epr_find_first_dsd: the product file handle must not be NULL");
        return 0;
    }

    l = 0;
    while (l < sph_length)
    {
        fgets(code_block, EPR_LINE_MAX_LENGTH, envisat_source_file);
        if (strncmp (code_block, "DS_NAME=\"", EPR_LENGTH_DS_NAME_IDENTIFIER) == 0)
        {
            return l;
        }
        else {
            l += strlen(code_block);
        }
    }
    return 0;
}

#define EPR_LENGTH_NUM_DSD_FIELD        20

/**
* Reads all dataset descriptions from an ENVISAT product file.
*
* @param product_id the file identifier, if <code>NULL</code> the function
*        immediately returns <code>NULL</code>.
* @return an array of dataset descriptions or <code>NULL</code> if an error occured.
*/
EPR_SPtrArray* epr_read_all_dsds(EPR_SProductId* product_id)
{
    EPR_SPtrArray* dsds_array = NULL;
    EPR_SRecord* dsd_record = NULL;
    const EPR_SField* field;
    EPR_SDSD* dsd;
    uint sph_length;
    char* code_block;
    int numread;
    uint dsd_number = 0;
    uint dsd_begin = 0;
    uint dsd_index;

    if (product_id == NULL) {
        epr_set_err(e_err_null_pointer,
            "epr_read_all_dsds: product_id must not be NULL");
        return NULL;
    }

    dsds_array = epr_create_ptr_array(32);

    field = epr_get_field(product_id->mph_record, "NUM_DSD");
    dsd_number = ((uint*) field->elems)[0];
    /*dsd_begin = epr_api.epr_head_size - dsd_number * EPR_DSD_SIZE;*/

    if (fseek(product_id->istream, EPR_MPH_SIZE, SEEK_SET) != 0) {
        epr_set_err(e_err_file_access_denied,
            "epr_read_all_dsds: file seek failed");
        if (dsds_array != NULL) {
            epr_free_ptr_array(dsds_array);
        }
        return NULL;
    }
    field = epr_get_field(product_id->mph_record, "SPH_SIZE");
    sph_length = ((uint*) field->elems)[0];
    dsd_begin = EPR_MPH_SIZE + (uint)epr_find_first_dsd(product_id->istream, sph_length);
    if (dsd_begin == EPR_MPH_SIZE) {
        epr_set_err(e_err_file_access_denied,
            "epr_read_all_dsds: no DS_NAME in SPH");
        if (dsds_array != NULL) {
            epr_free_ptr_array(dsds_array);
        }
        return NULL;
    }
    for(dsd_index = 0; dsd_index < dsd_number; dsd_index ++) {
        if (fseek(product_id->istream, dsd_begin + dsd_index * EPR_DSD_SIZE, SEEK_SET) != 0) {
            epr_set_err(e_err_file_access_denied,
                "epr_read_all_dsds: file seek failed");
            if (dsds_array != NULL) {
                epr_free_ptr_array(dsds_array);
            }
            return NULL;
        }
        code_block = epr_create_string(EPR_DSD_SIZE);
        numread = fread(code_block, 1, EPR_DSD_SIZE, product_id->istream);
        if ((uint)numread != EPR_DSD_SIZE) {
            epr_set_err(e_err_file_read_error,
                "epr_read_all_dsds: error in reading SPH from product data file");
            if (code_block != NULL) {
                epr_free_string(code_block);
                code_block = NULL;
            }
            if (dsds_array != NULL) {
                epr_free_ptr_array(dsds_array);
            }
            return NULL;
        }


        /* If this is NOT an empty DSD (empty DSD's seem to be quite 'normal')
        */
        if ((strlen(code_block) > 0) && (code_block[0] != ' ')) {

            dsd_record = epr_parse_header("dsd", code_block);

            dsd = epr_create_dsd(dsd_index);

            if (dsd == NULL) {
                epr_set_err(e_err_out_of_memory,
                    "epr_read_all_dsds: out of memory");
                if (code_block != NULL) {
                    epr_free_string(code_block);
                    code_block = NULL;
                }
                if (dsd_record != NULL) {
                    epr_free_record_info(dsd_record->info);
                    dsd_record->info = NULL;
                    epr_free_record(dsd_record);
                    dsd_record = NULL;
                }
                if (dsds_array != NULL) {
                    epr_free_ptr_array(dsds_array);
                }
                return NULL;
            }

            field = epr_get_field(dsd_record, "DS_NAME");
            dsd->ds_name = epr_clone_string((char*)field->elems);
            field = epr_get_field(dsd_record, "DS_TYPE");
            dsd->ds_type = epr_sub_string((char*)field->elems, 0, sizeof(uchar));
            field = epr_get_field(dsd_record, "FILENAME");
            dsd->filename = epr_clone_string((char*)field->elems);
            field = epr_get_field(dsd_record, "DS_OFFSET");
            dsd->ds_offset = (uint)((uint*) field->elems)[0];

            field = epr_get_field(dsd_record, "DS_SIZE");
            dsd->ds_size = (uint)((uint*) field->elems)[0];

            field = epr_get_field(dsd_record, "NUM_DSR");
            dsd->num_dsr = (uint)((uint*) field->elems)[0];

            field = epr_get_field(dsd_record, "DSR_SIZE");
            dsd->dsr_size = (uint)((uint*) field->elems)[0];

            epr_add_ptr_array_elem(dsds_array, dsd);

            if (dsd_record != NULL) {
            /* NOTE:dsd_record->info is not a shared object, it is NOT used by
            * multiple instances of a DSD record, and thus, we free it here!
                */
                epr_free_record_info(dsd_record->info);
                dsd_record->info = NULL;

                epr_free_record(dsd_record);
                dsd_record = NULL;
            } else {
                printf("%s\n", epr_get_last_err_message());
            }

        } else {
            epr_log(e_log_debug, "empty DSD seen (don't worry)");
        }

        epr_free_string(code_block);
        code_block = NULL;
    }

    rewind(product_id->istream);
    return dsds_array;
}

uint epr_get_num_dsds(const EPR_SProductId* product_id)
{
    return product_id->dsd_array->length;
}

EPR_SDSD* epr_get_dsd_at(const EPR_SProductId* product_id, uint dsd_index)
{
    EPR_SDSD* dsd = NULL;
    dsd = (EPR_SDSD*) epr_get_ptr_array_elem_at(product_id->dsd_array, dsd_index);
    return dsd;
}
