/*
 * $Id: epr_field.c,v 1.1.1.1 2004-10-28 19:22:23 norman Exp $
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


/**
 * Creates the field information of the given record and returns the poiter at it.
 *
 * @param data_type_id the data type identifier
 * @param description the field description
 * @param field_name the field name
 * @param num_elems the number of field elements
 * @param num_bytes the number of bytes in each element
 * @param more_count the number of the element repetition
 * @param unit the unit descrimtion (name)
 *
 * @return the the pointer at the field information, or <code>NULL</code> if the file
 *         it is not enough memory for some field_info element.
 */
EPR_SFieldInfo* epr_create_field_info(EPR_EDataTypeId data_type_id, char* description, char* field_name, uint num_elems, uint num_bytes, uint more_count, char* unit)
{
    EPR_SFieldInfo* field_info = NULL;
    uint data_type_size;

    field_info = (EPR_SFieldInfo*) calloc(1, sizeof (EPR_SFieldInfo));
    if (field_info == NULL)
    {
        epr_set_err(e_err_out_of_memory,
                    "epr_create_field_info: out of memory");
        return NULL;
    }

    epr_assign_string(&field_info->name, field_name);
    if (field_info->name == NULL)
    {
        epr_set_err(e_err_out_of_memory,
                    "epr_create_field_info: out of memory");
        return NULL;
    }

    if (description == NULL)
    {
        field_info->description = NULL;
    } else
        epr_assign_string(&field_info->description, description);

    if (unit == NULL)
    {
        field_info->unit = NULL;
    } else
        epr_assign_string(&field_info->unit, unit);

    data_type_size = epr_get_data_type_size(data_type_id);

    field_info->num_elems = num_elems;


    /* IMPORTANT !!!!it's one or the other !!!!*/
    field_info->tot_size = num_elems * num_bytes * more_count;
    /*for the case data_type_size=num_bytes*/
    /*this is not true for the case "spare" (1!=60)*/
    /*    field_info->tot_size = num_elems * data_type_size * more_count;*/


    field_info->data_type_id = data_type_id;

    return field_info;
}


/**
 * Frees the memory allocated by the given 'field_info'.
 *
 * <p> After calling this function the give field_info pointer
 * should not be used anymore.
 *
 * @param field_info the field_info to be released
 */
void epr_free_field_info(EPR_SFieldInfo* field_info)
{
    if (field_info == NULL)
        return;

    epr_free_string(field_info->name);
    field_info->name = NULL;

    epr_free_string(field_info->description);
    field_info->description = NULL;

    epr_free_string(field_info->unit);
    field_info->unit = NULL;

    field_info->num_elems = 0;
    field_info->tot_size = 0;
    field_info->data_type_id = 0;

    free(field_info);
}


/*
   Function:    epr_create_field
   Access:      public API
   Changelog:   2002/01/23  mp initial version
 */
/**
 * Creates a new field instance belongs to the given
 * record_info.
 *
 * @param the pointer at the field information.
 *        must not be <code>NULL</code>
 * @return the new field instance
 *         or <code>NULL</code> if an error occured.
 */
EPR_SField* epr_create_field(EPR_SFieldInfo* field_info)
{
    EPR_SField* field = NULL;
    uint data_type_size = 0;

    field = (EPR_SField*) calloc(1, sizeof(EPR_SField));
    if (field == NULL)
    {
        epr_set_err(e_err_out_of_memory,
                    "epr_create_field: out of memory");
        return NULL;
    }
    field->magic = EPR_MAGIC_FIELD;


    data_type_size = epr_get_data_type_size(field_info->data_type_id);
    if (data_type_size == 0)
    {
        epr_set_err(e_err_illegal_data_type,
                    "epr_create_field: illegal field_info data type identifier");
        return NULL;
    }
    field->info = field_info;

    if (field_info->num_elems == 0)
    {
        free(field);
        epr_set_err(e_err_invalid_value,
                    "epr_create_field: field_info->num_elems is zero");
        return NULL;
    }

    if (field_info->data_type_id == e_tid_spare) {
        field->elems = calloc(field_info->tot_size, 1 /*byte*/);
    } else {
        if (field_info->data_type_id == e_tid_string) {
			/*
			  Note that string always are considered as one single element,
			  so we get the total number of characters from tot_size.
			  Addidionally we reserve an extra character for the trailing zero (terminator),
			*/
            field->elems = calloc(field_info->tot_size + 1, sizeof (char));
        } else {
            field->elems = calloc(field_info->num_elems, data_type_size);
        }
    }
    if (field->elems == NULL)
    {
        free(field);
        epr_set_err(e_err_out_of_memory,
                    "epr_create_field: out of memory");
        return NULL;
    }

/*
	if (field->info->data_type_id == e_tid_string) {
		printf("string field: name=%s, num_elems=%d, tot_size=%d\n", field->info->name, field->info->num_elems, field->info->tot_size);
	}
*/
    return field;
}


/**
 * Frees the memory allocated by the given 'field'.
 *
 * <p> After calling this function the give field pointer
 * should not be used anymore.
 *
 * @param field the field to be released
 */
void epr_free_field(EPR_SField* field)
{
    if (field == NULL)
        return;

    /* Do NOT free field->info since many fields can
       share the same field->info! */
    field->info = NULL;

    if (field->elems != NULL)
    {
        free(field->elems);
        field->elems = NULL;
    }

    free(field);
}

/**
 * Gets a full field from the given record.
 *
 * <p> The field is hier identified through the given name.
 * <br>It contains the field info and all corresponding values.
 *
 * @param record the record identifier, must not be <code>NULL</code>
 * @param field_name the the name of required field, must not be <code>NULL</code>.
 *
 * @return the field or <code>NULL</code> if an error occured.
 */
const EPR_SField* epr_get_field(const EPR_SRecord* record, const char* field_name)
{
    EPR_SField* field;
    uint field_index;

    epr_clear_err();

    if (record == NULL) {
            epr_set_err(e_err_invalid_record_name,
                "epr_get_field: record must not be NULL");
            return NULL;
    }
    if (field_name == NULL) {
            epr_set_err(e_err_invalid_record_name,
                "epr_get_field: field_name must not be NULL");
            return NULL;
    }

    for (field_index = 0; field_index < record->num_fields; field_index++)  {
        field = record->fields[field_index];
        if (strcmp(field_name, field->info->name) == 0) {
            return field;
        }
    }
    epr_set_err(e_err_illegal_arg,
                    "epr_get_field: field not found");
    return NULL;
}


uint epr_get_field_num_elems(const EPR_SField* field)
{
    epr_clear_err();

    if (field == NULL) {
        epr_set_err(e_err_null_pointer,
                    "epr_get_field_num_elems: field-identifier must not be NULL");
        return 0;
    }
    return field->info->num_elems;
}

const char* epr_get_field_name(const EPR_SField* field)
{
    epr_clear_err();

    if (field == NULL) {
        epr_set_err(e_err_null_pointer,
                    "epr_get_field_name: field-identifier must not be NULL");
        return NULL;
    }
    return field->info->name;
}


EPR_EDataTypeId epr_get_field_type(const EPR_SField* field)
{
    epr_clear_err();

    if (field == NULL) {
        epr_set_err(e_err_null_pointer,
                    "epr_get_field_type: field-identifier must not be NULL");
        return e_tid_unknown;
    }
    return field->info->data_type_id;
}

const char* epr_get_field_unit(const EPR_SField* field)
{
    epr_clear_err();

    if (field == NULL) {
        epr_set_err(e_err_null_pointer,
                    "epr_get_field_unit: field-identifier must not be NULL");
        return NULL;
    }
    return field->info->unit;
}

const char* epr_get_field_description(const EPR_SField* field)
{
    epr_clear_err();

    if (field == NULL) {
        epr_set_err(e_err_null_pointer,
                    "epr_get_field_description: field-identifier must not be NULL");
        return NULL;
    }
    return field->info->description;
}
