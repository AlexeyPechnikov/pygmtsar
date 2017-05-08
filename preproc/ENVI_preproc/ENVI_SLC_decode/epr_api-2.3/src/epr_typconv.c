/*
 * $Id: epr_typconv.c,v 1.2 2009-03-27 10:25:54 sabine Exp $
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
#include "epr_field.h"

/*********************************** TYPE CONVERSION ***********************************/


/**
 * Interpretes a memory as a <code>char</code> value
 *
 * @param field the pointer at the array to convert
 * @param elem_index the index of the element
 * in the given array to convert
 *
 * @return the <code>char</code> typed element
 *         or <code>error_code</code> if an error occured.
 */
char epr_get_field_elem_as_char(const EPR_SField* field, uint elem_index)
{
    epr_clear_err();

    if (field == NULL) {
        epr_set_err(e_err_invalid_field_name,
                    "epr_get_field_elem_as_char: invalid field name");
        return (char)0;
    }
    if (elem_index >= field->info->num_elems) {
        epr_set_err(e_err_invalid_value,
                    "epr_get_field_elem_as_char: invalid elem_index parameter");
        return (char)0;
    }
    if (field->info->data_type_id != e_tid_char) {
        epr_set_err(e_err_invalid_data_format,
                    "epr_get_field_elem_as_char: invalid type");
        return (char)0;
    }
    return ((char*) field->elems)[elem_index];
}

/**
 * Interpretes a memory data as a <code>char</code> data
 *
 * @param field the pointer at the array to convert
 *
 * @return the <code>char</code> typed element
 *         or <code>NULL</code> if an error occured.
 */
const char* epr_get_field_elems_char(const EPR_SField* field)
{
    epr_clear_err();

    if (field == NULL) {
        epr_set_err(e_err_invalid_field_name,
                    "epr_get_field_elems_as_chars: invalid field name");
        return NULL;
    }
    if (field->info->data_type_id != e_tid_char) {
        epr_set_err(e_err_invalid_data_format,
                    "epr_get_field_elems_as_chars: invalid type");
        return NULL;
    }
    return (char*) field->elems;
}

/**
 * Interpretes a memory as a <code>uchar</code> value
 *
 * @param field the pointer at the array to convert
 * @param elem_index the index of the element
 * in the given array to convert
 *
 * @return the <code>uchar</code> typed element
 *         or <code>error_code</code> if an error occured.
 */
uchar epr_get_field_elem_as_uchar(const EPR_SField* field, uint elem_index)
{
    epr_clear_err();

    if (field == NULL) {
        epr_set_err(e_err_invalid_field_name,
                    "epr_get_field_elem_as_uchar: invalid field name");
        return (uchar)0;
    }
    if (elem_index >= field->info->num_elems) {
        epr_set_err(e_err_invalid_value,
                    "epr_get_field_elem_as_uchar: invalid elem_index parameter");
        return (uchar)0;
    }
    if (field->info->data_type_id != e_tid_uchar) {
        epr_set_err(e_err_invalid_data_format,
                    "epr_get_field_elem_as_uchar: invalid type");
        return (uchar)0;
    }
    return ((uchar*) field->elems)[elem_index];
}

/**
 * Interpretes a memory data as a <code>uchar</code> data
 *
 * @param field the pointer at the array to convert
 *
 * @return the <code>uchar</code> typed element
 *         or <code>NULL</code> if an error occured.
 */
const uchar* epr_get_field_elems_uchar(const EPR_SField* field)
{
    epr_clear_err();

    if (field == NULL) {
        epr_set_err(e_err_invalid_field_name,
                    "epr_get_field_elems_as_uchars: invalid field name");
        return NULL;
    }
    if (field->info->data_type_id != e_tid_uchar) {
        epr_set_err(e_err_invalid_data_format,
                    "epr_get_field_elems_as_uchars: invalid type");
        return NULL;
    }
    return (uchar*) field->elems;
}

/**
 * Interpretes a memory as a <code>short</code> value
 *
 * @param field the pointer at the array to convert
 * @param elem_index the index of the element
 * in the given array to convert
 *
 * @return the <code>short</code> typed element
 *         or <code>error_code</code> if an error occured.
 */
short epr_get_field_elem_as_short(const EPR_SField* field, uint elem_index)
{
    epr_clear_err();

    if (field == NULL) {
        epr_set_err(e_err_invalid_field_name,
                    "epr_get_field_elem_as_short: invalid field name");
        return (short)0;
    }
    if (elem_index >= field->info->num_elems) {
        epr_set_err(e_err_invalid_value,
                    "epr_get_field_elem_as_short: invalid elem_index parameter");
        return (short)0;
    }


    if (field->info->data_type_id == e_tid_uchar) {
        return (short)((uchar*) field->elems)[elem_index];
    }
    if (field->info->data_type_id == e_tid_char) {
        return (short)((char*) field->elems)[elem_index];
    }
    if (field->info->data_type_id == e_tid_short) {
        return (short)((short*) field->elems)[elem_index];
    }

    epr_set_err(e_err_invalid_data_format,
                "epr_get_field_elem_as_short: invalid type");

    return (short)0;
}

/**
 * Interpretes a memory data as a <code>short</code> data
 *
 * @param field the pointer at the array to convert
 *
 * @return the <code>short</code> typed element
 *         or <code>NULL</code> if an error occured.
 */
const short* epr_get_field_elems_short(const EPR_SField* field)
{
    epr_clear_err();

    if (field == NULL) {
        epr_set_err(e_err_invalid_field_name,
                    "epr_get_field_elems_as_shorts: invalid field name");
        return NULL;
    }
    if (field->info->data_type_id != e_tid_short) {
        epr_set_err(e_err_invalid_data_format,
                    "epr_get_field_elems_as_shorts: invalid type");
        return NULL;
    }
    return (short*) field->elems;
}

/**
 * Interpretes a memory as a <code>ushort</code> value
 *
 * @param field the pointer at the array to convert
 * @param elem_index the index of the element
 * in the given array to convert
 *
 * @return the <code>ushort</code> typed element
 *         or <code>error_code</code> if an error occured.
 */
ushort epr_get_field_elem_as_ushort(const EPR_SField* field, uint elem_index)
{
    epr_clear_err();

    if (field == NULL) {
        epr_set_err(e_err_invalid_field_name,
                    "epr_get_field_elem_as_ushort: invalid field name");
        return (ushort)0;
    }
    if (elem_index >= field->info->num_elems) {
        epr_set_err(e_err_invalid_value,
                    "epr_get_field_elem_as_ushort: invalid elem_index parameter");
        return (ushort)0;
    }

    if (field->info->data_type_id == e_tid_uchar) {
        return (ushort)((uchar*) field->elems)[elem_index];
    }
    if (field->info->data_type_id == e_tid_char) {
        return (ushort)((char*) field->elems)[elem_index];
    }
    if (field->info->data_type_id == e_tid_ushort) {
        return (ushort)((ushort*) field->elems)[elem_index];
    }

    epr_set_err(e_err_invalid_data_format,
                "epr_get_field_elem_as_ushort: invalid type");

    return (ushort)0;
}

/**
 * Interpretes a memory data as a <code>ushort</code> data
 *
 * @param field the pointer at the array to convert
 *
 * @return the <code>ushort</code> typed element
 *         or <code>NULL</code> if an error occured.
 */
const ushort* epr_get_field_elems_ushort(const EPR_SField* field)
{
    epr_clear_err();

    if (field == NULL) {
        epr_set_err(e_err_invalid_field_name,
                    "epr_get_field_elems_as_ushorts: invalid field name");
        return NULL;
    }
    if (field->info->data_type_id != e_tid_ushort) {
        epr_set_err(e_err_invalid_data_format,
                    "epr_get_field_elems_as_ushorts: invalid type");
        return NULL;
    }
    return (ushort*) field->elems;
}

/**
 * Interpretes a memory as a <code>int</code> value.
 *
 * <p> If an error occurs the method returns <code>0</code> (zero).
 * Whether an error really occured when zero is returned can by determined by
 * using the <code>epr_get_last_err_code</code> function.
 *
 * @param field the pointer at the array to convert
 * @param elem_index the index of the element
 * in the given array to convert
 *
 * @return a <code>int</code> value
 */
int epr_get_field_elem_as_int(const EPR_SField* field, uint elem_index)
{
    epr_clear_err();

    if (field == NULL) {
        epr_set_err(e_err_invalid_field_name,
                    "epr_get_field_elem_as_long: invalid field name");
        return (int)0;
    }
    if (elem_index >= field->info->num_elems) {
        epr_set_err(e_err_invalid_value,
                    "epr_get_field_elem_as_long: invalid elem_index parameter");
        return (int)0;
    }

    if (field->info->data_type_id == e_tid_uchar) {
        return (int)((uchar*) field->elems)[elem_index];
    }
    if (field->info->data_type_id == e_tid_char) {
        return (int)((char*) field->elems)[elem_index];
    }
    if (field->info->data_type_id == e_tid_ushort) {
        return (int)((ushort*) field->elems)[elem_index];
    }
    if (field->info->data_type_id == e_tid_short) {
        return (int)((short*) field->elems)[elem_index];
    }
    if (field->info->data_type_id == e_tid_int) {
        return (int)((int*) field->elems)[elem_index];
    }

    epr_set_err(e_err_invalid_data_format,
                "epr_get_field_elem_as_long: invalid type");

    return (int)0;
}

/**
 * Interpretes a memory data as a <code>int</code> data
 *
 * @param field the pointer at the array to convert
 *
 * @return the <code>int</code> typed element
 *         or <code>NULL</code> if an error occured.
 */
const int* epr_get_field_elems_int(const EPR_SField* field)
{
    epr_clear_err();

    if (field == NULL) {
        epr_set_err(e_err_invalid_field_name,
                    "epr_get_field_elems_as_longs: invalid field name");
        return NULL;
    }
    if (field->info->data_type_id != e_tid_int) {
        epr_set_err(e_err_invalid_data_format,
                    "epr_get_field_elems_as_longs: invalid type");
        return NULL;
    }
    return (int*) field->elems;
}

/**
 * Interpretes a memory as a <code>uint</code> value
 *
 * @param field the pointer at the array to convert
 * @param elem_index the index of the element
 * in the given array to convert
 *
 * @return the <code>uint</code> typed element
 *         or <code>error_code</code> if an error occured.
 */
uint epr_get_field_elem_as_uint(const EPR_SField* field, uint elem_index)
{
    epr_clear_err();

    if (field == NULL) {
        epr_set_err(e_err_invalid_field_name,
                    "epr_get_field_elem_as_uint: invalid field name");
        return (uint)0;
    }
    if (elem_index >= field->info->num_elems) {
        epr_set_err(e_err_invalid_value,
                    "epr_get_field_elem_as_uint: invalid elem_index parameter");
        return (uint)0;
    }

    if (field->info->data_type_id == e_tid_uint) {
        return (uint)((uint*) field->elems)[elem_index];
    }
    if (field->info->data_type_id == e_tid_uchar) {
        return (uint)((uchar*) field->elems)[elem_index];
    }
    if (field->info->data_type_id == e_tid_char) {
        return (uint)((char*) field->elems)[elem_index];
    }
    if (field->info->data_type_id == e_tid_ushort) {
        return (uint)((ushort*) field->elems)[elem_index];
    }
    if (field->info->data_type_id == e_tid_short) {
        return (uint)((short*) field->elems)[elem_index];
    }

    epr_set_err(e_err_invalid_data_format,
                "epr_get_field_elem_as_uint: invalid type");

    return (uint)0;
}

/**
 * Interpretes a memory data as a <code>uint</code> data
 *
 * @param field the pointer at the array to convert
 *
 * @return the <code>uint</code> typed element
 *         or <code>NULL</code> if an error occured.
 */
const uint* epr_get_field_elems_uint(const EPR_SField* field)
{
    epr_clear_err();

    if (field == NULL) {
        epr_set_err(e_err_invalid_field_name,
                    "epr_get_field_elems_as_uints: invalid field name");
        return NULL;
    }
    if (field->info->data_type_id != e_tid_uint) {
        epr_set_err(e_err_invalid_data_format,
                    "epr_get_field_elems_as_uints: invalid type");
        return NULL;
    }
    return (uint*) field->elems;
}

/**
 * Interpretes a memory as a <code>float</code> value
 *
 * @param field the pointer at the array to convert
 * @param elem_index the index of the element
 * in the given array to convert
 *
 * @return the <code>float</code> typed element
 *         or <code>error_code</code> if an error occured.
 */
float epr_get_field_elem_as_float(const EPR_SField* field, uint elem_index)
{
    epr_clear_err();

    if (field == NULL) {
        epr_set_err(e_err_invalid_field_name,
                    "epr_get_field_elem_as_float: invalid field name");
        return 0.0;

    }
    if (elem_index >= field->info->num_elems) {
        epr_set_err(e_err_invalid_value,
                    "epr_get_field_elem_as_float: invalid elem_index parameter");
        return 0.0;
    }

    if (field->info->data_type_id == e_tid_float) {
        return (float)((float*) field->elems)[elem_index];
    }
    if (field->info->data_type_id == e_tid_uchar) {
        return (float)((uchar*) field->elems)[elem_index];
    }
    if (field->info->data_type_id == e_tid_char) {
        return (float)((char*) field->elems)[elem_index];
    }
    if (field->info->data_type_id == e_tid_ushort) {
        return (float)((ushort*) field->elems)[elem_index];
    }
    if (field->info->data_type_id == e_tid_short) {
        return (float)((short*) field->elems)[elem_index];
    }
    if (field->info->data_type_id == e_tid_uint) {
        return (float)((uint*) field->elems)[elem_index];
    }
    if (field->info->data_type_id == e_tid_int) {
        return (float)((int*) field->elems)[elem_index];
    }

    epr_set_err(e_err_invalid_data_format,
                "epr_get_field_elems_as_float: invalid type");

    return 0.0;
}

/**
 * Interpretes a memory data as a <code>float</code> data
 *
 * @param field the pointer at the array to convert
 *
 * @return the <code>float</code> typed element
 *         or <code>NULL</code> if an error occured.
 */
const float* epr_get_field_elems_float(const EPR_SField* field)
{
    epr_clear_err();

    if (field == NULL) {
        epr_set_err(e_err_invalid_field_name,
                    "epr_get_field_elem_as_floats: invalid field name");
        return NULL;
    }
    if (field->info->data_type_id != e_tid_float) {
        epr_set_err(e_err_invalid_data_format,
                    "epr_get_field_elems_as_floats: invalid type");
        return NULL;
    }
    return (float*) field->elems;
}

/**
 * Interpretes a memory as a <code>double</code> value
 *
 * @param field the pointer at the array to convert
 * @param elem_index the index of the element
 * in the given array to convert
 *
 * @return the <code>double</code> typed element
 *         or <code>error_code</code> if an error occured.
 */
double epr_get_field_elem_as_double(const EPR_SField* field, uint elem_index)
{
    epr_clear_err();

    if (field == NULL) {
        epr_set_err(e_err_invalid_field_name,
                    "epr_get_field_elems_as_double: invalid field name");
        return 0.0;

    }
    if (elem_index >= field->info->num_elems) {
        epr_set_err(e_err_invalid_value,
                    "epr_get_field_elems_as_double: invalid elem_index parameter");
        return 0.0;
    }

    if (field->info->data_type_id == e_tid_double) {
        return (double)((double*) field->elems)[elem_index];
    }
    if (field->info->data_type_id == e_tid_float) {
        return (double)((float*) field->elems)[elem_index];
    }
    if (field->info->data_type_id == e_tid_uchar) {
        return (double)((uchar*) field->elems)[elem_index];
    }
    if (field->info->data_type_id == e_tid_char) {
        return (double)((char*) field->elems)[elem_index];
    }
    if (field->info->data_type_id == e_tid_ushort) {
        return (double)((ushort*) field->elems)[elem_index];
    }
    if (field->info->data_type_id == e_tid_short) {
        return (double)((short*) field->elems)[elem_index];
    }
    if (field->info->data_type_id == e_tid_uint) {
        return (double)((uint*) field->elems)[elem_index];
    }
    if (field->info->data_type_id == e_tid_int) {
        return (double)((int*) field->elems)[elem_index];
    }

    epr_set_err(e_err_invalid_data_format,
                "epr_get_field_elems_as_double: invalid type");

    return 0.0;
}

/**
 * Interpretes a memory data as a <code>double</code> data
 *
 * @param field the pointer at the array to convert
 *
 * @return the <code>double</code> typed element
 *         or <code>NULL</code> if an error occured.
 */
const double* epr_get_field_elems_double(const EPR_SField* field)
{
    epr_clear_err();

    epr_clear_err();
    if (field == NULL) {
        epr_set_err(e_err_invalid_field_name,
                    "epr_get_field_elems_as_doubles: invalid field name");
        return NULL;
    }
    if (field->info->data_type_id != e_tid_double) {
        epr_set_err(e_err_invalid_data_format,
                    "epr_get_field_elems_as_doubles: invalid type");
        return NULL;
    }
    return (double*) field->elems;
}

/**
 * Interpretes a memory data as a <code>short</code> data
 *
 * @param field the pointer at the array to convert
 * @param time the pointer at the time structure to get
 *
 * @return the time [days, seconds, microseconds]
 *         or <code>NULL</code> if an error occured.
 */
const EPR_STime* epr_get_field_elem_as_mjd(const EPR_SField* field)
{
    epr_clear_err();

    if (field == NULL) {
        epr_set_err(e_err_invalid_field_name,
                    "epr_get_field_elem_as_mjd: invalid field name");
        return NULL;

    }

    if (field->info->data_type_id != e_tid_time) {
        epr_set_err(e_err_invalid_data_format,
                    "epr_get_field_elem_as_mjd: invalid type");
        return NULL;
    }

	return (const EPR_STime*) field->elems;
}

/**
 * Interpretes a memory data as a string.
 *
 * @param field the pointer at the array to convert
 *
 * @return the <code>char</code> typed element
 *         or <code>NULL</code> if an error occured.
 */
const char* epr_get_field_elem_as_str(const EPR_SField* field)
{
    epr_clear_err();

    if (field == NULL) {
        epr_set_err(e_err_invalid_field_name,
                    "epr_get_field_elem_as_str: invalid field name");
        return NULL;
    }
    if (field->info->data_type_id != e_tid_string) {
        epr_set_err(e_err_invalid_data_format,
                    "epr_get_field_elem_as_str: invalid type");
        return NULL;
    }
    return (const char*) field->elems;
}




/**
 * Copies the data of the given field into the given buffer of <code>double</code>
 * elements. The actual number of elements copied is the minimum of the given
 * number of elements (the buffer's size) and the actual number of elements contained
 * in the field.
 * <p>If the actual field data type is not <code>double</code> the function automatically
 * performs the conversion.
 *
 * @param field the field from which to copy the elements
 * @param buffer the buffer in which to copy the data
 * @param num_elems the number of elements in the given buffer
 * @return the actual number of elements copied
 */
uint epr_copy_field_elems_as_doubles(const EPR_SField* field, double* buffer, uint num_elems)
{
    uint num_elems_min = 0;
    uint i;

    epr_clear_err();

    if (field == NULL) {
        epr_set_err(e_err_invalid_field_name,
                    "epr_get_field_elems_as_doubles: invalid field name");
        return num_elems_min;
    }

    num_elems_min = num_elems;
    if (field->info->num_elems < num_elems_min) {
        num_elems_min = field->info->num_elems;
    }

    if (field->info->data_type_id == e_tid_uchar) {
        for (i = 0; i < num_elems_min; i++) {
            buffer[i] = (double)((uchar*) field->elems)[i];
        }
    } else if (field->info->data_type_id == e_tid_char) {
        for (i = 0; i < num_elems_min; i++) {
            buffer[i] = (double)((char*) field->elems)[i];
        }
    } else if (field->info->data_type_id == e_tid_ushort) {
        for (i = 0; i < num_elems_min; i++) {
            buffer[i] = (double)((ushort*) field->elems)[i];
        }
    } else if (field->info->data_type_id == e_tid_short) {
        for (i = 0; i < num_elems_min; i++) {
            buffer[i] = (double)((short*) field->elems)[i];
        }
    } else if (field->info->data_type_id == e_tid_uint) {
        for (i = 0; i < num_elems_min; i++) {
            buffer[i] = (double)((uint*) field->elems)[i];
        }
    } else if (field->info->data_type_id == e_tid_int) {
        for (i = 0; i < num_elems_min; i++) {
            buffer[i] = (double)((int*) field->elems)[i];
        }
    } else if (field->info->data_type_id == e_tid_float) {
        for (i = 0; i < num_elems_min; i++) {
            buffer[i] = (double)((float*) field->elems)[i];
        }
    } else if (field->info->data_type_id == e_tid_double) {
        for (i = 0; i < num_elems_min; i++) {
            buffer[i] = (double)((double*) field->elems)[i];
        }
    } else {
        epr_set_err(e_err_invalid_data_format,
                    "epr_get_field_elems_as_double: invalid type");
        return 0;
    }

    return num_elems_min;
}

/**
 * Copies the data of the given field into the given buffer of <code>float</code>
 * elements. The actual number of elements copied is the minimum of the given
 * number of elements (the buffer's size) and the actual number of elements contained
 * in the field.
 * <p>If the actual field data type is not <code>float</code> the function automatically
 * performs the conversion.
 *
 * @param field the field from which to copy the elements
 * @param buffer the buffer in which to copy the data
 * @param num_elems the number of elements in the given buffer
 * @return the actual number of elements copied
 */
uint epr_copy_field_elems_as_floats(const EPR_SField* field, float* buffer, uint num_elems)
{
    uint num_elems_min = 0;
    uint i;

    epr_clear_err();

    if (field == NULL) {
        epr_set_err(e_err_invalid_field_name,
                    "epr_copy_field_elems_as_floats: invalid field name");
        return num_elems_min;
    }

    num_elems_min = num_elems;
    if (field->info->num_elems < num_elems_min) {
        num_elems_min = field->info->num_elems;
    }

    if (field->info->data_type_id == e_tid_uchar) {
        for (i = 0; i < num_elems_min; i++) {
            buffer[i] = (float)((uchar*) field->elems)[i];
        }
    } else if (field->info->data_type_id == e_tid_char) {
        for (i = 0; i < num_elems_min; i++) {
            buffer[i] = (float)((char*) field->elems)[i];
        }
    } else if (field->info->data_type_id == e_tid_ushort) {
        for (i = 0; i < num_elems_min; i++) {
            buffer[i] = (float)((ushort*) field->elems)[i];
        }
    } else if (field->info->data_type_id == e_tid_short) {
        for (i = 0; i < num_elems_min; i++) {
            buffer[i] = (float)((short*) field->elems)[i];
        }
    } else if (field->info->data_type_id == e_tid_uint) {
        for (i = 0; i < num_elems_min; i++) {
            buffer[i] = (float)((uint*) field->elems)[i];
        }
    } else if (field->info->data_type_id == e_tid_int) {
        for (i = 0; i < num_elems_min; i++) {
            buffer[i] = (float)((int*) field->elems)[i];
        }
    } else if (field->info->data_type_id == e_tid_float) {
        for (i = 0; i < num_elems_min; i++) {
            buffer[i] = (float)((float*) field->elems)[i];
        }
    } else {
        epr_set_err(e_err_invalid_data_format,
                    "epr_copy_field_elems_as_floats: invalid type");
        return 0;
    }

    return num_elems_min;
}

/**
 * Copies the data of the given field into the given buffer of <code>int</code>
 * elements. The actual number of elements copied is the minimum of the given
 * number of elements (the buffer's size) and the actual number of elements contained
 * in the field.
 * <p>If the actual field data type is not <code>int</code> the function automatically
 * performs the conversion.
 *
 * @param field the field from which to copy the elements
 * @param buffer the buffer in which to copy the data
 * @param num_elems the number of elements in the given buffer
 * @return the actual number of elements copied
 */
uint epr_copy_field_elems_as_longs(const EPR_SField* field, int* buffer, uint num_elems)
{
    uint num_elems_min = 0;
    uint i;

    epr_clear_err();

    if (field == NULL) {
        epr_set_err(e_err_invalid_field_name,
                    "epr_copy_field_elems_as_longs: invalid field name");
        return num_elems_min;
    }

    num_elems_min = num_elems;
    if (field->info->num_elems < num_elems_min) {
        num_elems_min = field->info->num_elems;
    }

    if (field->info->data_type_id == e_tid_uchar) {
        for (i = 0; i < num_elems_min; i++) {
            buffer[i] = (int)((uchar*) field->elems)[i];
        }
    } else if (field->info->data_type_id == e_tid_char) {
        for (i = 0; i < num_elems_min; i++) {
            buffer[i] = (int)((char*) field->elems)[i];
        }
    } else if (field->info->data_type_id == e_tid_ushort) {
        for (i = 0; i < num_elems_min; i++) {
            buffer[i] = (int)((ushort*) field->elems)[i];
        }
    } else if (field->info->data_type_id == e_tid_short) {
        for (i = 0; i < num_elems_min; i++) {
            buffer[i] = (int)((short*) field->elems)[i];
        }
    } else if (field->info->data_type_id == e_tid_int) {
        for (i = 0; i < num_elems_min; i++) {
            buffer[i] = (int)((int*) field->elems)[i];
        }
    } else {
        epr_set_err(e_err_invalid_data_format,
                    "epr_copy_field_elems_as_longs: invalid type");
        return 0;
    }

    return num_elems_min;
}

/**
 * Copies the data of the given field into the given buffer of <code>uint</code>
 * elements. The actual number of elements copied is the minimum of the given
 * number of elements (the buffer's size) and the actual number of elements contained
 * in the field.
 * <p>If the actual field data type is not <code>uint</code> the function automatically
 * performs the conversion.
 *
 * @param field the field from which to copy the elements
 * @param buffer the buffer in which to copy the data
 * @param num_elems the number of elements in the given buffer
 * @return the actual number of elements copied
 */
uint epr_copy_field_elems_as_uints(const EPR_SField* field, uint* buffer, uint num_elems)
{
    uint num_elems_min = 0;
    uint i;

    epr_clear_err();

    if (field == NULL) {
        epr_set_err(e_err_invalid_field_name,
                    "epr_copy_field_elems_as_uints: invalid field name");
        return num_elems_min;
    }

    num_elems_min = num_elems;
    if (field->info->num_elems < num_elems_min) {
        num_elems_min = field->info->num_elems;
    }

    if (field->info->data_type_id == e_tid_uchar) {
        for (i = 0; i < num_elems_min; i++) {
            buffer[i] = (uint)((uchar*) field->elems)[i];
        }
    } else if (field->info->data_type_id == e_tid_char) {
        for (i = 0; i < num_elems_min; i++) {
            buffer[i] = (uint)((char*) field->elems)[i];
        }
    } else if (field->info->data_type_id == e_tid_ushort) {
        for (i = 0; i < num_elems_min; i++) {
            buffer[i] = (uint)((ushort*) field->elems)[i];
        }
    } else if (field->info->data_type_id == e_tid_short) {
        for (i = 0; i < num_elems_min; i++) {
            buffer[i] = (uint)((short*) field->elems)[i];
        }
    } else if (field->info->data_type_id == e_tid_int) {
        for (i = 0; i < num_elems_min; i++) {
            buffer[i] = (uint)((uint*) field->elems)[i];
        }
    } else {
        epr_set_err(e_err_invalid_data_format,
                    "epr_copy_field_elems_as_uints: invalid type");
        return 0;
    }

    return num_elems_min;
}
/***************************************************************************************/
