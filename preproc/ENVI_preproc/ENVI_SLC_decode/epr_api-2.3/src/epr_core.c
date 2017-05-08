/*
 * $Id: epr_core.c,v 1.1.1.1 2004-10-28 19:22:22 norman Exp $
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
 * The one and only ENVISAT API instance.
 */
EPR_SAPI epr_api;

/*
   Function:    epr_str_to_data_type_id
   Access:      private API implementation helper
   Changelog:   2002/01/05  nf initial version
 */
/**
 * Converts the given string into a data type identifier.
 *
 * @param str the string to be converted.
 * @return the data type identifier represented by the given string.
 *         If the string can not be converted the special value
 *         <code>e_tid_unknown</code> is returned.
 */
EPR_EDataTypeId epr_str_to_data_type_id(const char* str)
{
    assert(str != NULL);
    if (epr_equal_names(str, "UChar") || epr_equal_names(str, "uchar"))
        return e_tid_uchar;
    else if (epr_equal_names(str, "AChar") || epr_equal_names(str, "SChar") || epr_equal_names(str, "char"))
        return e_tid_char;
    else if (epr_equal_names(str, "UShort") || epr_equal_names(str, "ushort"))
        return e_tid_ushort;
    else if (epr_equal_names(str, "SShort") || epr_equal_names(str, "short"))
        return e_tid_short;
    else if (epr_equal_names(str, "UInt") || epr_equal_names(str, "uint"))
        return e_tid_uint;
    else if (epr_equal_names(str, "SInt") || epr_equal_names(str, "int"))
        return e_tid_int;
    else if (epr_equal_names(str, "ULong") || epr_equal_names(str, "ulong"))
        return e_tid_uint;
    else if (epr_equal_names(str, "SLong") || epr_equal_names(str, "long"))
        return e_tid_int;
    else if (epr_equal_names(str, "Float") || epr_equal_names(str, "float"))
        return e_tid_float;
    else if (epr_equal_names(str, "Double") || epr_equal_names(str, "double"))
        return e_tid_double;
    else if (epr_equal_names(str, "@/types/UTC.dd") || epr_equal_names(str, "time"))
        return e_tid_time;
    else if (epr_equal_names(str, "String") || epr_equal_names(str, "string"))
        return e_tid_string;
    else if (epr_equal_names(str, "Spare") || epr_equal_names(str, "spare"))
        return e_tid_spare;
    else
        return e_tid_unknown;
}

/*
   Function:    epr_data_type_id_to_str
   Access:      public API implementation helper
   Changelog:   2003/07/10  nf initial version
 */
/**
 * Converts the given data type identifier to a string representing the C-type of the data type.
 *
 * @param data_type_id  the data type identifier.
 *         If the identifier can not be converted, an empty string
 *         <code>""</code> is returned.
 * @return the C-type  string
 */
const char* epr_data_type_id_to_str(EPR_EDataTypeId data_type_id)
{
    switch (data_type_id)
    {
        case e_tid_uchar:
            return "uchar";
        case e_tid_char:
            return "char";
        case e_tid_ushort:
            return "ushort";
        case e_tid_short:
            return "short";
        case e_tid_uint:
            return "uint";
        case e_tid_int:
            return "int";
        case e_tid_float:
            return "float";
        case e_tid_double:
            return "double";
        case e_tid_string:
            return "string";
        case e_tid_spare:
            return "spare";
        case e_tid_time:
            return "time";
        default:
            return "";
    }
}

/*
   Function:    epr_get_data_type_size
   Access:      private API implementation helper
   Changelog:   2002/01/24  nf initial version
 */
/**
 * Determines the length of the given data type identifier.
 *
 * @param data_type_id the data type identifier.
 * @return the the length of the data type identifier.
 *         If the data type identifier is unknown,
 *         <code>e_tid_unknown</code> is returned.
 */
uint epr_get_data_type_size(EPR_EDataTypeId data_type_id)
{
    switch (data_type_id)
    {
        case e_tid_uchar:
            return sizeof(uchar);
        case e_tid_char:
            return sizeof(char);
        case e_tid_ushort:
            return sizeof(ushort);
        case e_tid_short:
            return sizeof(short);
        case e_tid_uint:
            return sizeof(uint);
        case e_tid_int:
            return sizeof(int);
        case e_tid_float:
            return sizeof(float);
        case e_tid_double:
            return sizeof(double);
        case e_tid_string:
            return sizeof(char);
        case e_tid_spare:
            return sizeof(uchar);
        case e_tid_time:
            return sizeof (int)   /* days */
                +  sizeof (uint)  /* seconds */
                +  sizeof (uint); /* microseconds */
        default:
            return 0;
    }
}
/*
 * ===================== Logging ==============================
 */

/*
   Function:    epr_log
   Access:      private API implementation helper
   Changelog:   2002/01/05  mp initial version
 */
/**
 * Logs a message with the given log level.
 *
 * <p> The function calls this API's log handler if
 * it is not <code>NULL</code> and if the given log
 * level is greater than or equal to the global log level.
 *
 * @param log_level the log level (or message type) for the mesage
 * @param log_message the mesage
 */
void epr_log(EPR_ELogLevel log_level, const char* log_message)
{
    if (epr_api.log_handler != NULL && log_level >= epr_api.log_level)
        epr_api.log_handler(log_level, log_message);
}

/*
 * ===================== Error Handling ==============================
 */

/*
   Function:    epr_set_error
   Access:      private API implementation helper
   Changelog:   2002/01/05  mp initial version
 */
/**
 * Sets the given error code and the associated error message and
 * calls this API's error handler if it is not <code>NULL</code>.
 *
 * @param err_code the error code
 * @param err_message the error mesage
 */
void epr_set_err(EPR_EErrCode err_code, const char* err_message)
{
    epr_api.last_err_code = err_code;
    epr_assign_string(&epr_api.last_err_message, err_message);

    if (epr_api.log_handler != NULL)
    {
        epr_api.log_handler(e_log_error, err_message);
    }

    if (epr_api.err_handler != NULL)
    {
        epr_api.err_handler(err_code, err_message);
    }
}

/*
   Function:    epr_set_error
   Access:      private API implementation helper
   Changelog:   2002/01/05  mp initial version
 */
/**
 * Clears the last error. After calling this function, calling
 * <code>epr_get_last_err_code</code> returns <code>e_err_none</code> or zero and
 * <code>epr_get_last_err_message</code> returns <code>NULL</code>.
 */
void epr_clear_err()
{
    epr_api.last_err_code = e_err_none;
    epr_free_string(epr_api.last_err_message);
    epr_api.last_err_message = NULL;
}

/*
   Function:    epr_get_last_err_code
   Access:      private API implementation helper
   Changelog:   2002/01/05  nf initial version
 */
/**
 * Sets the error code of the error that occured during
 * the last API function call.
 *
 * @return the error code, <code>e_err_none</code> or zero if no error occured
 */
EPR_EErrCode epr_get_last_err_code()
{
    return epr_api.last_err_code;
}

/*
   Function:    epr_get_last_err_message
   Access:      private API implementation helper
   Changelog:   2002/01/05  nf initial version
 */
/**
 * Sets the error message of the error that occured during
 * the last API function call.
 *
 * @return the error message, <code>NULL</code> if no error occured
 */
const char* epr_get_last_err_message()
{
    return epr_api.last_err_message;
}


/**
 * Opens a file to read.
 *
 * @param file_path the path to the file.
 *
 * @return the file handle or
 *         <code>NULL</code> if an error occured.
 */
FILE* epr_open_file(char* file_path)
{
    FILE* db_file;
	epr_log(e_log_debug, "about to open file: ");
	epr_log(e_log_debug, file_path);

    db_file = fopen(epr_trim_string(file_path), "rb");
    if (db_file == NULL)
    {
        epr_log(e_log_debug, "open failed");
        if (errno == ENOENT)
        {
            epr_set_err(e_err_file_not_found,
                    "epr_open_file: file not found");
        }
        else {
            epr_set_err(e_err_file_access_denied,
                    "epr_open_file: file open failed");
        }
    }
    /*epr_log(e_log_debug, "open successful");*/
    return db_file;
}


/**
 * Converts the given string into an int number.
 *
 * @param str the string to be converted.
 * @return the int type number represented by the given string.
 *         If the string can not be converted,
 *         the value of <code>1</code> is returned.
 */
int epr_str_to_number(const char* str)
{
   char   *stopstring;
   int   l;
   assert(str != NULL);

   if (strcmp(str, "*") == 0) return 1;
   if (strcmp(str, "") == 0) return 1;

    errno = 0;
    l = strtol( str, &stopstring, 10 );

    if (errno == EDOM)
    {
        epr_set_err(e_err_illegal_conversion,
                "failed to convert string to integer: errno = EDOM");
        return -1;
    }
    if (errno == ERANGE)
    {
        epr_set_err(e_err_illegal_conversion,
                "failed to convert string to integer: errno = ERANGE");
        return -1;
    }

    return l;
}


/**
 * Converts the given string into a field length.
 *
 * The string can represent a single integer value or a sequence
 * of integer value and parameter references (names). Integers
 * and value are expected to be separated by the asterisk
 * character ('*'). E.g. the string "3 * 4 * num_pixels_across"
 * is represents a valid field length as long as the parameter
 * name 'num_pixels_across' is found in the given parameter table.
 *
 * @param count the string to be converted
 * @param product_id the Product identifier containing the values
 *                    for the Product
 * @return the field length computed from the given string or
 *         <code>(uint)-1</code> if an error occured.
 */
uint epr_parse_value_count(EPR_SProductId* product_id, const char* count)
{
    char seps[] = EPR_COUNT_SEPARATOR_ARRAY;
    char * token;
    char * str;
    uint c;
    int pos = 0;
    uint comes_from_convert;

    c = 1;

    str = epr_clone_string(count);
    if (str == NULL)
    {
        epr_set_err(e_err_out_of_memory,
                   "epr_parse_value_count: cannot allocate memory");
        return 16111;/*8999*/  /* @todo check: why in the hell 16111 ? */
    }

    while ((token = epr_str_tok(str, seps, &pos)) != NULL)
    {
        comes_from_convert = epr_str_to_number(token);
        if (comes_from_convert == 0 && (strcmp(token, "0") != 0))
        {
            /* check if "token" belongs to param_table*/
            comes_from_convert = epr_param_to_value(token, product_id->param_table);
            if (comes_from_convert == -1)
            {
                epr_set_err(e_err_illegal_conversion,
                   "epr_parse_value_count: parameter not found");
                return 16111; /* @todo check: why in the hell 16111 ? */
            }
        }
        c *= comes_from_convert;
        epr_free_string(token);
        token = NULL;
    }
    epr_free_string(str);
    return c;
}


/**
 * Finds in the param_table the value corresponding to its name.
 *
 * @param str the parameter name
 * @param param_table the pointer to param_table
 *
 * @return the value of the given name or
 *         <code>(uint)-1</code> if an error occured.
 */
uint epr_param_to_value(const char* str, EPR_SPtrArray* param_table)
{
    EPR_SParamElem* param_elem = NULL;
    int elem_index;
    int param_table_length;

    param_table_length = (int)epr_get_ptr_array_length(param_table);
    for (elem_index = 0; elem_index < param_table_length; elem_index++)
    {
        param_elem = (EPR_SParamElem*)epr_get_ptr_array_elem_at(param_table, elem_index);
        if (epr_equal_names(param_elem->param_name, str))
        {
            return param_elem->param_value;
        }
    }
    return (uint) -1; /*error*/
}


/**
 * Adapts path description to operating system.
 *
 * @param path the path to a file.
 *
 */
void epr_make_os_compatible_path(char* path)
{
    if (path != NULL)
    {
        char* pc = path;
        while (*pc != '\0')
        {
#ifdef WIN32
            if (*pc == '/')
                *pc = '\\';

#elif _M_MPPC
            if (*pc == '/')
                *pc = ':';	/* UK note: the colon is an old-style path separator of the Mac OS */
							/* possibly not used any more but supported for Classic compatibility */
							/* @to do: check whether the forward slash / should be used in Mac OS X */
#else
            if (*pc == '\\')
                *pc = '/';
#endif
            pc++;
        }
    }
}


epr_boolean epr_check_api_init_flag() {
    if (!epr_api.init_flag) {
        epr_set_err(e_err_api_not_initialized,
                    "epr_open_product: API not initialized (forgot to call epr_init_api?)");
        return FALSE;
    }
    return TRUE;
}
