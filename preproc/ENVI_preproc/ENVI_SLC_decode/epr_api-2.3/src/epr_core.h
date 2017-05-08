/*
 * $Id: epr_core.h,v 1.1.1.1 2004-10-28 19:22:22 norman Exp $
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

#ifndef EPR_CORE_H_INCL
#define EPR_CORE_H_INCL

#ifdef __cplusplus
extern "C"
{
#endif


#include <stdio.h> /* just to get the ANSI-C type FILE */

#include "epr_ptrarray.h"

struct EPR_API;
struct EPR_Parameter;

typedef struct EPR_API EPR_SAPI;
typedef struct EPR_Parameter EPR_SParameter;


#define EPR_ENVISAT_PRODUCT_MERIS        "MER"
#define EPR_ENVISAT_PRODUCT_ASAR         "ASA"
#define EPR_ENVISAT_PRODUCT_SAR          "SAR"
#define EPR_ENVISAT_PRODUCT_AATSR        "ATS"

#define EPR_LONGI_BAND_NAME               "longitude"

#define EPR_AATSR_LINES_PER_TIE_PT        32

#define EPR_MPH_SIZE                      1247
#define EPR_SPH_SIZE                      11622
#define EPR_DSD_SIZE                      280

#define EPR_PRODUCT_MPH_SIZE              1048 /*??????*/
#define EPR_PRODUCT_MAGIC_STR             "PRODUCT=\""
#define EPR_PRODUCT_ID_OFFSET             9

#define EPR_PRODUCT_TYPE_ID_STRLEN        10
#define EPR_LENGTH_NUM_DSD_IDENTIFIER     9
#define EPR_COUNT_SEPARATOR_ARRAY         ",*"
#define EPR_IRRELEVANCE_SYMBOL            '*'
#define EPR_FIELD_SEPARATOR_ARRAY         "|"
#define EPR_HEADER_SEPARATOR_ARRAY        "=<>"
#define EPR_HEADER_EXCEPTIONS_ARRAY       "eE"

#define EPR_LONGI_ABS_MAX                 180
#define EPR_LONGI_ABS_MIN                 -180

#define EPR_LINE_MAX_LENGTH               2000

#define EPR_BE_MAGIC_NUMBER  1162761801UL
#define EPR_LE_MAGIC_NUMBER  1230392901UL
#define EPR_LE_MAGIC_BYTE_0  'E'
#define EPR_LE_MAGIC_BYTE_1  'N'
#define EPR_LE_MAGIC_BYTE_2  'V'
#define EPR_LE_MAGIC_BYTE_3  'I'

#define EPR_ATS_NUM_PER_POINT_ACROSS_LOCAT      23
#define EPR_ATS_NUM_PER_POINT_ACROSS_SOLAR      11
#define EPR_ATS_LINE_LENGTH                     512

#define EPR_ASAR_NUM_PER_POINT_ACROSS_LOCAT     11

/**
 * The <code>EPR_API</code> structure is a container for all globally
 * required information related to to the ENVISAT product reader API.
 *
 * <p> A single global (but hidden) instance exists for this structure.
 */
struct EPR_API
{
    /**
     * The directory path to the record info database.
     */
    epr_boolean init_flag;

    /**
     * A boolean value indicating whether this code run's on a
     * little endian order machine or not.
     * <p><code>1</code> stands for little endian (LE),
     * <code>0</code> stands for big endian (BE).
     */
    int little_endian_order;

    /**
     * A unsigned int value indicating head length
     */
    uint epr_head_size;

    /**
     * The directory path to the record info database.
     */
    /*char* db_dir_path;*/

    /**
     * The current log level for the ENVISAT API.
     */
    EPR_ELogLevel log_level;

    /**
     * The log handler (function pointer) for the ENVISAT API.
     * Can be <code>NULL</code>.
     */
    EPR_FLogHandler log_handler;

    /**
     * The error code of the last error occured.
     */
    EPR_EErrCode last_err_code;

    /**
     * The error message of the last error occured.
     */
    char* last_err_message;

    /**
     * The error handler (function pointer) for the ENVISAT API.
     * Can be <code>NULL</code>.
     */
    EPR_FErrHandler err_handler;
};


/**
 * The one and only ENVISAT API instance.
 */
extern EPR_SAPI epr_api;


/*
 * ======================================================================
 */

/**
 * Frees the memory allocated by the given product file identifier.
 *
 * <p> After calling this function the give product file identifier pointer
 * gets invalid and should not be used anymore.
 *
 * @param record the product file identifier to be released,
 *        if <code>NULL</code> the function immediately returns
 */
void epr_free_product_id(EPR_SProductId* product_id);


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
void epr_log(EPR_ELogLevel log_level, const char* log_message);

/**
 * Sets the given error code and the associated error message and
 * calls this API's error handler if it is not <code>NULL</code>.
 *
 * @param err_code the error code
 * @param err_message the error mesage
 */
void epr_set_err(EPR_EErrCode err_code, const char* err_message);

/**
 * Reads the full main product header (MPH) of the ENVISAT product file
 * given by the given product identifier.
 *
 * @param product_id the product identifier, must not be <code>NULL</code>
 * @return a record representing the MPH of the specified product file
 *         or <code>NULL</code> if an error occured.
 */
EPR_SRecord* epr_read_mph(EPR_SProductId* product_id);

/**
 * Reads the full specific product header (SPH) of the ENVISAT product file
 * given by the given product identifier.
 *
 * @param product_id the product identifier, must not be <code>NULL</code>
 * @return a record representing the MPH of the specified product file
 *         or <code>NULL</code> if an error occured.
 */
EPR_SRecord* epr_read_sph(EPR_SProductId* product_id);

/**
 * Converts the given string into a data type identifier.
 *
 * @param str the string to be converted.
 * @return the data type identifier represented by the given string.
 *         If the string can not be converted the special value
 *         <code>e_tid_unknown</code> is returned.
 */
EPR_EDataTypeId epr_str_to_data_type_id(const char* str);

/*
 * Converts the given string into a field length.
 *
 * The string can represent a single integer value or a sequence
 * of integer value and parameter references (names). Integers
 * and value are expected to be separated by the asterisk
 * character ('*'). E.g. the string "3 * 4 * num_pixels_across"
 * is represents a valid field length as long as the parameter
 * name 'num_pixels_across' is found in the given parameter table.
 *
 * @param str the string to be converted
 * @param param_table the parameter table containing the values
 *                    for the parameter references in the string
 * @return the field length computed from the given string or
 *         <code>(uint)-1</code> if an error occured.
 */
/*uint epr_str_to_field_length(const char* str, EPR_SParamTable* param_table);*/

/**
 * Compares the two given names and returns <code>TRUE</code> if
 * they are equal ignoring the case of each letter.
 *
 * <p> This function is used to compare names throughout the
 * ENVISAT product reader API.
 *
 * @param name1 the first name, must not be NULL
 * @param name2 the second name, must not be NULL
 * @return  <code>TRUE</code> if the names are equal,
 *          <code>FALSE</code> otherwise
 */

char* epr_build_db_file_istream_name(EPR_SProductId* product_id, char* what);
FILE* epr_open_file(char* path_to_file);
int epr_str_to_number(const char* str);
uint epr_parse_value_count(EPR_SProductId* product_id, const char* str);
uint epr_param_to_value(const char* str, EPR_SPtrArray* param_table);
void epr_make_os_compatible_path(char* path);
epr_boolean epr_check_api_init_flag();

/*
void epr_make_image_header(EPR_SProductId* product_id, EPR_SDatasetId* dataset_id, EPR_SRecord* record);
int epr_make_image(EPR_SProductId* product_id, EPR_SDatasetId* dataset_id, EPR_SRecord* record);
*/

/**
 * Gets the element content to output FILE stream.
 *
 * @param field the pointer at the field to get out.
 * @param istream the identifier of the output file.
 */
void epr_output_element(const EPR_SField* field, uint field_index, uint element_index, FILE* istream);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* #ifndef EPR_CORE_H_INCL */
