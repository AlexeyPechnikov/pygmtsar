/*
 * $Id: epr_msph.c,v 1.1.1.1 2004-10-28 19:22:23 norman Exp $
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
 * Reads the full main product header (MPH) of the ENVISAT product file
 * by the given product identifier.
 *
 * @param product_id the product identifier, must not be <code>NULL</code>
 * @return a record representing the MPH of the specified product file
 *         or <code>NULL</code> if an error occured.
 */

EPR_SRecord* epr_read_mph(EPR_SProductId* product_id)
{
    EPR_SRecord* record = NULL;
    char* code_block;
    int numread;

    epr_clear_err();

    code_block = epr_create_string(EPR_MPH_SIZE);
    if (code_block == NULL)
    {
        epr_set_err(e_err_out_of_memory,
                    "epr_read_mph: out of memory");
        return NULL;
    }
    rewind(product_id->istream);
    numread = fread(code_block, 1, EPR_MPH_SIZE, product_id->istream);

    if (numread != EPR_MPH_SIZE)
    {
        epr_set_err(e_err_file_read_error,
            "epr_read_mph: wrong reading MPH from product data file");
        return NULL;
    }
    record = epr_parse_header("mph", code_block);
    if (record == NULL)
    {
        epr_set_err(e_err_invalid_record,
            "epr_read_mph: can not recognize the correct MPH from product data file");
    } else {
        epr_add_ptr_array_elem(product_id->record_info_cache, record->info);
    }

    epr_free_string(code_block);
    return record;
}

/**
 * Reads the full specific product header (SPH) of the ENVISAT product file
 * by the given product identifier.
 *
 * @param product_id the product identifier, must not be <code>NULL</code>
 * @return a record representing the MPH of the specified product file
 *         or <code>NULL</code> if an error occured.
 */
EPR_SRecord* epr_read_sph(EPR_SProductId* product_id)
{
    EPR_SRecord* sph_record = NULL;
    const EPR_SField* field;

    char* code_block;
    int numread;
    uint sph_length = 0;
    uint sph_without_dsd_length = 0;
    uint dsd_number = 0;

    epr_clear_err();

    if (product_id->mph_record == NULL) {
        product_id->mph_record = epr_read_mph(product_id);
        if (product_id->mph_record == NULL) {
            epr_set_err(e_err_file_read_error, "epr_read_sph: wrong MPH");
            return NULL;
        }
    }

    field = epr_get_field(product_id->mph_record, "SPH_SIZE");
    sph_length = ((uint*) field->elems)[0];
    if (sph_length == 0) {
        epr_set_err(e_err_invalid_value,
            "epr_read_sph: wrong MPH: SPH_SIZE must be > 0");
        return NULL;
    }
    field = epr_get_field(product_id->mph_record, "NUM_DSD");
    dsd_number = ((uint*) field->elems)[0];
    if (dsd_number == 0) {
        epr_set_err(e_err_invalid_value,
            "epr_read_sph: wrong MPH: NUM_DSD must be > 0");
        return NULL;
    }

    epr_api.epr_head_size = sph_length + EPR_MPH_SIZE;
    if (fseek(product_id->istream, EPR_MPH_SIZE, SEEK_SET) != 0) {
        epr_set_err(e_err_file_access_denied,
                    "epr_read_sph: file seek failed");
        return NULL;
    }

    sph_without_dsd_length = sph_length - dsd_number * EPR_DSD_SIZE;

    code_block = epr_create_string(sph_without_dsd_length);
    numread = fread(code_block, 1, sph_without_dsd_length, product_id->istream);
    if ((uint)numread != sph_without_dsd_length) {
        epr_set_err(e_err_file_read_error,
            "epr_read_sph: wrong reading SPH from product data file");
        return NULL;
    }

    sph_record = epr_parse_header("sph", code_block);
    if (sph_record == NULL) {
        epr_set_err(e_err_invalid_record,
            "epr_read_sph: can not recognize the correct SPH from product data file");
    } else {
        epr_add_ptr_array_elem(product_id->record_info_cache, sph_record->info);
    }

    epr_free_string(code_block);
    return sph_record;
}

void epr_store_header(const char* header_name, const char* ascii_source) {
    FILE* os;
    char fname[1024];
    sprintf(fname, "%s.txt", header_name);
    os=fopen(fname, "w");
    fprintf(os,"%s", ascii_source);
    fclose(os);
}


/**
 * Parses the header ASCII information.
 *
 * @param header_name name of the header ascii information;
 * @param ascii_source the header ascii information was read;
 * @param record the identifier of header ascii information.
 */
EPR_SRecord* epr_parse_header(const char* header_name, const char* ascii_source)
{
    EPR_SRecordInfo* record_info;
    EPR_SPtrArray* field_infos = NULL;
    EPR_SFieldInfo* field_info;
    EPR_SPtrArray* header_values = NULL;
    EPR_SRecord* record = NULL;
    EPR_EDataTypeId tp;
    char * code_block;
    char seps[] = EPR_HEADER_SEPARATOR_ARRAY;
    char * token_name;
    char * token_value;
    char * token_unit;
    char * h_name;
    int pos = 0;
    int pos_ascii = 0;
    uint num_bytes = 0;
    uint num_elems = 0;

    epr_clear_err();

    /* uncomment for debugging purpose */
    /* epr_store_header(header_name, ascii_source); */

    header_values = epr_create_ptr_array(16);
    field_infos = epr_create_ptr_array(16);
    h_name = epr_clone_string(header_name);

    while ((code_block = epr_str_tok(ascii_source, "\n", &pos_ascii)) != NULL) {
        /*if EMPTY code_block*/
        if ((strlen(code_block) > 0) && (code_block[0] == ' ')) {
            /* epr_log(e_log_info, "code_block is empty"); */
            if (code_block != NULL) {
                epr_free_string(code_block);
                code_block = NULL;
            }
            continue;
        }
        /*if '=' separator*/
        pos = 0;
        token_name = epr_str_tok(code_block, seps, &pos);
        if (pos == 1) {
            epr_set_err(e_err_invalid_keyword_name,
                        "epr_parse_header: invalid ascii header: keyword is empty");
            epr_free_string(token_name);
            if (code_block != NULL) {
                epr_free_string(code_block);
                code_block = NULL;
            }
            continue;
        }
        if (pos == (int)strlen(code_block) + 1) {
            epr_set_err(e_err_invalid_keyword_name,
                        "epr_parse_header: invalid ascii header: keyword not found");
            epr_free_string(token_name);
            if (code_block != NULL) {
                epr_free_string(code_block);
                code_block = NULL;
            }
            continue;
        }
        /*if STRING value*/
        if (code_block[pos] == '\"') {
            pos ++;
			/*
			  Note that strings always are considered as one single element,
			  so we get the total number of characters from tot_size.
			  Addidionally we reserve an extra character for the trailing zero (terminator),
			*/
            token_value = epr_strip_string_r(epr_str_tok(code_block, "\"", &pos));
            token_unit = NULL;
            tp = e_tid_string;
            num_bytes = (uint)strlen(token_value);
            num_elems = 1;
            epr_add_ptr_array_elem(header_values, token_value);
        } else {
            token_value = epr_str_tok(code_block, seps, &pos);
            if (token_value == NULL) {
                epr_set_err(e_err_invalid_value,
                            "epr_parse_header: invalid ascii header: value not found");
                token_value = epr_clone_string("");
                token_unit = NULL;
                tp = e_tid_uchar;
                num_bytes = 0;
                num_elems = 1;
                epr_add_ptr_array_elem(header_values, token_value);
            } else {
                /*if FLOAT-DOUBLE value*/
                if (strchr(token_value, '.') != NULL
                    || strchr(token_value, 'e') != NULL
                    || strchr(token_value, 'E') != NULL)
                {
                    epr_parse_double_token(header_values, token_value, &num_elems, &num_bytes, &tp);
                    token_unit = epr_str_tok(code_block, seps, &pos);

                    epr_free_string(token_value);
                    token_value = NULL;

                /*if INTEGER_LONG value*/
                } else if ((strlen(token_value) > 1)) {

                    epr_parse_int_token(header_values, token_value, &num_elems, &num_bytes, &tp);

                    epr_free_string(token_value);
                    token_value = NULL;

                    token_unit = epr_str_tok(code_block, seps, &pos);
                } else {
                    /*if CHAR value*/
                    if (strlen(token_value) > 1) {
                        epr_set_err(e_err_invalid_value,
                                    "epr_parse_header: invalid ascii header: illegal value");
                        token_value = epr_clone_string("");
                        token_unit = NULL;
                        tp = e_tid_uchar;
                        num_bytes = 0;
                        num_elems = 1;
                        epr_add_ptr_array_elem(header_values, token_value);
                        epr_free_string(token_name);
                        if (code_block != NULL) {
                            epr_free_string(code_block);
                            code_block = NULL;
                        }
                        continue;
                    } else {
                        token_unit = NULL;
                        tp = e_tid_uchar;
                        num_bytes = (uint)strlen(token_value);
                        num_elems = 1;
                        epr_add_ptr_array_elem(header_values, token_value);
                    }
                }
            }
        }
        field_info = epr_create_field_info(tp, h_name, token_name, num_elems, num_bytes, 1, token_unit);
        epr_add_ptr_array_elem(field_infos, field_info);
        epr_free_string(token_name);
        epr_free_string(token_unit);
        epr_free_string(code_block);
    }

    if (field_infos->length > 0) {
        record_info = epr_create_record_info(h_name, field_infos);
        record = epr_create_record_from_info(record_info);
        epr_set_header_field_values(record, header_values);
    }

    epr_free_char_ptr_array(header_values);

    epr_free_string(h_name);

    return record;
}


void epr_parse_string_token(EPR_SPtrArray* header_values, char* token_value, uint* num_elems, uint* num_bytes, EPR_EDataTypeId* tp)
{
    char exceptions[] = EPR_HEADER_EXCEPTIONS_ARRAY;
    char * token_value_o;
    char * tmp;
    uint pos_value = 0;
    int cyc = 0;

    pos_value = 0;
    *num_elems = 0;
    while ((tmp = epr_str_tok_tok(token_value + 1, "+-", exceptions, &pos_value)) != NULL) {
        cyc ++;
        token_value_o = epr_create_string(strlen(tmp) + 1);
        if (strlen(tmp) == strlen(token_value) - 1) {
            token_value_o[0] = token_value[0];
        } else if (pos_value < (uint)strlen(token_value) - 1) {
            token_value_o[0] = token_value[pos_value - strlen(tmp) - 1];
        } else {
            token_value_o[0] = token_value[pos_value - strlen(tmp)];
        }
        strcat(token_value_o, tmp);
        epr_add_ptr_array_elem(header_values, token_value_o);
        epr_free_string(tmp);
    }
    *num_bytes = sizeof(double);
    *tp = e_tid_double;
    *num_elems = cyc;
}

void epr_parse_double_token(EPR_SPtrArray* header_values, char* token_value, uint* num_elems, uint* num_bytes, EPR_EDataTypeId* tp)
{
    char exceptions[] = EPR_HEADER_EXCEPTIONS_ARRAY;
    char * token_value_o;
    char * tmp;
    uint pos_value = 0;
    int cyc = 0;

    pos_value = 0;
    *num_elems = 0;
    while ((tmp = epr_str_tok_tok(token_value + 1, "+-", exceptions, &pos_value)) != NULL) {
        cyc ++;
        token_value_o = epr_create_string(strlen(tmp) + 1);
        if (strlen(tmp) == strlen(token_value) - 1) {
            token_value_o[0] = token_value[0];
        } else if (pos_value < (uint)strlen(token_value) - 1) {
            token_value_o[0] = token_value[pos_value - strlen(tmp) - 1];
        } else {
            token_value_o[0] = token_value[pos_value - strlen(tmp)];
        }
        strcat(token_value_o, tmp);
        epr_add_ptr_array_elem(header_values, token_value_o);
        epr_free_string(tmp);
    }
    *num_bytes = sizeof(double);
    *tp = e_tid_double;
    *num_elems = cyc;
}


void epr_parse_int_token(EPR_SPtrArray* header_values, char* token_value, uint* num_elems, uint* num_bytes, EPR_EDataTypeId* tp)
{
    char * token_value_o;
    char * tmp;
    char * tmp_v;
    char * stopstring;
    int pos_value = 0;
    uint dlina;
    uint i;
    char value_buffer[32];
    int lmp;
    uint ulmp;
    int flag_int = 0;
    int flag_negative = 0;
    int cyc = 0;

    pos_value = 0;
    *num_elems = 0;
    flag_int = 0;
    flag_negative = 0;

    if (strchr(token_value, '-') != NULL) {
        flag_int = 1;
        *num_bytes = sizeof(int);
        *tp = e_tid_int;
    } else {
        *num_bytes = sizeof(uint);
        *tp = e_tid_uint;
    }

    while ((tmp = epr_str_tok(token_value + 1, "+-", &pos_value)) != NULL) {
        if (epr_if_no_letters(tmp) == 0) {
            epr_set_err(e_err_invalid_value,
                        "epr_parse_header: invalid ascii header: illegal value");
            cyc ++;
            tmp = epr_clone_string("-999999");
            *num_bytes = sizeof(int);
            *tp = e_tid_int;
            epr_add_ptr_array_elem(header_values, tmp);
        } else {
            cyc ++;
            token_value_o = epr_create_string(strlen(tmp) + 1);
            if (strlen(tmp) == strlen(token_value) - 1) {
                token_value_o[0] = token_value[0];
            } else if (pos_value < (int)strlen(token_value) - 1) {
                token_value_o[0] = token_value[pos_value - strlen(tmp) - 1];
            } else if (strlen(tmp) == 1) {
                if (cyc == 1)
                    token_value_o[0] = token_value[pos_value];
                else
                    token_value_o[0] = token_value[pos_value - 1];
            } else {
                token_value_o[0] = token_value[pos_value - strlen(tmp)];
            }
            strcat(token_value_o, tmp);
            dlina = (uint)strlen(token_value_o);
            tmp_v = epr_create_string(dlina);
            /*if int*/
            if (flag_int == 1) {
                lmp = strtol(token_value_o, &stopstring, 10);
                if (lmp != 0) {
                    tmp_v[0] = token_value_o[0];
                    for (i = 1; i < dlina; i ++) if (token_value_o[i] != '0') break;
                    if (token_value_o[0] == '+') strncpy(tmp_v + 0, token_value_o + i, dlina - i);
                    if (token_value_o[0] == '-') strncpy(tmp_v + 1, token_value_o + i, dlina - i);
                    sprintf(value_buffer, "%d", lmp);
                    /*if int value too large*/
                    if (strcmp(tmp_v, value_buffer) != 0)
                        epr_log(e_log_warning, "product header: int integer value out of range");
                }
            } else if (flag_int == 0) {
                ulmp = strtoul(token_value_o, &stopstring, 10);
                if (ulmp != 0UL) {
                    tmp_v[0] = token_value_o[0];
                    for (i = 1; i < dlina; i ++) if (token_value_o[i] != '0') break;
                    strncpy(tmp_v, token_value_o + i, dlina - i);
                    sprintf(value_buffer, "%u", ulmp);
                    /*if uint value too large*/
                    if (strcmp(tmp_v, value_buffer) != 0)
                        epr_log(e_log_warning, "product header: unsigned int integer value out of range");
                }
            }
            epr_free_string(tmp_v);
            epr_add_ptr_array_elem(header_values, token_value_o);
            epr_free_string(tmp);
        }
    }
    *num_elems = cyc;
}



/**
 * Fills the record for the header ASCII information.
 *
 * @param record to fill;
 * @param header_values values from the given product file;
 * @param record beeng filling.
 */
void epr_set_header_field_values(EPR_SRecord* record, EPR_SPtrArray* header_values)
{
    EPR_SFieldInfo* field_info;
    EPR_SField* field;
    uint ptr_index = 0;
    uint field_index;
    uint field_info_index;
    char * tmp;
    char * stopstring;

    assert(header_values != NULL);

    for (field_index = 0; field_index < record->num_fields; field_index++) {
        field = record->fields[field_index];
        field_info = (EPR_SFieldInfo*)epr_get_ptr_array_elem_at(record->info->field_infos, field_index);

        for (field_info_index = 0; field_info_index < field->info->num_elems; field_info_index++) {
            tmp = (char*)epr_get_ptr_array_elem_at(header_values, ptr_index);
            switch (field_info->data_type_id) {
                case e_tid_uchar:
                    *(((uchar*)field->elems) + field_info_index) = (uchar) tmp[field_info_index];
                    break;
                case e_tid_int:
                    *(((int*)field->elems) + field_info_index) = strtol(tmp, &stopstring, 10);
                    break;
                case e_tid_uint:
                    *(((uint*)field->elems) + field_info_index) = strtoul(tmp, &stopstring, 10);
                    break;
                case e_tid_string:;
                    /*epr_assign_string(&(char*)field->elems, tmp);*/
                    strncpy((char*)field->elems, tmp, field->info->tot_size);
                    break;
                case e_tid_double:
                    *(((double*)field->elems) + field_info_index) = strtod(tmp, &stopstring);
                    break;
                default:
                    epr_set_err(e_err_invalid_value,
                            "epr_set_header_field_values: internal error: illegal value type");
            }
            ptr_index ++;
        }
    }
}


uint epr_compare_param(EPR_SProductId* product_id)
{
    EPR_SDSD* dsd;
    uint of;

    epr_clear_err();

    if (product_id == NULL) {
        epr_set_err(e_err_invalid_product_id, "epr_compare_param: invalid product identifier");
        return 0UL;
    }

    of = epr_api.epr_head_size;
    dsd = (EPR_SDSD*)epr_get_ptr_array_elem_at(product_id->dsd_array, 0);
    if (dsd->ds_offset == epr_api.epr_head_size)
        return of;

    return 0UL;
}
