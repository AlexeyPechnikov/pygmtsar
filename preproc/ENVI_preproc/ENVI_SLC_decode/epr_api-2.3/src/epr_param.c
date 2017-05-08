/*
 * $Id: epr_param.c,v 1.1.1.1 2004-10-28 19:22:23 norman Exp $
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
 * Creates a new element (pair 'name-value') for the parameter table.
 *
 * @param param_name the name of the parameter,
 * @param param_value the value of this parameter,
 * @return the pointer at this element
 *         or <code>NULL</code> if an error occured.
*/
EPR_SParamElem* epr_create_param_elem(const char* param_name, int param_value)
{
    EPR_SParamElem* param_elem = NULL;
/*
    LINE_LENGTH=+02241<samples>
    -1
    /
    LINES_PER_TIE_PT=+064
    +1
    =
    */
    param_elem = (EPR_SParamElem*) calloc(1, sizeof (EPR_SParamElem));
    if (param_elem == NULL) {
        epr_set_err(e_err_out_of_memory,
                    "epr_create_field_info: out of memory");
        return NULL;
    }

    param_elem->param_name = epr_clone_string(param_name);
    param_elem->param_value = param_value;

    return param_elem;
}

EPR_SPtrArray* epr_create_param_table()
{
    EPR_SPtrArray* param_table = NULL;
    param_table = epr_create_ptr_array(16);
    return param_table;
}

/**
 * Frees the memory allocated by the given param_table.
 *
 * <p> After calling this function the give record_info pointer gets
 * invalid and should not be used anymore.
 *
 * @param param_table the table to be released, if <code>NULL</code>
 *        the function immediately returns
 */
void epr_free_param_table(EPR_SPtrArray* param_table)
{
    EPR_SParamElem* param_elem = NULL;
    int param_index = 0;

    if (param_table == NULL)
        return;

    for (param_index = 0; param_index < (int)param_table->length; param_index++) {
        param_elem = (EPR_SParamElem*)epr_get_ptr_array_elem_at(param_table, param_index);
        epr_free_param_elem(param_elem);
    }

    epr_free_ptr_array(param_table);
}


/**
 * Frees the memory allocated by the given param_elem.
 *
 * <p> After calling this function the give record_info pointer gets
 * invalid and should not be used anymore.
 *
 * @param param_table the table to be released, if <code>NULL</code>
 *        the function immediately returns
 */
void epr_free_param_elem(EPR_SParamElem* param_elem)
{

    if (param_elem == NULL)
        return;
    epr_free_string(param_elem->param_name);
    param_elem->param_name = NULL;
    param_elem->param_value = 0;

    free(param_elem);
    return;
}


int epr_set_dyn_dddb_params(EPR_SProductId* product_id)
{
    const EPR_SField* field;
    const EPR_SField* product_field;
    char* tmp;
    EPR_SParamElem* param_elem = NULL;

    uint line_length = 0;
    uint num_tie_points_across = 0;
    uint ntpa = 0;

    product_field = epr_get_field(product_id->mph_record, "PRODUCT");
    tmp = epr_sub_string((char*)product_field->elems, 0, 3);
    /* MERIS */
    if (strcmp(EPR_ENVISAT_PRODUCT_MERIS, tmp) == 0) {
        if (product_id->sph_record == NULL) {
            product_id->sph_record = epr_read_sph(product_id);
            if (product_id->sph_record == NULL) {
                epr_set_err(e_err_file_read_error, "epr_set_param: wrong SPH");
                epr_free_string(tmp);
                return 0;
            }
        }

        field = epr_get_field(product_id->sph_record, "LINE_LENGTH");
        if (field == NULL) {
            epr_set_err(e_err_invalid_value,
                "epr_set_param: wrong SPH: unable to read LINE_LENGTH");
            epr_free_string(tmp);
            return 0;
        }

        line_length = ((uint*) field->elems)[0];
        if (line_length == 0) {
            epr_set_err(e_err_invalid_value,
                "epr_set_param: wrong SPH: LINE_LENGTH must be > 0");
            epr_free_string(tmp);
            return 0;
        }

        field = epr_get_field(product_id->sph_record, "LINES_PER_TIE_PT");
        if (field == NULL) {
            epr_set_err(e_err_invalid_value,
                "epr_set_param: wrong SPH: unable to read LINES_PER_TIE_PT");
            epr_free_string(tmp);
            return 0;
        }

        num_tie_points_across = ((uint*) field->elems)[0];
        if (num_tie_points_across == 0) {
            epr_set_err(e_err_invalid_value,
                "epr_set_param: wrong SPH: LINES_PER_TIE_PT must be > 0");
            epr_free_string(tmp);
            return 0;
        }

        param_elem = epr_create_param_elem("sceneRasterWidth", line_length);
        epr_add_ptr_array_elem(product_id->param_table, param_elem);

        ntpa = ((line_length - 1) / num_tie_points_across) + 1;

        param_elem = epr_create_param_elem("tiePointGridWidth", ntpa);
        epr_add_ptr_array_elem(product_id->param_table, param_elem);

    }

    /* AATSR does NOT have any dynamic parameters in DDDB */

    /* ASAR */
    else if ((strcmp(EPR_ENVISAT_PRODUCT_ASAR, epr_sub_string((char*)product_field->elems, 0, 3)) == 0) ||
             (strcmp(EPR_ENVISAT_PRODUCT_SAR, epr_sub_string((char*)product_field->elems, 0, 3)) == 0)) {

        field = epr_get_field(product_id->sph_record, "LINE_LENGTH");
        if (field == NULL) {
            epr_set_err(e_err_invalid_value,
                "epr_set_param: wrong SPH: unable to read LINE_LENGTH");
            epr_free_string(tmp);
            return 0;
        }

        line_length = ((uint*) field->elems)[0];
        if (line_length == 0) {
            epr_set_err(e_err_invalid_value,
                "epr_set_param: wrong SPH: LINE_LENGTH must be > 0");
            epr_free_string(tmp);
            return 0;
        }

        param_elem = epr_create_param_elem("sceneRasterWidth", line_length);
        epr_add_ptr_array_elem(product_id->param_table, param_elem);

        param_elem = epr_create_param_elem("tiePointGridWidth", EPR_ASAR_NUM_PER_POINT_ACROSS_LOCAT);
        epr_add_ptr_array_elem(product_id->param_table, param_elem);
    }

    epr_free_string(tmp);
    return 1;
}
