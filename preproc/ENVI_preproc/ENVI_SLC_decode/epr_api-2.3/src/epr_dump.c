/*
 * $Id: epr_dump.c,v 1.2 2009-03-27 10:25:54 sabine Exp $
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


/****************************** RESULTS OUTPUT ******************************/

/**
 * Dumps the record content to stdout.
 *
 * @param record the pointer to the record to be printed out.
 */
void epr_dump_record(const EPR_SRecord* record)
{
    epr_print_record(record, stdout);
}

/**
 * Dumps the record content to an output FILE stream.
 *
 * @param record the pointer to the record to be printed out.
 * @param ostream the identifier of the output file.
 */
void epr_print_record(const EPR_SRecord* record, FILE* ostream)
{
    uint field_index;
    EPR_SField* field = NULL;

    epr_clear_err();

    for (field_index = 0; field_index < record->num_fields; field_index++)
    {
        field = record->fields[field_index];
        epr_print_field(field, ostream);
    }
}

/**
 * Dumps the field content to stdout.
 *
 * @param field the pointer to the field to be printed out
 */
void epr_dump_field(const EPR_SField* field)
{
    epr_print_field(field, stdout);
}

/**
 * Gets the field content to output FILE stream.
 *
 * @param field the pointer to the field to be printed out
 * @param ostream the identifier of the output file.
 */
void epr_print_field(const EPR_SField* field, FILE* ostream)
{
    uint i;

    epr_clear_err();

    fprintf(ostream, "%s = ", field->info->name);
    if (field->info->data_type_id == e_tid_string)
    {
		fprintf(ostream, "\"%s\"", (const char*) field->elems);
/*
        fprintf(ostream, "\"");
        for (i = 0; i < field->info->num_elems; i++)
        {
            fprintf(ostream, "%c", ((char*) field->elems)[i]);
        }
        fprintf(ostream, "\"");
*/
    }
    else if (field->info->data_type_id == e_tid_time)
    {
        EPR_STime* time = (EPR_STime*) field->elems;
        fprintf(ostream, "{d=%d, j=%d, m=%d}", time->days, time->seconds, time->microseconds);
    }
    else {
		if (field->info->num_elems > 1) {
			fprintf(ostream, "{");
		}
        for (i = 0; i < field->info->num_elems; i++)
        {
            if (i > 0)
                fprintf(ostream, ", ");
            switch (field->info->data_type_id)
            {
            case e_tid_uchar:
                fprintf(ostream, "%u", ((uchar*) field->elems)[i]);
                break;
            case e_tid_char:
                fprintf(ostream, "%d", ((char*) field->elems)[i]);
                break;
            case e_tid_ushort:
                fprintf(ostream, "%u", ((ushort*) field->elems)[i]);
                break;
            case e_tid_short:
                fprintf(ostream, "%d", ((short*) field->elems)[i]);
                break;
            case e_tid_uint:
                fprintf(ostream, "%u", ((uint*) field->elems)[i]);
                break;
            case e_tid_int:
                fprintf(ostream, "%d", ((int*) field->elems)[i]);
                break;
            case e_tid_float:
                fprintf(ostream, "%f", ((float*) field->elems)[i]);
                break;
            case e_tid_double:
                fprintf(ostream, "%f", ((double*) field->elems)[i]);
                break;
            default:
                fprintf(ostream, "<<unknown data type>>");
            }
        }
		if (field->info->num_elems > 1) {
	        fprintf(ostream, "}");
		}
    }
    fprintf(ostream, "\n");
}


/**
 * Dumps the element content to stdout.
 *
 * @param record the pointer to the element to be printed out.
 */
void epr_dump_element(const EPR_SRecord* record, uint field_index, uint element_index)
{
    epr_print_element(record, field_index, element_index, stdout);
}


/**
 * Dumps the element content to an output FILE stream..
 *
 * @param record the pointer to the element to be written out.
 */
void epr_print_element(const EPR_SRecord* record, uint field_index, uint element_index, FILE* ostream)
{
    EPR_SField* field = NULL;

    epr_clear_err();

    if (field_index >= record->num_fields)
    {
        epr_set_err(e_err_illegal_arg,
                    "epr_print_element: element_index too large");
        return;
    }

    field = record->fields[field_index];

    if (element_index >= field->info->num_elems)
    {
        epr_set_err(e_err_illegal_arg,
                    "epr_print_element: element_index too large");
        return;
    }

    fprintf(ostream, "%s [%d][%d] = ", field->info->name, field_index, element_index);
    if (field->info->data_type_id == e_tid_string)
    {
        fprintf(ostream, "\"");
        fprintf(ostream, "%c", ((char*) field->elems)[element_index]);
        fprintf(ostream, "\"");
    }
    else if (field->info->data_type_id == e_tid_time)
    {
        EPR_STime* time = (EPR_STime*) field->elems;
        fprintf(ostream, "{d=%d, j=%d, m=%d}", time->days, time->seconds, time->microseconds);
    }
    else {
        fprintf(ostream, "{ ");
        switch (field->info->data_type_id)
        {
            case e_tid_uchar:
                fprintf(ostream, "%u", ((uchar*) field->elems)[element_index]);
                break;
            case e_tid_char:
                fprintf(ostream, "%d", ((char*) field->elems)[element_index]);
                break;
            case e_tid_ushort:
                fprintf(ostream, "%u", ((ushort*) field->elems)[element_index]);
                break;
            case e_tid_short:
                fprintf(ostream, "%d", ((short*) field->elems)[element_index]);
                break;
            case e_tid_uint:
                fprintf(ostream, "%u", ((uint*) field->elems)[element_index]);
                break;
            case e_tid_int:
                fprintf(ostream, "%d", ((int*) field->elems)[element_index]);
                break;
            case e_tid_float:
                fprintf(ostream, "%f", ((float*) field->elems)[element_index]);
                break;
            case e_tid_double:
                fprintf(ostream, "%f", ((double*) field->elems)[element_index]);
                break;
            default:
                fprintf(ostream, "<<unknown data type>>");
        }
        fprintf(ostream, " }");
    }
    fprintf(ostream, " [%s]\n", field->info->unit);
}
