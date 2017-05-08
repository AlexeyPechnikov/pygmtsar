/*
 * $Id: epr_swap.c,v 1.1.1.1 2004-10-28 19:22:23 norman Exp $
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


/*
 * Function: byte_swap_short.c
*/
/**
 *
 * Swaps bytes within NUMBER_OF_SWAPS two-byte words,
 *   starting at address BUFFER.
 *
 * @param buffer the one element typed buffer
 * to convert for a little endian order machine
 *
 * @param number_of_swaps number of elements to convert
 *
 */
void byte_swap_short(short *buffer, uint number_of_swaps)
{
   short* temp = buffer;
   uint swap_loop;

   for (swap_loop = 0, temp = buffer; swap_loop < number_of_swaps; swap_loop++, temp++) {
      *temp = (short)(((*temp & 0x00ff) << 8) |
                      ((*temp & 0xff00) >> 8));
   }
}


/*
   Function: byte_swap_int.c
*/
/**
 *
 *  Swaps bytes within NUMBER_OF_SWAPS four-byte words,
 *     starting at address BUFFER.
 *
 *
 */
void byte_swap_int(int *buffer, uint number_of_swaps)
{
   int *temp = buffer;
   uint swap_loop;

   for (swap_loop = 0, temp = buffer; swap_loop < number_of_swaps; swap_loop++, temp++) {
      *temp = ((*temp & 0x000000ff) << 24) |
              ((*temp & 0x0000ff00) << 8)  |
              ((*temp & 0x00ff0000) >> 8)  |
              ((*temp & 0xff000000) >> 24);
   }
}


/*
   Function: byte_swap_short.c
*/
/**
 *
 * Swaps bytes within NUMBER_OF_SWAPS two-byte words,
 *   starting at address BUFFER.
 *
 * @param buffer the one element typed buffer
 * to convert for a little endian order machine
 *
 * @param number_of_swaps number of elements to convert
 *
 */
void byte_swap_ushort(ushort* buffer, uint number_of_swaps)
{
   byte_swap_short((short*) buffer, number_of_swaps);
}

/*
 *  Function: byte_swap_uint.c
 */
/**
 *
 * Swaps bytes within NUMBER_OF_SWAPS four-byte words,
 *     starting at address BUFFER.
 *
 * @param buffer the one element typed buffer
 * to convert for a little endian order machine
 *
 * @param number_of_swaps number of elements to convert
 *
 */
void byte_swap_uint(uint* buffer, uint number_of_swaps)
{
   byte_swap_int((int*) buffer, number_of_swaps);
}

/*
 *  Function: byte_swap_int.c
 */
/**
 *
 * Swaps bytes within NUMBER_OF_SWAPS four-byte words,
 *     starting at address BUFFER.
 *
 * @param buffer the one element typed buffer
 * to convert for a little endian order machine
 *
 * @param number_of_swaps number of elements to convert
 *
 */
void byte_swap_float(float* buffer, uint number_of_swaps)
{
   byte_swap_int((int*) buffer, number_of_swaps);
}

/**
 * A boolean value indicating whether this code run's on a
 * little endian order machine or not.
 * <p><code>1</code> stands for little endian (LE),
 * <code>0</code> stands for big endian (BE).
 */
/*
   Function:    epr_is_little_endian_order
   Access:      public API
   Changelog:   2002/02/04 nf nitial version
 */
/**
 * Returns a oolean value indicating whether this code run's on a
 * little endian order machine or not.
 * <p><code>1</code> stands for little endian (LE), <code>0/code> otherwise
 */
int epr_is_little_endian_order()
{
    uint le_value = EPR_LE_MAGIC_NUMBER;
    return (((uchar*)(&le_value))[0] == EPR_LE_MAGIC_BYTE_0)
        && (((uchar*)(&le_value))[1] == EPR_LE_MAGIC_BYTE_1)
        && (((uchar*)(&le_value))[2] == EPR_LE_MAGIC_BYTE_2)
        && (((uchar*)(&le_value))[3] == EPR_LE_MAGIC_BYTE_3);
}

/*
   Function:    epr_is_big_endian_order
   Access:      public API
   Changelog:   2002/02/04  nf nitial version
 */
/**
 * Returns a oolean value indicating whether this code run's on a
 * little endian order machine or not.
 * <p><code>1</code> stands for little endian (BE), <code>0/code> otherwise
 */
int epr_is_big_endian_order()
{
    uint be_value = EPR_BE_MAGIC_NUMBER;
    return (((uchar*)(&be_value))[0] == EPR_LE_MAGIC_BYTE_0)
        && (((uchar*)(&be_value))[1] == EPR_LE_MAGIC_BYTE_1)
        && (((uchar*)(&be_value))[2] == EPR_LE_MAGIC_BYTE_2)
        && (((uchar*)(&be_value))[3] == EPR_LE_MAGIC_BYTE_3);
}


/*
   Function:    epr_swap_endian_order
   Access:      public API
   Changelog:   2002/02/04  mp nitial version
 */
/**
 * Converts bytes for a little endian order machine
 *
 * @param field the pointer at data reading in
 *
 */
void epr_swap_endian_order(const EPR_SField* field)
{
    switch (field->info->data_type_id) {
        case e_tid_uchar:
        case e_tid_char:
        case e_tid_string:
            /* no conversion required */
            break;
        case e_tid_time:
            byte_swap_uint((uint*)field->elems, 3);
            break;
        case e_tid_spare:
            /* no conversion required */
            break;
        case e_tid_ushort:
            byte_swap_ushort((ushort*) field->elems, field->info->num_elems);
            break;
        case e_tid_short:
            byte_swap_short((short*) field->elems, field->info->num_elems);
            break;
        case e_tid_uint:
            byte_swap_uint((uint*) field->elems, field->info->num_elems);
            break;
        case e_tid_int:
            byte_swap_int((int*) field->elems, field->info->num_elems);
            break;
        case e_tid_float:
            byte_swap_float((float*) field->elems, field->info->num_elems);
            break;
        case e_tid_double:
            epr_set_err(e_err_invalid_data_format,
                    "epr_swap_endian_order: DOUBLE type was not yet processed");
            break;
        default:
            epr_set_err(e_err_invalid_data_format,
                    "epr_swap_endian_order: unknown data type");
    }
}
