/*
 * $Id: epr_field.h,v 1.1.1.1 2004-10-28 19:22:23 norman Exp $
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

#ifndef EPR_FIELD_H_INCL
#define EPR_FIELD_H_INCL

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h> /* just to get the ANSI-C type FILE */

/**
 * The <code>EPR_FieldInfo</code> structure contains meta information
 * about a particular record field.
 */
struct EPR_FieldInfo
{
    /**
     * This field's name.
     */
    char* name;

    /**
     * This field's internal data type.
     */
    EPR_EDataTypeId data_type_id;

    /**
     * The number of data elements contained in this field (field-width).
     */
    uint num_elems;

    /**
     * This field's unit. Optional, can be NULL.
     */
    char* unit;

    /**
     * This field's description. Optional, can be NULL.
     */
    char* description;

    /**
     * The total size in bytes of all data elements of a field.
     * <code>tot_size</code> is a derived variable, it is computed at
     * runtime and not stored in the DSD-DB.
     */
    uint tot_size;
};

EPR_SFieldInfo* epr_create_field_info(EPR_EDataTypeId data_type_id, char* description, char* field_name, uint num_elems, uint num_bytes, uint more_count, char* unit);
EPR_SField* epr_create_field(EPR_SFieldInfo* field_info);
void epr_free_field_info(EPR_SFieldInfo* field_info);
void epr_free_field(EPR_SField* field);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
/* #ifndef EPR_FIELD_H_INCL */
