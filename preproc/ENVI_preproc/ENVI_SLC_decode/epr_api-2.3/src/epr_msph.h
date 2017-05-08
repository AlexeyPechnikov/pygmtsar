/*
 * $Id: epr_msph.h,v 1.1.1.1 2004-10-28 19:22:23 norman Exp $
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

#ifndef EPR_MSPH_H_INCL
#define EPR_MSPH_H_INCL

#ifdef __cplusplus
extern "C" 
{
#endif

#include <stdio.h> /* just to get the ANSI-C type FILE */

/*void epr_read_mph(EPR_SProductId* product_id);*/
EPR_SRecord* epr_parse_header(const char* header_name, const char* ascii_source);
void epr_set_header_field_values(EPR_SRecord* record, EPR_SPtrArray* header_values);
uint epr_compare_param(EPR_SProductId* product_id);

void epr_parse_double_token(EPR_SPtrArray* header_values, char* token_value, uint* value_number, uint* l, EPR_EDataTypeId* tp);
void epr_parse_int_token(EPR_SPtrArray* header_values, char* token_value, uint* value_number, uint* l, EPR_EDataTypeId* tp);
#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* EPR_MSPH_H_INCL */
