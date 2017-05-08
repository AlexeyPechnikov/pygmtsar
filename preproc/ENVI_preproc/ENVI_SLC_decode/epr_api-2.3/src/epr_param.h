/*
 * $Id: epr_param.h,v 1.1.1.1 2004-10-28 19:22:23 norman Exp $
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

#ifndef EPR_PARAM_H_INCL
#define EPR_PARAM_H_INCL

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h> /* just to get the ANSI-C type FILE */

/**
 * The <code>EPR_ParamElem</code> structure contains meta information
 * how many time is a correspondent parameter represents.
 */
struct EPR_ParamElem
{
    /**
     * This is parameter name.
     */
    char* param_name;

    /**
     * The number of data elements represented of this parameter.
     */
    uint param_value;
};

EPR_SPtrArray* epr_create_param_table(void);
EPR_SParamElem* epr_create_param_elem(const char* param_name, int param_value);
int epr_set_dyn_dddb_params(EPR_SProductId* product_id);

void epr_free_param_table(EPR_SPtrArray* param_table);
void epr_free_param_elem(EPR_SParamElem* param_elem);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
/* #ifndef EPR_PARAM_H_INCL */
