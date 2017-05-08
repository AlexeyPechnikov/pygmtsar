/*
 * $Id: epr_dsd.h,v 1.1.1.1 2004-10-28 19:22:23 norman Exp $
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

#ifndef EPR_DSD_H_INCL
#define EPR_DSD_H_INCL

#ifdef __cplusplus
extern "C"
{
#endif

#include "epr_ptrarray.h"

#include <stdio.h> /* just to get the ANSI-C type FILE */

/**
 * Opens dsd for a dataset description,
 * obtained from an ENVISAT product file.
 *
 * @param dsd_index the number of dsd (zero-based), emrty dsd inclusive
 *
 * @return the the pointer at the dsd information.
 */
EPR_SDSD* epr_create_dsd(int dsd_index);

/**
 * Release the memory allocated through a dataset description,
 * obtained from an ENVISAT product file.
 *
 * @param dsd the file identifier, if <code>NULL</code> the function
 *        immediately returns zero.
 * @return zero for success, an error code otherwise
 * @see epr_free_dsd_id
 */
void epr_free_dsd(EPR_SDSD* dsd);


/**
 * Reads a dataset description from an ENVISAT product file.
 *
 * @param envisat_source_file the handle of the given ENVISAT product file,
 *        must not be <code>NULL</code>
 * @param pos number of the dataset description in ENVISAT product file,
 * @return a new dataset description or <code>NULL</code> if an error occured.
 */
EPR_SDSD* epr_read_each_dsd(FILE* envisat_source_file, int* pos);

/**
 * Reads all dataset descriptions from an ENVISAT product file.
 *
 * @param product_id the file identifier, if <code>NULL</code> the function
 *        immediately returns <code>NULL</code>.
 * @return an array of dataset descriptions or <code>NULL</code> if an error occured.
 */
EPR_SPtrArray* epr_read_all_dsds(EPR_SProductId* product_id);

/**
 * Finds the first dataset description from an ENVISAT product file.
 *
 * @param envisat_source_file the handle of the given ENVISAT product file,
 *        must not be <code>NULL</code>
 * @param sph_length [bytes] the length of SPH part from the given ENVISAT product file,
 *        must not be <code>NULL</code>
 * @return the offset to first founded dsd or <code>0</code> if not found.
 */
uint epr_find_first_dsd(FILE* envisat_source_file, uint sph_length);

int epr_detect_meris_iodd_version(EPR_SProductId* product_id);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
/* #ifndef EPR_DSD_H_INCL */
