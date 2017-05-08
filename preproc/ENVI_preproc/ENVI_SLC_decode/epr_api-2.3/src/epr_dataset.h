/*
 * $Id: epr_dataset.h,v 1.1.1.1 2004-10-28 19:22:22 norman Exp $
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

#ifndef EPR_DATASET_H_INCL
#define EPR_DATASET_H_INCL

#ifdef __cplusplus
extern "C"
{
#endif

/**
 * @todo add docu here...
 */
EPR_SDatasetId* epr_create_dataset_id(EPR_SProductId* product_id,
                                      const EPR_SDSD* dsd,
                                      const char* dataset_name,
                                      const struct RecordDescriptor* record_info_ref,
                                      const char* dsd_name,
                                      const char* description);
/**
 * Release the memory allocated through a dataset ID.
 *
 * @param product_id the file identifier, if <code>NULL</code> the function
 *        immediately returns
 * @see epr_create_dataset_id
 */
void epr_free_dataset_id(EPR_SDatasetId* dataset_id);


EPR_SPtrArray* epr_create_dataset_ids(EPR_SProductId* product_id);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
/* #ifndef EPR_DATASET_H_INCL */

