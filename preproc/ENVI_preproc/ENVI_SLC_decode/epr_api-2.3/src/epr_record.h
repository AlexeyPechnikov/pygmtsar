/*
 * $Id: epr_record.h,v 1.1.1.1 2004-10-28 19:22:23 norman Exp $
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

#ifndef EPR_RECORD_H_INCL
#define EPR_RECORD_H_INCL

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h> /* just to get the ANSI-C type FILE */

/**
 * The <code>EPR_RecordInfo</code> structure contains meta information
 * about a particular record.
 */
struct EPR_RecordInfo
{
    /**
     * The name of the dataset to which this record belongs to.
     */
    char* dataset_name;

    /**
     * The array of field info pointers.
     */
    EPR_SPtrArray* field_infos;
    /**
     * The total size in bytes of all data elements of all fields
     * of a record in a product file.
     * <code>tot_size</code> is a derived variable, it is computed at
     * runtime and not stored in the DSD-DB.
     */
    uint tot_size;
};

/**
 * Returns information about the structure of the records contained in
 * a dataset specified by the given <code>dataset_id</code>.
 *
 * @param dataset_id the the dataset identifier
 *
 * @return the the pointer for the record structure information.
 */
EPR_SRecordInfo* epr_get_record_info(EPR_SDatasetId* dataset_id);
EPR_SRecordInfo* epr_read_record_info(EPR_SProductId* product_id, EPR_SDatasetId* dataset_id);
void epr_read_sub_record_info(EPR_SProductId* product_id, FILE* db_file_istream, const char* parent_name, EPR_SPtrArray* field_infos);
EPR_SRecordInfo* epr_create_record_info(const char* dataset_name, EPR_SPtrArray* field_infos);
EPR_SRecord* epr_create_record_from_info(EPR_SRecordInfo* record_info);
char* epr_get_record_info_path(EPR_SProductId* product_id, const char* dataset_name);
void epr_free_record_info(EPR_SRecordInfo* record_info);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
/* #ifndef EPR_RECORD_H_INCL */
