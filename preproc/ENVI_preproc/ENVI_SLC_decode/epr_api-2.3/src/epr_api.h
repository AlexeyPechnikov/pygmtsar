/*
 * $Id: epr_api.h,v 1.3 2009-03-27 10:25:54 sabine Exp $
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

#ifndef EPR_API_H_INCL
#define EPR_API_H_INCL


#ifdef __cplusplus
extern "C"
{
#endif

/* to make the FILE structure available */
#include <stdio.h>

/* to make dynamic arrays available*/
#include "epr_ptrarray.h"

#define EPR_PRODUCT_API_NAME_STR         "ENVISAT Product Reader API"
#define EPR_PRODUCT_API_VERSION_STR      "2.3"

/* needed by Doxygen */
/** \mainpage
 * \htmlinclude doxygen_main_content.html
 */

/**
 * The <code>EPR_DataTypeId</code> enumeration lists all possible data
 * types for field elements in ENVISAT dataset records.
 */
enum EPR_DataTypeId
{
    /** The ID for unknown types. */
    e_tid_unknown = 0,
    /** An array of unsigned 8-bit integers, C type is <code>uchar*</code> */
    e_tid_uchar   = 1,
    /** An array of signed 8-bit integers, C type is <code>char*</code> */
    e_tid_char    = 2,
    /** An array of unsigned 16-bit integers, C type is <code>ushort*</code> */
    e_tid_ushort  = 3,
    /** An array of signed 16-bit integers, C type is <code>short*</code> */
    e_tid_short   = 4,
    /** An array of unsigned 32-bit integers, C type is <code>uint*</code> */
    e_tid_uint    = 5,
    /** An array of signed 32-bit integers, C type is <code>int*</code> */
    e_tid_int     = 6,
    /** An array of 32-bit floating point numbers, C type is <code>float*</code> */
    e_tid_float   = 7,
    /** An array of 64-bit floating point numbers, C type is <code>double*</code> */
    e_tid_double  = 8,
    /** A zero-terminated ASCII string, C type is <code>char*</code> */
    e_tid_string  = 11,
	/** An array of unsigned character, C type is <code>uchar*</code> */
    e_tid_spare  = 13,
    /** A time (MJD) structure, C type is <code>EPR_Time</code> */
    e_tid_time    = 21
};


/**
 * The <code>EPR_ErrCode</code> enumeration lists all possible error
 * codes for the ENVISAT product reader API.
 */
enum EPR_ErrCode
{
    /* Not an error */
    e_err_none                 =    0,

    /* Low level errors */
    e_err_null_pointer         =    1,
    e_err_illegal_arg          =    2,
    e_err_illegal_state        =    3,
    e_err_out_of_memory        =    4,
    e_err_index_out_of_range   =    5,
    e_err_illegal_conversion   =	6,
    e_err_illegal_data_type	   =	7,

    /* I/O errors */
    e_err_file_not_found       =  101,
    e_err_file_access_denied   =  102,
    e_err_file_read_error      =  103,
    e_err_file_write_error     =  104,
    e_err_file_open_failed     =  105,
    e_err_file_close_failed    =  106,

    /* API related errors */
    e_err_api_not_initialized  =  201,
    e_err_invalid_product_id   =  203,
    e_err_invalid_record	   =  204,
    e_err_invalid_band		   =  205,
    e_err_invalid_raster       =  206,
    e_err_invalid_dataset_name =  207,
    e_err_invalid_field_name   =  208,
	e_err_invalid_record_name  =  209,
	e_err_invalid_product_name =  210,
    e_err_invalid_band_name    =  211,
	e_err_invalid_data_format  =  212,
	e_err_invalid_value        =  213,
	e_err_invalid_keyword_name =  214,
	e_err_unknown_endian_order =  216,

    /* Bitmask term errors */
	e_err_flag_not_found       =  301,


    /* DDDB errors */
    e_err_invalid_ddbb_format  =  402
};


/**
 * The <code>EPR_LogLevel</code> enumeration lists possible log levels
 * for the ENVISAT product reader API.
 */
enum EPR_LogLevel
{
    e_log_debug   = -1,
    e_log_info    =  0,
    e_log_warning =  1,
    e_log_error   =  2
};

enum EPR_SampleModel
{
    e_smod_1OF1 = 0,
    e_smod_1OF2 = 1,
    e_smod_2OF2 = 2,
    e_smod_3TOI = 3,
    e_smod_2TOF = 4
};

enum EPR_ScalingMethod
{
    e_smid_non = 0,
    e_smid_lin = 1,
    e_smid_log = 2
};

struct EPR_ProductId;
struct EPR_DatasetId;
struct EPR_BandId;
struct EPR_Record;
struct EPR_RecordInfo;
struct EPR_Field;
struct EPR_FieldInfo;
struct EPR_ProductInfo;
struct EPR_DSD;
struct EPR_Raster;
struct EPR_DatasetRef;
struct EPR_Flag;
struct EPR_BandId;
struct EPR_ParamElem;
struct EPR_Time;

typedef enum   EPR_DataTypeId      EPR_EDataTypeId;
typedef enum   EPR_ErrCode         EPR_EErrCode;
typedef enum   EPR_LogLevel        EPR_ELogLevel;
typedef enum   EPR_SampleModel     EPR_ESampleModel;
typedef enum   EPR_ScalingMethod   EPR_EScalingMethod;
typedef struct EPR_ProductId       EPR_SProductId;
typedef struct EPR_DatasetId       EPR_SDatasetId;
typedef struct EPR_BandId		   EPR_SBandId;
typedef struct EPR_Record          EPR_SRecord;
typedef struct EPR_RecordInfo      EPR_SRecordInfo;
typedef struct EPR_Field           EPR_SField;
typedef struct EPR_FieldInfo       EPR_SFieldInfo;
typedef struct EPR_DSD             EPR_SDSD;
typedef struct EPR_Raster          EPR_SRaster;
typedef struct EPR_FlagDef         EPR_SFlagDef;
typedef struct EPR_ParamElem	   EPR_SParamElem;
typedef struct EPR_Time            EPR_STime;
typedef struct EPR_DatasetRef      EPR_SDatasetRef;
typedef struct EPR_BitmaskTerm     EPR_SBitmaskTerm;
typedef struct EPR_FlagSet         EPR_SFlagSet;
typedef void (*EPR_FErrHandler)(EPR_EErrCode err_code, const char* err_message);
typedef void (*EPR_FLogHandler)(EPR_ELogLevel log_level, const char* log_message);


typedef int            epr_boolean;
typedef unsigned char  uchar;
typedef unsigned short ushort;
typedef unsigned int   uint;
typedef unsigned long  ulong;


typedef int EPR_Magic;

#define EPR_MAGIC_PRODUCT_ID     0xCAFFEE64
#define EPR_MAGIC_DATASET_ID     0xEFEABDCA
#define EPR_MAGIC_BAND_ID        0xFEC21ABD
#define EPR_MAGIC_RECORD         0x7BABACAE
#define EPR_MAGIC_FIELD          0xBA0BABBA
#define EPR_MAGIC_RASTER         0x0BABA0EB
#define EPR_MAGIC_FLAG_DEF       0xCABA11AD

#define TRUE   1
#define FALSE  0

#define EPR_PRODUCT_ID_STRLEN    48


/*************************************************************************/
/******************************** STRUCTURES *****************************/
/*************************************************************************/

/**
 * The <code>EPR_ProductId</code> structure contains information
 * about an ENVISAT product file which has been opened with the
 * <code>epr_open_product()</code> function.
 *
 * @see epr_open_product
 */
struct EPR_ProductId
{
    /**
     * The magic number for this structure.
     * IMPORTANT: This must always be the first member of this structure.
     */
    EPR_Magic magic;

    /**
     * The file's path including the file name.
     */
    char* file_path;

    /**
     * The input stream as returned by the ANSI C <code>fopen</code>
     * function for the given file path.
     */
    FILE* istream;

    /**
     * The total size in bytes of the product file.
     */
    uint  tot_size;

    /**
     * The total width of product's scene raster in pixels.
     */
    uint  scene_width;

    /**
     * The total height of product's scene raster in pixels.
     */
    uint  scene_height;

    /**
     * The product identifier string obtained from the MPH
     * parameter 'PRODUCT'.
     * <p>The first 10 characters of this string identify the
     * the product type, e.g. "MER_1P__FR" for a MERIS Level 1b
     * full resolution product. The rest of the string decodes
     * product instance properties.
     */
    char id_string[EPR_PRODUCT_ID_STRLEN + 1];

    /**
     * The record representing the main product header (MPH).
     */
    EPR_SRecord* mph_record;

    /**
     * The record representing the specific product header (SPH).
     */
    EPR_SRecord* sph_record;

    /**
     * An array containing all (!) DSDs read from the product's
     * specific product header (SPH).
     */
    EPR_SPtrArray* dsd_array;

    /**
     * Cache for record infos. Contains all record infos read
     * from the database for this file so far.
     *
     * The reason for caching record infos on a per-file-base is that
     * a some record infos instances can contain file related content
     * such as the number of pixels in a measurecment dataset record
     * (MDSR).
     */
    EPR_SPtrArray* record_info_cache;

    /**
     * A table containing dynamic field info parameters.
     * Dynamic field info parameters are created at runtime because
     * the are derived from the product file contents and can
     * not be staically stored in the record info database.
     */
    EPR_SPtrArray* param_table;

    /**
     * Contains and array of all dataset IDs for the product (type EPR_SDatasetId*)
     */
    EPR_SPtrArray* dataset_ids;

    /**
     * Contains and array of all band IDs for the product (type EPR_SBandId*)
     */
    EPR_SPtrArray* band_ids;

	/**
     * For MERIS L1b and RR and FR to provide backward compatibility
     */
	int meris_iodd_version;
};




/**
 * The <code>EPR_DatasetId</code> structure contains information
 * about a dataset within an ENVISAT product file which has been opened with the
 * <code>epr_open_product()</code> API function.
 *
 * A new <code>EPR_DatasetId</code> instance can be obtained with the
 * <code>epr_get_dataset_id()</code> or <code>epr_get_dataset_id_at()</code> API functions.
 *
 * @see epr_open_product
 * @see epr_get_dataset_id
 * @see epr_get_dataset_id_at
 */
struct EPR_DatasetId
{
    /**
     * The magic number for this structure.
     * IMPORTANT: This must always be the first member of this structure.
     */
    EPR_Magic magic;

    /**
     * The ID of the product to which this dataset belongs to.
     */
    EPR_SProductId* product_id;

    /**
     * The name as presented to the user in a dsd selection dialog
     */
    char* dsd_name;

    /**
     * The dataset descriptor obtained from the current product.
     */
    const EPR_SDSD* dsd;

    /**
     * The name as presented to the user in a dataset selection dialog
     */
    char* dataset_name;

    /**
     * The record descriptor found in the DDDB for this dataset.
     */
	const struct RecordDescriptor* record_descriptor;

    /**
     * The record info which describes a record of this dataset.
     */
    EPR_SRecordInfo* record_info;


    /**
     * A short description of the band's contents
     */
    char* description;
};


/**
 * The <code>EPR_DSD</code> structure contains information
 * about the propertier of a dataset properties and its location
 * within an ENVISAT product file
 *
 * @see epr_read_each_dsd
 */
struct EPR_DSD
{
    /**
     * The magic number for this structure.
     * IMPORTANT: This must always be the first member of this structure.
     */
    EPR_Magic magic;

    /**
     * The index of this DSD (zero-based)
     */
    int index;

    /**
     * The dataset name.
     */
    char* ds_name;

    /**
     * The dataset type descriptor.
     */
    char* ds_type;

    /**
     * The filename in the DDDB with the description of this dataset.
     */
    char* filename;

    /**
     * The offset of dataset-information the product file.
     */
    uint ds_offset;

    /**
     * The size of dataset-information in dataset product file.
     */
    uint ds_size;

    /**
     * The number of dataset records for the given dataset name.
     */
    uint num_dsr;

    /**
     * The size of dataset record for the given dataset name.
     */
    uint dsr_size;
};


/**
 * The <code>EPR_Record</code> structure represents a
 * record instance read from an ENVISAT dataset.
 * A record is composed of multiple fields.
 *
 * @see EPR_Field
 */
struct EPR_Record
{
    /**
     * The magic number for this structure.
     * IMPORTANT: This must always be the first member of this structure.
     */
    EPR_Magic magic;

    /**
     * The corresponding record info for this record (a 'soft' pointer).
     */
    EPR_SRecordInfo* info;

    /**
     * The number of fields contained in this record.
     * The value is always equal <code>info->field_infos->length</code> and is
     * provided here for convenience only.
     */
    uint num_fields;

    /**
     * The record fields. An array of <code>EPR_Field*</code>
     * of length <code>info->num_fields</code>
     */
    EPR_SField** fields;
};

/**
 * Represents a field within a record. A field is composed of
 * one or more data elements of one of the types defined in the
 * in <code>field_info</code>.
 *
 * @see EPR_Record
 */
struct EPR_Field
{
    /**
     * The magic number for this structure.
     * IMPORTANT: This must always be the first member of this structure.
     */
    EPR_Magic magic;

    /**
     * The corresponding field info for this field.
     * <strong>supplierCardinality</strong> 1
     */
    EPR_SFieldInfo* info;

    /**
     * The elements of this field.
     *
     * In order to use the data, this member must be casted to one
     * of the following array data types:
     *
     * -# <code>unsigned char*</code> - array of unsigned 8-bit integer elements
     * -# <code>char*</code> - array of signed 8-bit integer fields
     * -# <code>unsigned short*</code> - array of unsigned 16-bit integer elements
     * -# <code>short*</code> - array of signed 16-bit integer elements
     * -# <code>unsigned int*</code> - array of unsigned 32-bit integer elements
     * -# <code>int*</code> - array of signed 32-bit integer elements
     * -# <code>float*</code> - array of signed 32-bit floating point elements
     * -# <code>double*</code> - array of signed 64-bit floating point elements
     * -# <code>EPR_STime*</code> - array of MJD elements
     *
     * Dedicated access routine are available
     * The element type is given by <code>info->data_type_id</code> and the array length by
     * <code>info->num_elems</code>.
     */
    void* elems;
};

/**
 * Represents a raster in which data will be stored.
 *
 * All 'size' parameter are in PIXEL.
 */
struct EPR_Raster
{
    /**
     * The magic number for this structure.
     * IMPORTANT: This must always be the first member of this structure.
     */
    EPR_Magic magic;

    /**
     * The data type of the band's pixel values.
     * <br>All corresponding with EPR_DataTypeId types are possible
     */
    EPR_EDataTypeId data_type;

    /**
     * The size in byte of a single element (sample) of this raster's buffer.
     */
    uint elem_size;

    /**
     * The width of the source .
     */
    uint source_width;

     /**
     * The height of the source.
     */
    uint source_height;

     /**
     * The sub-sampling for the across-track direction in pixel.
     */
    uint source_step_x;

     /**
     * The sub-sampling for the along-track direction in pixel.
     */
    uint source_step_y;

     /**
     * The width of the raster in pixel.
     * <br>raster_width  = (source_width  - 1) / source_step_x + 1
     */
    uint raster_width;

     /**
     * The height of the raster in pixel.
     * <br>raster_height = (source_height - 1) / source_step_y + 1
     */
    uint raster_height;

     /**
     * The elements of this raster.
     * <br>Its volume is <b>raster_width * raster_height * sizeof(data_type) in bytes</b>.
     */
    void* buffer;
};


/**
 * The <code>EPR_DatasetRef</code> structure represents the information from <code>dddb</code>
 * <br>with the reference to data name (in dddb), field-name and index
 * of the element in field-array, in which (by name) searchable values are located.
 * Example for the search for a scaling_offset information:
 * This information for the <code>reflec_10</code> is described with the <code>Scaling_Factor_GADS.22.10</code>
 * In <code>dataset_id</code> the searched ENVISAT product name (e.g. <code>MER_RR__2P</code>) is located.
 * <br>In the corresponding file (e.g. <code>/product/MER_RR__2P.dd</code>) the path,
 * how to find that information will be decribed.
 * In that file, in the field number <code>22</code> there is an information about the location
 * of the searched value in the ENVISAT product data.
 *
 * @see EPR_SDatasetId
 */
struct EPR_DatasetRef
{
    EPR_SDatasetId* dataset_id;
    int             field_index; /* -1 if not used */
    int             elem_index;  /* -1 if not used */
};

/**
 * Represents a flag-field within a flag-record.
 *
 */
struct EPR_FlagDef
{
    /**
     * The magic number for this structure.
     * IMPORTANT: This must always be the first member of this structure.
     */
    EPR_Magic magic;

    /**
     * The flag name.
     */
    char* name;

     /**
     * The bit mask describing this flag
     */
    uint bit_mask;

     /**
     * The flag description.
     */
    char* description;
};


/**
 * The <code>EPR_BandId</code> structure contains information
 * about a band within an ENVISAT product file which has been opened with the
 * <code>epr_open_product()</code> API function.
 *
 * A new <code>EPR_BandId</code> instance can be obtained with the
 * <code>epr_get_band_id()</code> API function.
 *
 * @see epr_open_product
 * @see epr_get_band_id
 */
struct EPR_BandId
{
    /**
     * The magic number for this structure.
     * IMPORTANT: This must always be the first member of this structure.
     */
    EPR_Magic magic;

    /**
     * The ID of the product to which this band belongs to.
     */
    EPR_SProductId* product_id;

    /**
     * The name as presented to the user in a band selection dialog
     * (also known as spectral subset)
     */
    char* band_name;

    /**
     * The (zero-based) spectral band index. -1 if this is not a spectral band.
     */
    int spectr_band_index;

    /**
     * The reference of the source dataset containing the raw data used to
     * create the band's pixel values. The external format used in the DDDB
     * is <code>MDS-name.field</code>, where <code>field</code> is a one-based
     * index (field=1 corresponds to the first field)
     */
     EPR_SDatasetRef dataset_ref;

    /**
     * The sample model operation applied to the source dataset for getting the
     * correct samples from the MDS (for example MERIS L2).
     * Possible values are:
     * - <code>*</code> --> no operation (direct copy)
     * - <code>1OF2</code>  --> first byte of 2-byte interleaved MDS
     * - <code>2OF2</code>  --> second byte of 2-byte interleaved MDS
     * - <code>0123</code>  --> combine 3-bytes interleaved to 4-byte integer
     */
     EPR_ESampleModel sample_model;

    /**
     * The data type of the band's pixel values. Possible values are:
     * - <code>*</code>     --> the datatype remains unchanged.
     * - <code>uint8_t</code>  --> 8-bit unsigned integer
     * - <code>uint32_t</code>  --> 32-bit unsigned integer
     * - <code>Float</code>  --> 32-bit IEEE floating point
     */
     EPR_EDataTypeId data_type;

    /**
     * The scaling method which must be applied to the raw source data in order
     * to get the 'real' pixel values in geo-physical units.
     * Possible values are:
     * - <code>*</code>            --> no scaling applied
     * - <code>Linear_Scale</code> --> linear scaling applied: y = offset + scale * x
     * - <code>Log_Scale</code>    --> logarithmic scaling applied: y = log10(offset + scale * x)
     */
     EPR_EScalingMethod scaling_method;

    /**
     * The scaling offset. Possible values are:<br>
     * - <code>*</code> --> no offset provided (implies scaling_method=*)
     * - <code>const</code> --> a floating point constant
     * - <code>GADS.field[.field2]</code> --> value is provided in
     * global annotation dataset with name <code>GADS</code> in field
     * <code>field</code>. Optionally a second element index
     *                  for multiple-element fields can be given too
     */
     float scaling_offset;

    /**
     * The scaling factor.
     * Possible values are:
     * - <code>*</code>                   --> no factor provided (implies scaling_method=*)
     * - <code>const</code>               --> a floating point constant
     * - <code>GADS.field[.field2]</code> --> value is provided in global annotation dataset with name
     *                                        <code>GADS</code> in field <code>field</code>. Optionally
     *                                        a second element index for multiple-element fields can be
     *                                        given too
     */
     float scaling_factor;

    /**
     * A bit-mask expression used to filter valid pixels. All others are set to zero.
     */
     char* bm_expr;

    /**
     * The flag coding is a list of EPR_SFlag instances. It determines each of the flags
     * used in this band (= flags dataset).
     * Each flag has a name, a bit-index and a description.
     */
    EPR_SPtrArray* flag_coding;

    /**
     * The geophysical unit for the band's pixel values
     */
    char* unit;

    /**
     * A short description of the band's contents
     */
    char* description;

    /**
     * If true (=1) lines will be mirrored (flipped) after read into a raster
     * in order to ensure a pixel ordering in raster X direction from
     * WEST to EAST.
     */
    epr_boolean lines_mirrored;
};

/**
 * Represents a binary time value field in ENVISAT records.
 *
 * <p> Refer to ENVISAT documentation for the exact definition of
 * this data type.
 */
struct EPR_Time
{
    int  days;
    uint seconds;
    uint microseconds;
};



/*************************************************************************/
/********************************* FUNCTIONS *****************************/
/*************************************************************************/

/*
 * ============================ (1) Initialisation ==========================
 */

/**
 * @defgroup INIT API Initialisation
 * @{
 */

/**
 * Initializes the ENVISAT product reader API.
 *
 *
 * @param log_level the log level. All logging messages with a log level lower
 *        than the given one, will be supressed
 * @param log_handler the log handler function pointer which
 *        will be used for logging, can be <code>NULL</code>,
 *        if logging shall be disabled
 * @param err_handler the new error handler (function pointer),
 *         can be <code>NULL</code>, if errors shall not be reported
 * @return zero for success, an error code otherwise
 *
 * @author Norman Fomferra
 */
int epr_init_api(EPR_ELogLevel   log_level,
                 EPR_FLogHandler log_handler,
                 EPR_FErrHandler err_handler);


/**
 * Closes the ENVISAT product reader API by releasing all
 * resources allocated by the API.
 *
 * @author Norman Fomferra
 */
void epr_close_api();
/** @} */


/*
 * ============================ (2) Logging ============================
 */

/**
 * @defgroup LOGGING Logging
 * @{
 */

/**
 * Sets the log level for the ENVISAT API. All logging
 * messages with a log level lower than the given one, will
 * be supressed, thus the log handler will not be called
 * for such messages.
 *
 * @param log_level the new log level. All logging messages with a log level lower
 *        than the given one, will be supressed
 * @return zero for success, an error code otherwise
 */
int epr_set_log_level(EPR_ELogLevel log_level);

/**
 * Sets the log handler for the ENVISAT API.
 *
 * @param log_handler the log handler function pointer which
 *        will be used for logging, can be NULL, if logging shall
 *        be disabled
 *
 * @see epr_log_message
 */
void epr_set_log_handler(EPR_FLogHandler log_handler);

/**
 * A default implementation for a logging function to be passed into
 * the <code>epr_init()</code> function. The function writes to
 * <code>stdout</code>, the format is: <i>log_level date time log_message</i>.
 *
 * @param log_level the log level
 * @param log_message the log message
 */
void epr_log_message(EPR_ELogLevel log_level, const char* log_message);

/** @} */

/*
 * ========================= (3) Error Handling ==========================
 */

/**
 * @defgroup ERROR Error Handling
 * @{
 */

/**
 * Sets the error handler for the ENVISAT API.
 *
 * @param err_handler the new error handler (function pointer),
 *         can be NULL, if errors shall not be reported
 */
void epr_set_err_handler(EPR_FErrHandler err_handler);

/**
 * Gets the error code of the error that occured during
 * the last API function call.
 *
 * @return the error code, <code>e_err_none</code> or zero if no error occured
 */
EPR_EErrCode epr_get_last_err_code();

/**
 * Gets the error message of the error that occured during
 * the last API function call.
 *
 * @return the error message, <code>NULL</code> if no error occured
 */
const char* epr_get_last_err_message();

/**
 * Clears the last error. After calling this function, calling
 * <code>epr_get_last_err_code</code> returns <code>e_err_none</code> or zero and
 * <code>epr_get_last_err_message</code> returns <code>NULL</code>.
 */
void epr_clear_err();

/** @} */

/*
 * ========================== (4) Input / Output ============================
 */

/**
 * @defgroup IO Input / Output
 * @{
 */

/*
 * ======================= (4.1) Product File Access ==========================
 */

/**
 * @defgroup ProductIO Product IO
 * @{
 */

/**
 * Opens the ENVISAT product file with the given file path,
 * <br>reads MPH, SPH and all DSDs, <br>organized the table with
 * parameter of line length and tie points number;
 * <br>returns a file identifier for the product.
 *
 * <p>The ENVISAT product reader API must be initialized before.
 *
 * @param product_file_path the path to the ENVISAT product file
 * @return the product identifier, or <code>NULL</code> if the file
 *         could not be opened. <code>epr_get_error_code()</code> should
 *         be called in this case in order to obtain the error code.
 */
EPR_SProductId* epr_open_product(const char* product_file_path);

/**
 * Closes the ENVISAT product file determined by the given product identifier.
 *
 * @param product_id the product identifier, if <code>NULL</code> the function
 *        immediately returns zero.
 * @return zero for success, an error code otherwise
 */
int epr_close_product(EPR_SProductId* product_id);
/** @} */

/** @} */

/*
 * ================= (4.2) Writing to a file or standard output =================
 */

/**
 * @ingroup IO
 * @defgroup WtFoSO Writing to a file or standard output
 * This group of functions is for writing an object to a file
 * or standard output.
 *
 * An object can be:
 * - record
 * - field
 * - field element
 *
 * If <code>FILE* istream</code> is given, the ASCII file will be outputed,
 * else printed to standard output device.
 *
 * <p>In case <i>record and/or field</i>:
 * @param record the record, must not be <code>NULL</code>
 * @param field the field, must not be <code>NULL</code>
 *
 * <p>In case <i>field element</i>:
 * @param record the record, must not be <code>NULL</code>
 * @param field_index the index of field in the given record
 * @param element_index the index of element in the given field
 *
 * @param ostream the identifier of the output file.
 */
/** @{ */
void epr_print_record(const EPR_SRecord* record, FILE* ostream);
void epr_print_field(const EPR_SField* field, FILE* ostream);
void epr_print_element(const EPR_SRecord* record, uint field_index, uint element_index, FILE* ostream);
void epr_dump_record(const EPR_SRecord* record);
void epr_dump_field(const EPR_SField* field);
void epr_dump_element(const EPR_SRecord* record, uint field_index, uint element_index);
/** @} */

/*
 * ======================= (5) Basic Data Access =========================
 */

/**
 * @defgroup DA Basic Data Access
 * @{
 */

/**
 * Gets the product's scene width in pixels.
 *
 * @param product_id the product identifier, must not be <code>NULL</code>
 * @return the product's total scene width in pixels, or <code>0</code>
 *         if an error occured.
 */
uint epr_get_scene_width(const EPR_SProductId* product_id);

/**
 * Gets the product's scene height in pixels.
 *
 * @param product_id the product identifier, must not be <code>NULL</code>
 * @return the product's total scene height in pixels, or <code>0</code>
 *         if an error occured.
 */
uint epr_get_scene_height(const EPR_SProductId* product_id);

/** @} */

/*
 * ============================ (5.1) Dataset ==============================
 */

/**
 * @ingroup DA
 * @defgroup DATASET Dataset Access
 * @{
 */

/**
 * Gets the number off all datasets contained in a product.
 *
 * @param product_id the product identifier, must not be <code>NULL</code>
 * @return the number off all dataset
 */
uint epr_get_num_datasets(EPR_SProductId* product_id);

/**
 * Gets the dataset_id at the specified position within the product.
 *
 * @param product_id the product identifier, must not be <code>NULL</code>
 * @param index the index identifying the position of the dataset, starting with 0,
 * must not be negative
 * @return the requested dataset_id
 */
EPR_SDatasetId* epr_get_dataset_id_at(EPR_SProductId* product_id, uint index);

/**
 * Gets the dataset_id coresponding to the specified dataset name.
 *
 * @param product_id the product identifier, must not be <code>NULL</code>
 * @param dataset_name the dataset name, must not be <code>NULL</code>
 * @return the requested dataset_id
 */
EPR_SDatasetId* epr_get_dataset_id(EPR_SProductId* product_id, const char* dataset_name);

/**
 * Gets the name of the dataset for the given dataset ID.
 *
 * @param dataset_id the dataset identifier, must not be <code>NULL</code>
 * @return the name of the dataset.
 */
const char* epr_get_dataset_name(EPR_SDatasetId* dataset_id);

/**
 * Gets the name of the dsd for the given dataset ID.
 *
 * @param dataset_id the dataset identifier, must not be <code>NULL</code>
 * @return the name of the dsd.
 */
const char* epr_get_dsd_name(const EPR_SDatasetId* dataset_id);

/**
 * Gets the MPH record from the given <code>product_id</code>.
 *
 * @param product_id the product identifier, must not be <code>NULL</code>
 * @return the MPH record or <code>NULL</code> if an error occured.
 */
EPR_SRecord* epr_get_mph(const EPR_SProductId* product_id);

/**
 * Gets the SPH record from the given <code>product_id</code>.
 *
 * @param product_id the product identifier, must not be <code>NULL</code>
 * @return the SPH record or <code>NULL</code> if an error occured.
 */
EPR_SRecord* epr_get_sph(const EPR_SProductId* product_id);

/**
 * Gets the dataset descriptor (DSD) for the dataset specified by <code>dataset_id</code>.
 *
 * @param dataset_id the dataset identifier, must not be <code>NULL</code>
 * @return the pointer at the dsd or <code>NULL</code> if an error occured.
 */
const EPR_SDSD* epr_get_dsd(const EPR_SDatasetId* dataset_id);

/**
 * Gets the number of records of the dataset specified by <code>dataset_id</code>.
 *
 * @param dataset_id the dataset identifier, must not be <code>NULL</code>
 * @return the number of records or <code>0</code> if an error occured.
 */
uint epr_get_num_records(const EPR_SDatasetId* dataset_id);


uint epr_get_num_dsds(const EPR_SProductId* product_id);
EPR_SDSD* epr_get_dsd_at(const EPR_SProductId* product_id, uint dsd_index);

/** @} */

/*
 * ================================= (5.2) Records ============================
 */

/**
 * @ingroup DA
 * @defgroup REC Record Access
 * @{
 */

/**
 * Creates a new, empty record with a structure compatible with the dataset specified
 * by dataset_id. Such a record is typically used in subsequent calls to epr_read_record.
 *
 * @param dataset_id the dataset identifier, must not be <code>NULL</code>
 * @return the new record instance
 *         or <code>NULL</code> if an error occured.
 */
EPR_SRecord* epr_create_record(EPR_SDatasetId* dataset_id);

/**
 * Reads a record of a dataset specified by dataset_id.
 * <p>
 * The record is identified through the given dataset identifier and the given
 * zero-based record index. In order to reduce memory reallocation, a
 * record (pre-) created by the function <code>epr_create_record</code>
 * can be passed to this function. Data is then read into this given record.
 * If no record (<code>NULL</code>) is given, the function initiates a new
 * one. In both cases, the record in which the data is read into will be
 * returned.
 *
 * @param dataset_id the dataset identifier, must not be <code>NULL</code>
 * @param record_index the zero-based record index
 * @param record a pre-created record to reduce memory reallocation,
 *        can be <code>NULL</code> to let the function allocate a new record
 * @return the record in which the data has been read into
 *         or <code>NULL</code> if an error occured.
 */
EPR_SRecord* epr_read_record(EPR_SDatasetId* dataset_id,
                             uint record_index,
                             EPR_SRecord* record);

/**
 * Frees the memory allocated through the given record.
 *
 * <p> After calling this function the given record pointer becomes
 * invalid and should not be used anymore.
 *
 */
void epr_free_record(EPR_SRecord* record);

/** @} */

/*
 * =========================== (5.3) Field Access =============================
 */

/**
 * @ingroup DA
 * @defgroup FA Field Access
 * @{
 */

/**
 * Gets a field from the given record.
 *
 * <p> The field is here identified through the given name.
 * <br>It contains the field info and all corresponding values.
 *
 * @param record the record identifier, must not be <code>NULL</code>
 * @param field_name the the name of required field, must not be <code>NULL</code>.
 * @return the field or <code>NULL</code> if an error occured.
 */
const EPR_SField* epr_get_field(const EPR_SRecord* record, const char* field_name);

/**
 * Gets the number of fields contained in the given record.
 *
 * @param record the record to be analysed, must not be <code>NULL</code>
 * @return the fields number or <code>0</code> if an error occured.
 */
uint epr_get_num_fields(const EPR_SRecord* record);

/**
 * Gets a field at the specified position within the record.
 *
 * @param record the record from the field shall be returned,
 *        must not be <code>NULL</code>
 * @param field_index the zero-based index (position within record) of the field
 * @return the field or <code>NULL</code> if an error occured.
 */
const EPR_SField* epr_get_field_at(const EPR_SRecord* record, uint field_index);

/**
 * Gets the unit of the field.
 *
 * @param field the field from which the unit shall be returned, must not be <code>NULL</code>
 * @return the field unit or <code>NULL</code> if an error occured.
 */
const char* epr_get_field_unit(const EPR_SField* field);

/**
 * Gets the description of the field.
 *
 * @param field field from which the description shall be returned, must not be <code>NULL</code>
 *
 * @return the field description or <code>NULL</code> if an error occured.
 */
const char* epr_get_field_description(const EPR_SField* field);

/**
 * Gets the number of elements of the field.
 *
 * @param field field to be analysed, must not be <code>NULL</code>
 *
 * @return the number of elements of the field or <code>0</code> if an error occured.
 */
uint epr_get_field_num_elems(const EPR_SField* field);

/**
 * Gets the name of the field.
 *
 * @param field field to be analysed, must not be <code>NULL</code>
 *
 * @return the field name or <code>NULL</code> if an error occured.
 */
const char* epr_get_field_name(const EPR_SField* field);

/**
 * Gets the type of the field.
 *
 * @param field field to be analysed, must not be <code>NULL</code>
 *
 * @return the field type or <code>0</code> if an error occured.
 */
EPR_EDataTypeId epr_get_field_type(const EPR_SField* field);

/** @} */

/*
 * ========================= (5.4) Single Element Access =========================
 */

/**
 * @ingroup DA
 * @defgroup FSEA Field Single Element Access
 * This group of functions is for getting the elements of a field as a typed value.
 * <br> Typed value means that the returned variable is of a certain variable type,
 * e.g. such as short or float.
 * <br> The type is located in the field info. See epr_get_field_type.
 * <br> One field must have one type only.
 *
 * @param field the field, must not be <code>NULL</code>
 * @param elem_index the zero-based index of element to be returned, must not be negative.
 *
 * @return the typed value from given field
 */
/*@{*/
char epr_get_field_elem_as_char(const EPR_SField* field, uint elem_index);
uchar epr_get_field_elem_as_uchar(const EPR_SField* field, uint elem_index);
short epr_get_field_elem_as_short(const EPR_SField* field, uint elem_index);
ushort epr_get_field_elem_as_ushort(const EPR_SField* field, uint elem_index);
int epr_get_field_elem_as_int(const EPR_SField* field, uint elem_index);
uint epr_get_field_elem_as_uint(const EPR_SField* field, uint elem_index);
float epr_get_field_elem_as_float(const EPR_SField* field, uint elem_index);
double epr_get_field_elem_as_double(const EPR_SField* field, uint elem_index);
const EPR_STime* epr_get_field_elem_as_mjd(const EPR_SField* field);
const char* epr_get_field_elem_as_str(const EPR_SField* field);
/*@}*/

/*
 * =========================== (5.5) Array Element Access =============================
 */

/**
 * @ingroup DA
 * @defgroup FAEA Field Array Element Access * This group of functions is for getting an array of field elements
 * of a certain data type.
 * <p>If the given field is not of the expected type, the functions
 * return <code>NULL</code>.
 *
 * @param field the field, must not be <code>NULL</code>
 * @return the data array of the expected type or  <code>NULL</code>
 *         if the field is of a different type
 */
/*@{*/
const char* epr_get_field_elems_char(const EPR_SField* field);
const uchar* epr_get_field_elems_uchar(const EPR_SField* field);
const short* epr_get_field_elems_short(const EPR_SField* field);
const ushort* epr_get_field_elems_ushort(const EPR_SField* field);
const int* epr_get_field_elems_int(const EPR_SField* field);
const uint* epr_get_field_elems_uint(const EPR_SField* field);
const float* epr_get_field_elems_float(const EPR_SField* field);
const double* epr_get_field_elems_double(const EPR_SField* field);
/*@}*/

/**
 * @ingroup DA
 * @defgroup CFE Copy Field Elems * This group of functions is for copying the data of the given field into the given
 * buffer of elements by selected type. The actual number of elements copied is the
 * minimum of the given number of elements (the buffer's size) and the actual number
 * of elements contained
 * in the field.
 * <p>If the actual field data type is not selected type, the function automatically
 * performs the conversion.
 *
 * @param field the field from which to copy the elements
 * @param buffer the buffer in which to copy the data
 * @param num_elems the number of elements in the given buffer
 * @return the actual number of elements copied
 */
/*@{*/
uint epr_copy_field_elems_as_ints(const EPR_SField* field, int* buffer, uint num_elems);
uint epr_copy_field_elems_as_uints(const EPR_SField* field, uint* buffer, uint num_elems);
uint epr_copy_field_elems_as_floats(const EPR_SField* field, float* buffer, uint num_elems);
uint epr_copy_field_elems_as_doubles(const EPR_SField* field, double* buffer, uint num_elems);
/*@}*/



/*
 * ======================== (6) Geophysical Data Access =========================
 */

/**
 * @defgroup GDA Geopysical Data Access
 * @{
 */

 /** @} */

/*
 * ================================== (6.1) Raster ===============================
 */

/**
 * @ingroup GDA
 * @defgroup RASTER Raster Data Access
 * @{
 */

/**
 * Creates a raster which is compatible with the data type contained in the band
 * identified by band_id. The created raster is used to read the data in it
 * (see epr_read_band_raster).
 * <p>
 * The raster is defined on the grid of the product, from which the data are read. Spatial
 * subsets and undersampling are possible) through the parameter of the function.
 * <p>
 * The concept of defining the raster is such: A certain portion of the ENVISAT product
 * will be read into the raster. This is called the source. The complete ENVISAT product can be
 * much greater than the source. One can move the raster over the complete ENVISAT product and
 * read in turn different parts - always of the size of the source - of it into the raster.
 * The source is specified by the parameter source_height and source_width.
 * <p>
 * A typical example is a processing in blocks. Lets say, a block
 * has 64x32 pixel. Then, my source has a width of 64 pixel and a height of 32 pixel. Another
 * example is a processing of complete image lines. Then, my source has a widths of the complete
 * product (for example 1121 for a MERIS RR product), and a height of 1). One can loop over all blocks
 * read into the raster and process it.
 * <p>
 * In addition, it is possible to defined a subsampling step for a raster. This means, that the
 * source is not read 1:1 into the raster, but that only every 2nd or 3rd pixel is read. This step
 * can be set differently for the across track (source_step_x) and along track (source_step_y) directions.
 *
 * @param band_id the band identifier. The raster will be compatible with the data type
 *        of that band; must not be <code>NULL</code>
 * @param source_width the width (across track dimension) of the source to be read into the raster. See text above.
 * @param source_height the height (along track dimension) of the source to be read into the raster. See text above.
 * @param source_step_x the subsampling step across track of the source when reading into the raster. See text above.
 * @param source_step_y the subsampling step along track of the source when reading into the raster. See text above.
 *
 * @return the new raster instance
 *         or <code>NULL</code> if an error occured.
 */
EPR_SRaster* epr_create_compatible_raster(EPR_SBandId* band_id,
                                          uint source_width,
                                          uint source_height,
                                          uint source_step_x,
                                          uint source_step_y);

/**
 * Creates a raster of the specified data type. This function can be used to create any type of raster,
 * e.g. for later use as a bit-mask.
 *
 * @param data_type the type of the data to stored in the raster, must not be <code>NULL</code>
 * @param source_width the width (across track dimension) of the source to be read into the raster. See description of epr_create_compatible_raster.
 * @param source_height the height (along track dimension) of the source to be read into the raster. See description of epr_create_compatible_raster.
 * @param source_step_x the subsampling step across track of the source when reading into the raster. See description of epr_create_compatible_raster.
 * @param source_step_y the subsampling step along track of the source when reading into the raster. See description of epr_create_compatible_raster.
 * @return the new raster instance
 *         or <code>NULL</code> if an error occured.
 */
EPR_SRaster* epr_create_raster(EPR_EDataTypeId data_type,
                               uint source_width,
                               uint source_height,
                               uint source_step_x,
                               uint source_step_y);


/**
 * Creates a raster to be used for reading bitmasks. The raster returned always is of type <code>byte</code>.
 *
 * @param source_width the width (across track dimension) of the source to be read into the raster. See description of epr_create_compatible_raster.
 * @param source_height the height (along track dimension) of the source to be read into the raster. See description of epr_create_compatible_raster.
 * @param source_step_x the subsampling step across track of the source when reading into the raster. See description of epr_create_compatible_raster.
 * @param source_step_y the subsampling step along track of the source when reading into the raster. See description of epr_create_compatible_raster.
 * @return the new raster instance
 *         or <code>NULL</code> if an error occured.
 */
EPR_SRaster* epr_create_bitmask_raster(uint source_width,
                                       uint source_height,
                                       uint source_step_x,
                                       uint source_step_y);

/**
 * Reads (geo-)physical values of the given band of the specified source-region.
 * <p> The source-region is a defined part of the whole ENVISAT product image, which shall be read into
 * a raster. In this routine the co-ordinates are specified, where the source-region to be read starts.
 * The dimension of the region and the sub-sampling are attributes of the raster into which the data are
 * read.
 *
 * @param band_id the identified of the band to be read into the raster.
 * @param offset_x across-track source co-ordinate in pixel co-ordinates (zero-based) of the upper right corner of the source-region
 * @param offset_y along-track source co-ordinate in pixel co-ordinates (zero-based) of the upper right corner of the source-region
 * @param raster the identifier to given raster information and raster buffer
 *
 * @return zero for success, and error code otherwise
 *
 * @see epr_create_compatible_raster
 * @see epr_create_rater
 */
int epr_read_band_raster(EPR_SBandId* band_id,
                         int offset_x,
                         int offset_y,
                         EPR_SRaster* raster);


/**
 * @todo 1 se/nf - doku
 */
uint epr_get_raster_elem_size(const EPR_SRaster* raster);

/**
 * @todo 1 se/nf - doku
 */
void* epr_get_raster_elem_addr(const EPR_SRaster* raster, uint offset);

/**
 * @todo 1 se/nf - doku
 */
void* epr_get_raster_pixel_addr(const EPR_SRaster* raster, uint x, uint y);

/**
 * @todo 1 se/nf - doku
 */
void* epr_get_raster_line_addr(const EPR_SRaster* raster, uint y);


/**
 * Gets the raster's scene width in pixels.
 *
 * @param raster the raster identifier, must not be <code>NULL</code>
 * @return the raster's total scene width in pixels, or <code>0</code>
 *         if an error occured.
 */
uint epr_get_raster_width(EPR_SRaster* raster);

/**
 * Gets the raster's scene height in pixels.
 *
 * @param raster the product identifier, must not be <code>NULL</code>
 * @return the raster's total scene height in pixels, or <code>0</code>
 *         if an error occured.
 */
uint epr_get_raster_height(EPR_SRaster* raster);



/**
 * Gets the number of all bands contained in a product.
 *
 * @param product_id the source product ID, must not be <code>NULL</code>
 * @return the number off all bands
 */
uint epr_get_num_bands(EPR_SProductId* product_id);

/**
 * Gets the band ID at the specified position within the product.
 *
 * @param product_id the source product ID, must not be <code>NULL</code>
 * @param index the index identifying the position of the band, starting with 0,
 * must not be negative
 * @return the requested band ID, or <code>NULL</code> if not found
 */
EPR_SBandId* epr_get_band_id_at(EPR_SProductId* product_id, uint index);

/**
 * Gets the band ID corresponding to the specified name.
 *
 * @param product_id the source product ID, must not be <code>NULL</code>
 * @param band_name the name of the band, must not be <code>NULL</code>
 * @return the requested band ID, or <code>NULL</code> if not found
 */
EPR_SBandId* epr_get_band_id(EPR_SProductId* product_id, const char* band_name);

/**
 * Gets the name of the band for the given band ID.
 *
 * @param band_id the band identifier, must not be <code>NULL</code>
 * @return the name of the band.
 */
const char* epr_get_band_name(EPR_SBandId* band_id);

/**
 * Release the memory allocated through a raster.
 *
 * @param raster the raster to be released.
 */
void epr_free_raster(EPR_SRaster* raster);

/** @} */

/*
 * ============================ (6.2) Single Pixel Access ========================
 */

/**
 * @ingroup GDA
 * @defgroup SPA Single Pixel Access
 * @{
 */

/**
 * This group of functions is for getting the values of the elements of a raster
 * (i.e. pixel) in a type-safe way. <br>
 *
 * @param raster the raster which contains the pixel, must not be <code>NULL</code>
 * @param x the (zero-based) X co-ordinate of the pixel
 * @param y the (zero-based) Y co-ordinate of the pixel
 *
 * @return the typed value at the given co-ordinate.
 */
uint epr_get_pixel_as_uint(const EPR_SRaster* raster, int x, int y);
int epr_get_pixel_as_int(const EPR_SRaster* raster, int x, int y);
float epr_get_pixel_as_float(const EPR_SRaster* raster, int x, int y);
double epr_get_pixel_as_double(const EPR_SRaster* raster, int x, int y);
/*@}*/

/*
 * ================================= (7) Bitmasks ==========================
 */

/**
 * @defgroup BM Bitmask
 * @{
 */

/**
 * Calculates a bit-mask, composed of flags of the given product and combined as described in the
 * given bit-mask expression, for the a certain dimension and sub-sampling as defined in the
 * given raster.
 * <p>
 *
 * @param product_id Identifier of the ENVISAT product for which the bit-mask shall be created.
 *                 This is used by the function to retreive the needed flags.
 * @param bm_expr A string holding the logical expression for the defintion of the bit-mask.
 *                 In a bit-mask expression, any number of the flag-names (found in the DDDB) can
 *                 be composed with "(", ")", "NOT", "AND", "OR". Valid bit-mask expression are for example: <br>
 *                 "flags.LAND OR flags.CLOUD" or "NOT flags.WATER AND flags.TURBID_S".
 * @param offset_x across-track co-ordinate in pixel co-ordinates (zero-based) of the upper right corner of the source-region
 * @param offset_y along-track co-ordinate in pixel co-ordinates (zero-based) of the upper right corner of the source-region
 * @param raster the raster for the bit-mask. The data type of the raster must be either e_tid_uchar or e_tid_char.
 *
 * @return zero for success, an error code otherwise
 *
 * @see create_band_raster
 */
int epr_read_bitmask_raster(EPR_SProductId* product_id,
                            const char* bm_expr,
                            int offset_x,
                            int offset_y,
                            EPR_SRaster* raster);

/** @} */

/*
 * ================================= (8) Utility functions ==========================
 */

/**
 * @defgroup UTILS Utility Functions
 * @{
 */

/**
 * Gets the size in bytes for an element of the given data type.
 */
uint epr_get_data_type_size(EPR_EDataTypeId data_type_id);

/**
 * Gets the 'C' data type string for the given data type.
 */
const char* epr_data_type_id_to_str(EPR_EDataTypeId data_type_id);

/** @} */

#ifdef __cplusplus
} /* extern "C" */
#endif
#endif /* #ifndef EPR_API_H_INCL */
