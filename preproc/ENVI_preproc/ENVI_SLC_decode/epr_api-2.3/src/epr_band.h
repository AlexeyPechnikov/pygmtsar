/*
 * $Id: epr_band.h,v 1.1.1.1 2004-10-28 19:22:22 norman Exp $
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

#ifndef EPR_BAND_H_INCL
#define EPR_BAND_H_INCL

#ifdef __cplusplus
extern "C"
{
#endif

#include "epr_ptrarray.h"

#include <stdio.h> /* just to get the ANSI-C type FILE */


/**
 * Converts the given string into a scaling method identifier.
 *
 * @param str the string to be converted.
 * @return the scaling method identifier represented by the given string.
 *         If the string is equal of '*' the value
 *         <code>e_non_smid</code> is returned.
 */
EPR_EScalingMethod epr_str_to_scaling_method(const char* str);


/**
 * Converts the given string into a sample offset identifier.
 *
 * @param str the string to be converted.
 * @return the sample offset identifier represented by the given string.
 *         If the string is equal of '*' the value
 *         <code>e_none_samoff</code> is returned.
 */
EPR_ESampleModel epr_str_to_sample_offset(const char* str);


/**
 * Gets the dataset_id, field_index and elem_index
 *
 * @param product_id the the product file identifier
 * @param str the string with the name, separator ('.') and indexes.
 * @return the dataset_id, field_index and elem_index (-1 if no).
 *    <code>NULL</code> if correspondent dataset name was not found.
 */
EPR_SDatasetRef epr_get_ref_struct(EPR_SProductId* product_id, const char* str);


/**
 * Gets the scaling factor by the given dataset_id, field_index, elem_index
 *
 * @param product_id the the product file identifier
 * @param str the string with the name, separator ('.') and indexes.
 * @return the dataset_id, field_index and elem_index (-1 if no).
 *    <code>NULL</code> if correspondent dataset name was not found.
 */
float epr_get_scaling_factor(EPR_SProductId* product_id,  const char* str);
float epr_get_scaling_params(EPR_SProductId* product_id,  const char* str);

/**
 * Reads the measurement data and converts its in physical values.
 *
 * @param band_id the information about properties and quantities of ENVISAT data.
 * @param offset_x X-coordinate in pixel co-ordinates (zero-based) of the upper right corner raster to search
 * @param offset_y Y-coordinate in pixel co-ordinates (zero-based) of the upper right corner raster to search
 * @param raster the instance to the buffer information was used
 *
 * @return zero for success, and error code otherwise
 */
int epr_read_band_measurement_data(EPR_SBandId* band_id, int offset_x, int offset_y, EPR_SRaster* raster);

/**
 * Reads the annotation data and converts its in physical values.
 *
 * @param band_id the information about properties and quantities of ENVISAT data.
 * @param offset_x X-coordinate in pixel co-ordinates (zero-based) of the upper right corner raster to search
 * @param offset_y Y-coordinate in pixel co-ordinates (zero-based) of the upper right corner raster to search
 * @param raster the instance to the buffer information was used
 *
 * @return zero for success, and error code otherwise
 */
int epr_read_band_annotation_data(EPR_SBandId* band_id,
                       int offset_x,
                       int offset_y,
                       EPR_SRaster* raster);


/**
 * This group of functions is for scaling the field element for a physical measurement values.
 * <br> The type is located in the field info.
 * <br> One field must have one type only.
 *
 * @param sourceArray the sourse array identifier (to be scaled)
 * @param band_id the band ID with the information about the field's physical properties
 * @param xo [PIXEL] X-coordinate (0-bazed) of the upper right corner raster to search
 * @param raster_width [PIXEL] the width of the raster is been research
 * @param s_x X-step to get the next raster to search
 * @param raster_buffer [BYTE] the memory buffer to save information was scaled
 * @param raster_pos shows the point of filling of the array raster_buffer
 *
 */
/*@{*/
void decode_line_uchar_1_of_1_to_float   (void* sourceArray, EPR_SBandId* band_id, int xo, int raster_width, int s_x, void* raster_buffer, int raster_pos);
void decode_line_ushort_1_of_1_to_float  (void* sourceArray, EPR_SBandId* band_id, int xo, int raster_width, int s_x, void* raster_buffer, int raster_pos);
void decode_line_short_1_of_1_to_float   (void* sourceArray, EPR_SBandId* band_id, int xo, int raster_width, int s_x, void* raster_buffer, int raster_pos);
void decode_line_short_1_of_2_to_float   (void* sourceArray, EPR_SBandId* band_id, int xo, int raster_width, int s_x, void* raster_buffer, int raster_pos);
void decode_line_short_2_of_2_to_float   (void* sourceArray, EPR_SBandId* band_id, int xo, int raster_width, int s_x, void* raster_buffer, int raster_pos);
void decode_line_char_1_of_1_to_float    (void* sourceArray, EPR_SBandId* band_id, int xo, int raster_width, int s_x, void* raster_buffer, int raster_pos);
void decode_line_uchar_1_of_1_to_uchar   (void* sourceArray, EPR_SBandId* band_id, int xo, int raster_width, int s_x, void* raster_buffer, int raster_pos);
void decode_line_uchar_1_of_2_to_uchar   (void* sourceArray, EPR_SBandId* band_id, int xo, int raster_width, int s_x, void* raster_buffer, int raster_pos);
void decode_line_uchar_2_of_2_to_uchar   (void* sourceArray, EPR_SBandId* band_id, int xo, int raster_width, int s_x, void* raster_buffer, int raster_pos);
void decode_line_ushort_1_of_1_to_ushort (void* sourceArray, EPR_SBandId* band_id, int xo, int raster_width, int s_x, void* raster_buffer, int raster_pos);
void decode_line_uchar_1_of_2_to_float   (void* sourceArray, EPR_SBandId* band_id, int xo, int raster_width, int s_x, void* raster_buffer, int raster_pos);
void decode_line_uchar_2_of_2_to_float   (void* sourceArray, EPR_SBandId* band_id, int xo, int raster_width, int s_x, void* raster_buffer, int raster_pos);
void decode_line_uchar_2_to_f_to_float   (void* sourceArray, EPR_SBandId* band_id, int xo, int raster_width, int s_x, void* raster_buffer, int raster_pos);
void decode_line_uchar_3_to_i_to_uint   (void* sourceArray, EPR_SBandId* band_id, int xo, int raster_width, int s_x, void* raster_buffer, int raster_pos);
/*@}*/

/**
 * This group of functions is for scaling the field element for a physical annotation values.
 * <br> The type is located in the field info.
 * <br> One field must have one type only.
 *
 * @param sourceArray the sourse array identifier (to be scaled)
 * @param band_id the band ID with the information about the field's physical properties
 * @param raster_buffer [BYTE] the memory buffer to save information was scaled
 * @param nel number of element to scale
 *
 */
/*@{*/
void transform_array_short_to_float (void* sourceArray, EPR_SBandId* band_id, float* raster_buffer, uint nel);
void transform_array_ushort_to_float(void* sourceArray, EPR_SBandId* band_id, float* raster_buffer, uint nel);
void transform_array_int_to_float  (void* sourceArray, EPR_SBandId* band_id, float* raster_buffer, uint nel);
void transform_array_uint_to_float (void* sourceArray, EPR_SBandId* band_id, float* raster_buffer, uint nel);
/*@}*/

/**
 * This group of functions is for mirroring the scaled line of a physical MERIS values.
 *
 * @param raster_buffer [BYTE] the memory buffer to be Y-mirrored
 * @param raster_width [PIXEL] the width of the raster is been Y-mirrored
 * @param raster_height [PIXEL] the height of the raster is been Y-mirrored
 *
 */
/*@{*/
void mirror_float_array  (float*  raster_buffer, uint raster_width, uint raster_height);
void mirror_uchar_array  (uchar*  raster_buffer, uint raster_width, uint raster_height);
void mirror_ushort_array (ushort* raster_buffer, uint raster_width, uint raster_height);
void mirror_uint_array  (uint*  raster_buffer, uint raster_width, uint raster_height);
/*@}*/

/**
 * Two dimenzional interpolation
 *
 * @param wi the interpolation point location in [0,1] in "horizontal" direction
 * @param wj the interpolation point location in [0,1] in "vertical" direction
 * @param x00 the first point in "horizontal" direction
 * @param x10 the second point in "horizontal" direction
 * @param x01 the first point in "vertical" direction
 * @param x11 the second point in "vertical" direction
 *
 * @return float interpolated value
 */
float epr_interpolate2D(float wi, float wj, float x00, float x10, float x01, float x11);

/**
 * Computes the physical values for the annotation data.
 *
 * @param sa_beg the float array of tie points "before" Y-coordinate of the point to search
 * @param sa_end the float array of tie points "after" Y-coordinate of the point to search
 * @param samples_per_tie_pt the "distance" between two neighbour tie point (in scan-line direction)
 * @param num_elems number of elements in one tie point scan-line
 * @param band_id the information about properties and quantities of ENVISAT data.
 * @param xo [PIXEL] X-coordinate (0-bazed) of the upper right corner raster to search
 * @param y_mod [PIXEL] relativ location of the point is been researched (in fly direction)
 * @param raster_width [PIXEL] the width of the raster is been researched
 * @param s_x [PIXEL] X-step to get the next point (in source coordinates) to search
 * @param raster_buffer the float user array to be filled with physical values
 * @param raster_pos the actual "filled" position in raster_buffer-array
 *
 */
void decode_tiepoint_band(float* sa_beg,
						  float* sa_end,
						  uint samples_per_tie_pt,
						  uint num_elems,
						  EPR_SBandId* band_id,
						  int xo,
						  float scan_offset_x,
						  float y_mod,
						  int raster_width,
						  int s_x,
						  float* raster_buffer,
						  int raster_pos);
typedef void (*EPR_FLineDecoder)(void* sourceArray,
                    EPR_SBandId* band_id,
                    int xo,
                    int raster_width,
                    int s_x,
                    void* raster_buffer,
                    int raster_pos);

/**
 * Selects the line decode function, depended on measurement data type.
 */
EPR_FLineDecoder select_line_decode_function(EPR_EDataTypeId band_daty, EPR_ESampleModel band_smod, EPR_EDataTypeId daty_id);

typedef void (*EPR_FArrayTransformer)(void* sourceArray,
                                     EPR_SBandId* band_id,
                                     float* raster_buffer,
                                     uint nel);

/**
 * Selects the transform array function, dependent on annotation data type.
 */
EPR_FArrayTransformer select_transform_array_function(EPR_EDataTypeId band_daty, EPR_EDataTypeId daty_id);

/**
 * Masks the band information out.
 * The band information will be masked dependent on bit mask filter for the same
 * selected area (described in raster).
 *
 * @param raster selected and physically processed the ENVISAT product data band information
 * @param bm_raster selected the ENVISAT flag bit mask filter
 */
void epr_zero_invalid_pixels(EPR_SRaster* raster, EPR_SRaster* bm_raster);

/**
 * Release the memory allocated through a band ID.
 *
 * @param band_id the band identifier to be released.
 */
void epr_free_band_id(EPR_SBandId* band_id);



EPR_SPtrArray* epr_create_band_ids(EPR_SProductId* product_id);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
/* #ifndef EPR_BAND_H_INCL */
