/* provides structures to read SAR tapes*/

/*
include files modified from the rceos.c programs
written by C. Tomassini & F. Lorenna

other format information from:
from CERS (RAW) CCT format specifications STD-TM#92-767F
   Canada Centre for Remote Sensing (CCRS)
   Surveys, Mapping and Remote Sensing Sector
   Energy, Mines and Resources Canada

   Table 6.1.2.2 "SARLEADER" FILE POINTER RECORD CONTENTS
   page 6.

	R. J. Mellors 
	July 1997, IGPP-SIO

from esa annex A (Document ER-IS-EPO-GS-5902.I)
     Issue 2.1.1

        Paul F. Jamason
        25-FEB-1997, IGPP-SIO

*/

struct sarleader_dss_binary {
	int	record_seq_no;
	char	record_subtype_code1;
	char	record_type_code1;
	char	record_subtype_code2;
	char	record_subtype_code3;
	int	record_length;
};

/* ccrs raw data set summary record format */
#define SARLEADER_DSS_RCS "%4c%4c%16c%32c%32c%16c%16c%16c%16c%16c%16c"\
"%16c%16c%16c%16c%16c%16c%16c%16c%8c%8c%16c%16c%16c%4c%4c%16c%32c%8c%8c"\
"%8c%8c%8c%8c%8c%16c%2c%16c%16c%16c%16c%16c%16c%16c%16c%16c%16c%16c%8c%8c"\
"%16c%16c%16c%4c%4c%32c%8c%12c%16c%16c%16c%32c%16c%16c%4c%16c%32c%16c%32c"\
"%8c%8c%16c%8c%8c%32c%32c%32c%16c%16c%16c%16c%16c%16c%32c%32c%16c%16c%16c"\
"%32c%16c%16c%16c%16c%16c%16c%16c%8c%8c%16c%16c%16c%16c%16c%16c%16c%16c%8c"\
"%4c%4c%16c%16c%16c%32c%16c%16c%16c%24c%24c%24c%546c"

/* ccrs raw data set summary corresponding log file output */
#define SARLEADER_DSS_RVL(SP)\
(SP)->dss_rec_seq_num,\
(SP)->chan_ind,\
(SP)->reserved1 ,\
(SP)->scene_number ,\
(SP)->input_scene_center_time,\
(SP)->spare1,\
(SP)->center_lat,\
(SP)->center_long,\
(SP)->center_heading,\
(SP)->ellipsoid_designator,\
(SP)->ellipsoid_semimajor_axis,\
(SP)->ellipsoid_semiminor_axis,\
(SP)->earth_constant,\
(SP)->spare2,\
(SP)->ellipsoid_j2,\
(SP)->ellipsoid_j3,\
(SP)->ellipsoid_j4,\
(SP)->spare,\
(SP)->reserved_new,\
(SP)->scene_centre_line_number,\
(SP)->scene_centre_pixel_number,\
(SP)->scene_length,\
(SP)->scene_width,\
(SP)->spare3,\
(SP)->nchan,\
(SP)->spare4,\
(SP)->mission_identifier,\
(SP)->sensor_id_and_mode,\
(SP)->orbit_number,\
(SP)->lat_nadir_center,\
(SP)->long_nadir_center,\
(SP)->heading_nadir_center,\
(SP)->clock_angel,\
(SP)->incidence_angle_center,\
(SP)->radar_freq,\
(SP)->radar_wavelength,\
(SP)->motion_compensation,\
(SP)->range_pulse_code_specifier,\
(SP)->range_pulse_amplitude_const,\
(SP)->range_pulse_amplitude_lin,\
(SP)->range_pulse_amplitude_quad,\
(SP)->range_pulse_amplitude_cube,\
(SP)->range_pulse_amplitude_quart,\
(SP)->range_pulse_phase_const,\
(SP)->range_pulse_phase_lin,\
(SP)->range_pulse_phase_quad,\
(SP)->range_pulse_phase_cube,\
(SP)->range_pulse_phase_quart,\
(SP)->chirp_extraction_index,\
(SP)->spare5,\
(SP)->sampling,\
(SP)->range_gate_early_edge_start_image,\
(SP)->range_pulse_length,\
(SP)->reserved2,\
(SP)->range_compressed_flag,\
(SP)->reserved3,\
(SP)->quantisation_in_bits,\
(SP)->quantizer_descriptor,\
(SP)->dc_bias_i,\
(SP)->dc_bias_q,\
(SP)->gain_imbalance,\
(SP)->spare6,\
(SP)->reserved4,\
(SP)->antenna_mech_bor,\
(SP)->reserved5,\
(SP)->nominal_prf,\
(SP)->reserved6,\
(SP)->satelite_encoded_binary_time,\
(SP)->satelite_clock_time,\
(SP)->satelite_clock_increment,\
(SP)->spare7,\
(SP)->processing_facility_identifier,\
(SP)->processing_system_id,\
(SP)->processing_version_id,\
(SP)->reserved7,\
(SP)->product_type_id,\
(SP)->alg_id,\
(SP)->nlooks_az,\
(SP)->neff_looks_range,\
(SP)->bandwidth_look_az,\
(SP)->bandwidth_look_range,\
(SP)->total_look_bandwidth_az,\
(SP)->total_look_bandwidth_range,\
(SP)->w_func_designator_az,\
(SP)->w_func_designator_range,\
(SP)->data_input_source,\
(SP)->nom_res_3db_range,\
(SP)->nom_res_az,\
(SP)->reserved8,\
(SP)->a_track_dop_freq_const_early_image,\
(SP)->a_track_dop_freq_lin_early_image,\
(SP)->a_track_dop_freq_quad_early_image,\
(SP)->spare8,\
(SP)->c_track_dop_freq_const_early_image,\
(SP)->c_track_dop_freq_lin_early_image,\
(SP)->c_track_dop_freq_quad_early_image,\
(SP)->time_direction_along_pixel,\
(SP)->time_direction_along_line,\
(SP)->a_track_dop_freq_rate_const_early_image,\
(SP)->a_track_dop_freq_rate_lin_early_image,\
(SP)->a_track_dop_freq_rate_quad_early_image,\
(SP)->spare9,\
(SP)->c_track_dop_freq_rate_const_early_image,\
(SP)->c_track_dop_freq_rate_lin_early_image,\
(SP)->c_track_dop_freq_rate_quad_early_image,\
(SP)->spare10,\
(SP)->line_content_indicator,\
(SP)->clut_lock_flag,\
(SP)->autofocussing_flag,\
(SP)->line_spacing,\
(SP)->pixel_spacing_range,\
(SP)->range_compression_designator,\
(SP)->spare11,\
(SP)->zero_dop_range_time_f_pixel,\
(SP)->zero_dop_range_time_c_pixel,\
(SP)->zero_dop_range_time_l_pixel,\
(SP)->zero_dop_az_time_f_pixel,\
(SP)->zero_dop_az_time_c_pixel,\
(SP)->zero_dop_az_time_l_pixel,\
(SP)->gec_local_use_segment

struct sarleader_dss {
	char    dss_rec_seq_num[4];   	/*dss record sequence number (1)*/
	char    chan_ind[4];            /*sar channel indicator (1)*/
	char    reserved1[16] ;         /* scene identifier*/
	char    scene_number[32] ;
	char    input_scene_center_time[32];
	char    spare1[16];
	char    center_lat[16];
	char    center_long[16];
	char    center_heading[16];
	char    ellipsoid_designator[16];
	char    ellipsoid_semimajor_axis[16];
	char    ellipsoid_semiminor_axis[16];
	char    earth_constant[16];
	char    spare2[16];
	char    ellipsoid_j2[16];
	char    ellipsoid_j3[16];
	char    ellipsoid_j4[16];
	char    spare[16];
	char    reserved_new[16];
	char    scene_centre_line_number[8];
	char    scene_centre_pixel_number[8];
	char    scene_length[16];
	char    scene_width[16];
	char    spare3[16];
	char    nchan[4];
	char    spare4[4];
	char    mission_identifier[16];
	char    sensor_id_and_mode[32];
	char    orbit_number[8];
	char    lat_nadir_center[8];
	char    long_nadir_center[8];
	char    heading_nadir_center[8];
	char    clock_angel[8];
	char    incidence_angle_center[8];
	char    radar_freq[8];
	char    radar_wavelength[16];
	char    motion_compensation[2];
	char    range_pulse_code_specifier[16];
	char    range_pulse_amplitude_const[16];
	char    range_pulse_amplitude_lin[16];
	char    range_pulse_amplitude_quad[16];
	char    range_pulse_amplitude_cube[16];
	char    range_pulse_amplitude_quart[16];
	char    range_pulse_phase_const[16];
	char    range_pulse_phase_lin[16];
	char    range_pulse_phase_quad[16];
	char    range_pulse_phase_cube[16];
	char    range_pulse_phase_quart[16];
	char    chirp_extraction_index[8];
	char    spare5[8];
	char    sampling[16];
	char    range_gate_early_edge_start_image[16];
	char    range_pulse_length[16];
	char    reserved2[4];
	char    range_compressed_flag[4];
	char    reserved3[32];
	char    quantisation_in_bits[8];
	char    quantizer_descriptor[12];
	char    dc_bias_i[16];
	char    dc_bias_q[16];
	char    gain_imbalance[16];
	char    spare6[32];
	char    reserved4[16];
	char    antenna_mech_bor[16];
	char    reserved5[4];
	char    nominal_prf[16];
	char    reserved6[32];
	char    satelite_encoded_binary_time[16];
	char    satelite_clock_time[32];
	char    satelite_clock_increment[8];
	char    spare7[8];
	char    processing_facility_identifier[16];
	char    processing_system_id[8];
	char    processing_version_id[8];
	char    reserved7[32];
	char    product_type_id[32];
	char    alg_id[32];
	char    nlooks_az[16];
	char    neff_looks_range[16];
	char    bandwidth_look_az[16];
	char    bandwidth_look_range[16];
	char    total_look_bandwidth_az[16];
	char    total_look_bandwidth_range[16];
	char    w_func_designator_az[32];
	char    w_func_designator_range[32];
	char    data_input_source[16];
	char    nom_res_3db_range[16];
	char    nom_res_az[16];
	char    reserved8[32];
	char    a_track_dop_freq_const_early_image[16];
	char    a_track_dop_freq_lin_early_image[16];
	char    a_track_dop_freq_quad_early_image[16];
	char    spare8[16];
	char    c_track_dop_freq_const_early_image[16];
	char    c_track_dop_freq_lin_early_image[16];
	char    c_track_dop_freq_quad_early_image[16];
	char    time_direction_along_pixel[8];
	char    time_direction_along_line[8];
	char    a_track_dop_freq_rate_const_early_image[16];
	char    a_track_dop_freq_rate_lin_early_image[16];
	char    a_track_dop_freq_rate_quad_early_image[16];
	char    spare9[16];
	char    c_track_dop_freq_rate_const_early_image[16];
	char    c_track_dop_freq_rate_lin_early_image[16];
	char    c_track_dop_freq_rate_quad_early_image[16];
	char    spare10[16];
	char    line_content_indicator[8];
	char    clut_lock_flag[4];
	char    autofocussing_flag[4];
	char    line_spacing[16];
	char    pixel_spacing_range[16];
	char    range_compression_designator[16];
	char    spare11[32];
	char    zero_dop_range_time_f_pixel[16];
	char    zero_dop_range_time_c_pixel[16];
	char    zero_dop_range_time_l_pixel[16];
	char    zero_dop_az_time_f_pixel[25];
	char    zero_dop_az_time_c_pixel[24];
	char    zero_dop_az_time_l_pixel[24];
	char    gec_local_use_segment[546];
} ;

/* dpaf raw data set summary record format */
#define SARLEADER_DPAF_DSS_RCS "%4c%4c%16c%32c%32c%16c%16c%16c%16c%16c%16c"\
"%16c%16c%16c%16c%16c%16c%16c%16c%8c%8c%16c%16c%16c%4c%4c%16c%32c%8c%8c"\
"%8c%8c%8c%8c%8c%16c%2c%16c%16c%16c%16c%16c%16c%16c%16c%16c%16c%16c%8c%8c"\
"%16c%16c%16c%4c%4c%32c%8c%12c%16c%16c%16c%32c%16c%16c%4c%16c%32c%16c%32c"\
"%8c%8c%16c%8c%8c%32c%32c%32c%16c%16c%16c%16c%16c%16c%32c%32c%16c%16c%16c"\
"%32c%16c%16c%16c%16c%16c%16c%16c%8c%8c%16c%16c%16c%16c%16c%16c%16c%16c%8c"\
"%4c%4c%16c%16c%16c%32c%16c%16c%16c%24c%24c%24c"

/* dpaf raw data set summary corresponding log file output */
#define SARLEADER_DPAF_DSS_RVL(SP)\
(SP)->dss_rec_seq_num,\
(SP)->chan_ind,\
(SP)->reserved1 ,\
(SP)->scene_number ,\
(SP)->input_scene_center_time,\
(SP)->spare1,\
(SP)->center_lat,\
(SP)->center_long,\
(SP)->center_heading,\
(SP)->ellipsoid_designator,\
(SP)->ellipsoid_semimajor_axis,\
(SP)->ellipsoid_semiminor_axis,\
(SP)->earth_constant,\
(SP)->spare2,\
(SP)->ellipsoid_j2,\
(SP)->ellipsoid_j3,\
(SP)->ellipsoid_j4,\
(SP)->spare,\
(SP)->reserved_new,\
(SP)->scene_centre_line_number,\
(SP)->scene_centre_pixel_number,\
(SP)->scene_length,\
(SP)->scene_width,\
(SP)->spare3,\
(SP)->nchan,\
(SP)->spare4,\
(SP)->mission_identifier,\
(SP)->sensor_id_and_mode,\
(SP)->orbit_number,\
(SP)->lat_nadir_center,\
(SP)->long_nadir_center,\
(SP)->heading_nadir_center,\
(SP)->clock_angel,\
(SP)->incidence_angle_center,\
(SP)->radar_freq,\
(SP)->radar_wavelength,\
(SP)->motion_compensation,\
(SP)->range_pulse_code_specifier,\
(SP)->range_pulse_amplitude_const,\
(SP)->range_pulse_amplitude_lin,\
(SP)->range_pulse_amplitude_quad,\
(SP)->range_pulse_amplitude_cube,\
(SP)->range_pulse_amplitude_quart,\
(SP)->range_pulse_phase_const,\
(SP)->range_pulse_phase_lin,\
(SP)->range_pulse_phase_quad,\
(SP)->range_pulse_phase_cube,\
(SP)->range_pulse_phase_quart,\
(SP)->chirp_extraction_index,\
(SP)->spare5,\
(SP)->sampling,\
(SP)->range_gate_early_edge_start_image,\
(SP)->range_pulse_length,\
(SP)->reserved2,\
(SP)->range_compressed_flag,\
(SP)->reserved3,\
(SP)->quantisation_in_bits,\
(SP)->quantizer_descriptor,\
(SP)->dc_bias_i,\
(SP)->dc_bias_q,\
(SP)->gain_imbalance,\
(SP)->spare6,\
(SP)->reserved4,\
(SP)->antenna_mech_bor,\
(SP)->reserved5,\
(SP)->nominal_prf,\
(SP)->reserved6,\
(SP)->satelite_encoded_binary_time,\
(SP)->satelite_clock_time,\
(SP)->satelite_clock_increment,\
(SP)->spare7,\
(SP)->processing_facility_identifier,\
(SP)->processing_system_id,\
(SP)->processing_version_id,\
(SP)->reserved7,\
(SP)->product_type_id,\
(SP)->alg_id,\
(SP)->nlooks_az,\
(SP)->neff_looks_range,\
(SP)->bandwidth_look_az,\
(SP)->bandwidth_look_range,\
(SP)->total_look_bandwidth_az,\
(SP)->total_look_bandwidth_range,\
(SP)->w_func_designator_az,\
(SP)->w_func_designator_range,\
(SP)->data_input_source,\
(SP)->nom_res_3db_range,\
(SP)->nom_res_az,\
(SP)->reserved8,\
(SP)->a_track_dop_freq_const_early_image,\
(SP)->a_track_dop_freq_lin_early_image,\
(SP)->a_track_dop_freq_quad_early_image,\
(SP)->spare8,\
(SP)->c_track_dop_freq_const_early_image,\
(SP)->c_track_dop_freq_lin_early_image,\
(SP)->c_track_dop_freq_quad_early_image,\
(SP)->time_direction_along_pixel,\
(SP)->time_direction_along_line,\
(SP)->a_track_dop_freq_rate_const_early_image,\
(SP)->a_track_dop_freq_rate_lin_early_image,\
(SP)->a_track_dop_freq_rate_quad_early_image,\
(SP)->spare9,\
(SP)->c_track_dop_freq_rate_const_early_image,\
(SP)->c_track_dop_freq_rate_lin_early_image,\
(SP)->c_track_dop_freq_rate_quad_early_image,\
(SP)->spare10,\
(SP)->line_content_indicator,\
(SP)->clut_lock_flag,\
(SP)->autofocussing_flag,\
(SP)->line_spacing,\
(SP)->pixel_spacing_range,\
(SP)->range_compression_designator,\
(SP)->spare11,\
(SP)->zero_dop_range_time_f_pixel,\
(SP)->zero_dop_range_time_c_pixel,\
(SP)->zero_dop_range_time_l_pixel,\
(SP)->zero_dop_az_time_f_pixel,\
(SP)->zero_dop_az_time_c_pixel,\
(SP)->zero_dop_az_time_l_pixel

struct sarleader_dpaf_dss {
	char    dss_rec_seq_num[4];   	/*dss record sequence number (1)*/
	char    chan_ind[4];            /*sar channel indicator (1)*/
	char    reserved1[16] ;         /* scene identifier*/
	char    scene_number[32] ;
	char    input_scene_center_time[32];
	char    spare1[16];
	char    center_lat[16];
	char    center_long[16];
	char    center_heading[16];
	char    ellipsoid_designator[16];
	char    ellipsoid_semimajor_axis[16];
	char    ellipsoid_semiminor_axis[16];
	char    earth_constant[16];
	char    spare2[16];
	char    ellipsoid_j2[16];
	char    ellipsoid_j3[16];
	char    ellipsoid_j4[16];
	char    spare[16];
	char    reserved_new[16];
	char    scene_centre_line_number[8];
	char    scene_centre_pixel_number[8];
	char    scene_length[16];
	char    scene_width[16];
	char    spare3[16];
	char    nchan[4];
	char    spare4[4];
	char    mission_identifier[16];
	char    sensor_id_and_mode[32];
	char    orbit_number[8];
	char    lat_nadir_center[8];
	char    long_nadir_center[8];
	char    heading_nadir_center[8];
	char    clock_angel[8];
	char    incidence_angle_center[8];
	char    radar_freq[8];
	char    radar_wavelength[16];
	char    motion_compensation[2];
	char    range_pulse_code_specifier[16];
	char    range_pulse_amplitude_const[16];
	char    range_pulse_amplitude_lin[16];
	char    range_pulse_amplitude_quad[16];
	char    range_pulse_amplitude_cube[16];
	char    range_pulse_amplitude_quart[16];
	char    range_pulse_phase_const[16];
	char    range_pulse_phase_lin[16];
	char    range_pulse_phase_quad[16];
	char    range_pulse_phase_cube[16];
	char    range_pulse_phase_quart[16];
	char    chirp_extraction_index[8];
	char    spare5[8];
	char    sampling[16];
	char    range_gate_early_edge_start_image[16];
	char    range_pulse_length[16];
	char    reserved2[4];
	char    range_compressed_flag[4];
	char    reserved3[32];
	char    quantisation_in_bits[8];
	char    quantizer_descriptor[12];
	char    dc_bias_i[16];
	char    dc_bias_q[16];
	char    gain_imbalance[16];
	char    spare6[32];
	char    reserved4[16];
	char    antenna_mech_bor[16];
	char    reserved5[4];
	char    nominal_prf[16];
	char    reserved6[32];
	char    satelite_encoded_binary_time[16];
	char    satelite_clock_time[32];
	char    satelite_clock_increment[8];
	char    spare7[8];
	char    processing_facility_identifier[16];
	char    processing_system_id[8];
	char    processing_version_id[8];
	char    reserved7[32];
	char    product_type_id[32];
	char    alg_id[32];
	char    nlooks_az[16];
	char    neff_looks_range[16];
	char    bandwidth_look_az[16];
	char    bandwidth_look_range[16];
	char    total_look_bandwidth_az[16];
	char    total_look_bandwidth_range[16];
	char    w_func_designator_az[32];
	char    w_func_designator_range[32];
	char    data_input_source[16];
	char    nom_res_3db_range[16];
	char    nom_res_az[16];
	char    reserved8[32];
	char    a_track_dop_freq_const_early_image[16];
	char    a_track_dop_freq_lin_early_image[16];
	char    a_track_dop_freq_quad_early_image[16];
	char    spare8[16];
	char    c_track_dop_freq_const_early_image[16];
	char    c_track_dop_freq_lin_early_image[16];
	char    c_track_dop_freq_quad_early_image[16];
	char    time_direction_along_pixel[8];
	char    time_direction_along_line[8];
	char    a_track_dop_freq_rate_const_early_image[16];
	char    a_track_dop_freq_rate_lin_early_image[16];
	char    a_track_dop_freq_rate_quad_early_image[16];
	char    spare9[16];
	char    c_track_dop_freq_rate_const_early_image[16];
	char    c_track_dop_freq_rate_lin_early_image[16];
	char    c_track_dop_freq_rate_quad_early_image[16];
	char    spare10[16];
	char    line_content_indicator[8];
	char    clut_lock_flag[4];
	char    autofocussing_flag[4];
	char    line_spacing[16];
	char    pixel_spacing_range[16];
	char    range_compression_designator[16];
	char    spare11[32];
	char    zero_dop_range_time_f_pixel[16];
	char    zero_dop_range_time_c_pixel[16];
	char    zero_dop_range_time_l_pixel[16];
	char    zero_dop_az_time_f_pixel[24];
	char    zero_dop_az_time_c_pixel[24];
	char    zero_dop_az_time_l_pixel[24];
} ;


