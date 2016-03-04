/* provides structures to read SAR tapes*/
/* modified from the rceos programs by
 C. Tomassini & F. Lorenna */

/*
also from:
 from CERS (RAW) CCT format specifications STD-TM#92-767F
   Canada Centre for Remote Sensing (CCRS)
   Surveys, Mapping and Remote Sensing Sector
   Energy, Mines and Resources Canada

	R. J. Mellors 
	July 1997, IGPP-SIO
*/

#define SARDATA_HEADER_WCS "*********** SARDATA **********\n"\
"n_record			==>	%.6s\n"\
"record_length			==>	%.6s\n"\
"blank1				==>	%.24s\n"\
"n_bits_per_sample		==>	%.4s\n"\
"n_samples_data_group		==>	%.4s\n"\
"n_bytes_per_data_group		==>	%.4s\n"\
"justification			==>	%.4s\n"\
"n_channels			==>	%.4s\n"\
"n_lines			==>	%.8s\n"\
"n_left_border_pixels		==>	%.4s\n"\
"n_data_groups			==>	%.8s\n"\
"n_right_border_pixels		==>	%.4s\n"\
"n_top_border_lines		==>	%.4s\n"\
"bottom_border_lines		==>	%.4s\n"\
"interleaving			==>	%.4s\n"\
"n_phys_records			==>	%.2s\n"\
"n_phy_records_multi		==>	%.2s\n"\
"n_bytes_prefix_data		==>	%.4s\n"\
"n_bytes_SARdata_record		==>	%.8s\n"\
"n_bytes_suffix_data		==>	%.4s\n"\
"blank2				==>	%.48s\n"\
"blank3				==>	%.28s\n"\
"blank4				==>	%.32s\n"\
"data_format_type_id		==>	%.28s\n"\
"data_format_type_code		==>	%.4s\n"\
"n_left_fill_bits_pixel		==>	%.4s\n"\
"n_right_fill_bits_pixel	==>	%.4s\n"\
"max_data_range_pixel		==>	%.8s\n"\
"blank_filler			==>	%.10s\n"

#define SARDATA_HEADER_RCS "%6c%6c%24c%4c%4c%4c%4c%4c%8c%4c%8c%4c%4c%4c%4c%2c%2c%4c%8c%4c%48c%28c%32c%28c%4c%4c%4c%8c%10196c"

#define SARDATA_HEADER_RVL(SP)\
(SP)->n_records,\
(SP)->record_length,\
(SP)->blank1,\
(SP)->n_bits_per_sample,\
(SP)->n_samples_data_group,\
(SP)->n_bytes_per_data_group,\
(SP)->justification,\
(SP)->n_channels,\
(SP)->n_lines,\
(SP)->n_left_border_pixels,\
(SP)->n_data_groups,\
(SP)->n_right_border_pixels,\
(SP)->n_top_border_lines,\
(SP)->bottom_border_lines,\
(SP)->interleaving,\
(SP)->n_phys_records,\
(SP)->n_phy_records_multi,\
(SP)->n_bytes_prefix_data,\
(SP)->n_bytes_SARdata_record,\
(SP)->n_bytes_suffix_data,\
(SP)->blank2,\
(SP)->blank3,\
(SP)->blank4,\
(SP)->data_format_type_id,\
(SP)->data_format_type_code,\
(SP)->n_left_fill_bits_pixel,\
(SP)->n_right_fill_bits_pixel,\
(SP)->max_data_range_pixel,\
(SP)->blank_filler

struct sardata_header {
	char n_records[7];
	char record_length[7];
	char blank1[25];
	char n_bits_per_sample[5];
	char n_samples_data_group[5];
	char n_bytes_per_data_group[5];
	char justification[5];
	char n_channels[5];
	char n_lines[9];
	char n_left_border_pixels[5];
	char n_data_groups[9];
	char n_right_border_pixels[5];
	char n_top_border_lines[5];
	char bottom_border_lines[5];
	char interleaving[5];
	char n_phys_records[3];
	char n_phy_records_multi[3];
	char n_bytes_prefix_data[5]; 
 	char n_bytes_SARdata_record[9];
	char n_bytes_suffix_data[5];
	char blank2[49];
	char blank3[29];
	char blank4[33];
	char data_format_type_id[29];
	char data_format_type_code[5];
 	char n_left_fill_bits_pixel[5];
 	char n_right_fill_bits_pixel[5];
	char max_data_range_pixel[5];
	char blank_filler[10197];
};

#define SARDATA_REC_WCS "*********** SARDATA **********\n"\
"line_no		==>	%1d\n"\
"rec_index		==>	%1d\n"\
"left_fill_pixels	==>	%1d\n"\
"data_pixels		==>	%1d\n"\
"right_fill_pixels	==>	%1d\n"\
"sensor_update_flag			==>	%1d\n"\
"year			==>	%1d\n"\
"day_of_year		==>	%1d\n"\
"msecs_of_day		==>	%1d\n"\
"chan_ind		==>	%1d\n"\
"chan_code		==>	%1d\n"\
"trans_polar		==>	%1d\n"\
"rec_polar		==>	%1d\n"\
"prf			==>	%1d\n"\
"spare			==>	%1d\n"\
"range_flag		==>	%1d\n"\
"chirp_type		==>	%1d\n"\
"chirp_length		==>	%1d\n"\
"chirp_const_coeff	==>	%1d\n"\
"chirp_lin_coef		==>	%1d\n"\
"chirp_quad		==>	%1d\n"\
"spare1			==>	%1x\n"\
"spare2			==>	%1x\n"\
"nought_gain		==>	%1d\n"\
"rec_gain		==>	%1d\n"\
"ant_elec_angle		==>	%1d\n"\
"ant_mech_angle		==>	%1d\n"\
"ant_elec_squint_angle		==>	%1d\n"\
"ant_mech_squint_angle		==>	%1d\n"\
"slant_first_sample		==>	%1d\n"\
"sample_delay		==>	%1d\n"\
"spare3			==>	%1x\n"\
"pltformref		==>	%1x\n"\
"ICU_onboard_time	==>	%1d\n"\
"task			==>	%1d\n"\
"image_format_counter	==>	%1d\n"\
"sample_window_start_time		==>	%1d\n"\
"pulse_repetition_interval		==>	%1d\n"\
"cal_gain		==>	%1d\n"\
"longspare			==>	%1x\n"\
"calpulse36		==>	%1x\n"

#define SARDATA_REC_RVL(SP)\
(SP)->line_no,\
(SP)->rec_index,\
(SP)->left_fill_pixels,\
(SP)->data_pixels,\
(SP)->right_fill_pixels,\
(SP)->sensor_update_flag,\
(SP)->year,\
(SP)->day_of_year,\
(SP)->msecs_of_day,\
(SP)->chan_ind,\
(SP)->chan_code,\
(SP)->trans_polar,\
(SP)->rec_polar,\
(SP)->prf,\
(SP)->spare,\
(SP)->range_flag,\
(SP)->chirp_type,\
(SP)->chirp_length,\
(SP)->chirp_const_coeff,\
(SP)->chirp_lin_coef,\
(SP)->chirp_quad,\
(SP)->spare1,\
(SP)->spare2,\
(SP)->nought_gain,\
(SP)->rec_gain,\
(SP)->ant_elec_angle,\
(SP)->ant_mech_angle,\
(SP)->ant_elec_squint_angle,\
(SP)->ant_mech_squint_angle,\
(SP)->slant_first_sample,\
(SP)->sample_delay,\
(SP)->spare3,\
(SP)->pltformref,\
(SP)->ICU_onboard_time,\
(SP)->task,\
(SP)->image_format_counter,\
(SP)->sample_window_start_time,\
(SP)->pulse_repetition_interval,\
(SP)->cal_gain,\
(SP)->longspare,\
(SP)->calpulse36

struct sardata_rec {
	int line_no;
	int rec_index;
	int left_fill_pixels;	 
	int data_pixels;
	int right_fill_pixels;
	int sensor_update_flag;
	int  year;
	int  day_of_year;
	int  msecs_of_day;
	short  chan_ind;
	short  chan_code;
	short  trans_polar;
	short  rec_polar;
	int  prf;
	int  spare;
	short range_flag;
	short chirp_type;
	int  chirp_length;
	int  chirp_const_coeff;
	int  chirp_lin_coef;
	int  chirp_quad;
	char spare1[4];
	char spare2[4];
	int  nought_gain;
	int  rec_gain;
	int  ant_elec_angle;
	int  ant_mech_angle;
	int  ant_elec_squint_angle;
	int  ant_mech_squint_angle;
  	int  slant_first_sample;
	int  sample_delay; 
	char spare3[4];
	char pltformref[64];
	int ICU_onboard_time;
	short task;
	int image_format_counter;
	short sample_window_start_time;
	short pulse_repetition_interval;
	short cal_gain;
	char longspare[120];
	char calpulse36[72];
};


/*#define SARDATA_DPAF_REC_WCS "*********** SARDATA **********\n"\
"line_no		==>	%1d\n"\
"rec_index		==>	%1d\n"\
"left_fill_pixels	==>	%1d\n"\
"data_pixels		==>	%1d\n"\
"right_fill_pixels	==>	%1d\n"\ 
"reserved_a	        ==>	%1x\n"\
"spare_a   	        ==>	%1x\n"\
"spare_b   	        ==>	%1x\n"\
"reserved_b	        ==>	%1x\n"\
"spare_c   	        ==>	%1x\n"\
"platform_reference_info                ==>     %1x\n"\
"sensor_facil_aux_info  ==>	%1x\n"\
"IDHT_general_ctr       ==>	%1x\n"\
"packet_ctr	        ==>	%1d\n"\
"subcommutation_ctr     ==>	%1d\n"\
"IDHT_gen_hdr_src_pack  ==>	%1x\n"\
"fixed_code	        ==>	%1d\n"\
"OGRC_OBRC_flag	        ==>	%1d\n"\
"ICU_onboard_time	==>	%1d\n"\
"task			==>	%1d\n"\
"image_format_counter	==>	%1d\n"\
"sample_window_start_time		==>	%1d\n"\
"pulse_repetition_interval		==>	%1d\n"\
"cal_gain		==>	%1d\n"\
"longspare			==>	%1x\n"\
"Calpulse36		==>	%1x\n"

#define SARDATA_DPAF_REC_RVL(SP)\
(SP)->line_no,\
(SP)->rec_index,\
(SP)->left_fill_pixels,\
(SP)->data_pixels,\
(SP)->right_fill_pixels,\
(SP)->reserved_a,\
(SP)->spare_a,\
(SP)->spare_b,\
(SP)->reserved_b,\
(SP)->spare_c,\
(SP)->platform_reference_info,\
(SP)->sensor_facil_aux_info,\
(SP)->IDHT_general_ctr,\
(SP)->packet_ctr,\
(SP)->subcommutation_ctr,\
(SP)->IDHT_gen_hdr_src_pack,\
(SP)->fixed_code,\
(SP)->OGRC_OBRC_flag,\
(SP)->ICU_onboard_time,\
(SP)->task,\
(SP)->image_format_counter,\
(SP)->sample_window_start_time,\
(SP)->pulse_repetition_interval,\
(SP)->cal_gain,\
(SP)->longspare,\
(SP)->calpulse36

struct sardata_dpaf_rec {
	int line_no;
	int rec_index;
	int left_fill_pixels;	 
	int data_pixels;
  	int right_fill_pixels; 
        char reserved_a[52];
	char spare_a[4];
	char spare_b[4];
	char reserved_b[32];
	char spare_c[4];
	char platform_reference_info[64];
	char sensor_facil_aux_info[220];
	char IDHT_general_ctr[10];
	int packet_ctr;
	int subcommutation_ctr;
	char IDHT_gen_hdr_src_pack[8];
	int fixed_code;
	int OGRC_OBRC_flag;
	int ICU_onboard_time;
	short task;
	int image_format_counter;
	short sample_window_start_time;
	short pulse_repetition_interval;
	short cal_gain;
	char longspare[120];
	char calpulse36[72];
};

*/

#define SARDATA_DPAF_REC_WCS "*********** SARDATA **********\n"\
"line_no		==>	%1d\n"\
"rec_index		==>	%1d\n"\
"left_fill_pixels	==>	%1d\n"\
"data_pixels		==>	%1d\n"\
"right_fill_pixels	==>	%1d\n" 

#define SARDATA_DPAF_REC_RVL(SP)\
(SP)->line_no,\
(SP)->rec_index,\
(SP)->left_fill_pixels,\
(SP)->data_pixels,\
(SP)->right_fill_pixels

struct sardata_dpaf_rec {
	int line_no;
	int rec_index;
	int left_fill_pixels;	 
	int data_pixels;
  	int right_fill_pixels; 
};

