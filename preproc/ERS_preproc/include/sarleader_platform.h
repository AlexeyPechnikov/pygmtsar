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

#define PLATFORM_RCS "%128c%4c%4c%4c%4c%4c%22c%22c%64c%22c%16c%16c%16c%48c"

#define PLATFORM_RVL(SP)\
(SP)->reserved1,\
(SP)->num_data_points,\
(SP)->year_of_data_points,\
(SP)->month_of_data_points,\
(SP)->day_of_data_points,\
(SP)->day_of_data_points_in_year,\
(SP)->sec_of_day_of_data,\
(SP)->data_points_time_gap,\
(SP)->ref_coord_sys,\
(SP)->greenwhich_mean_hour_angle,\
(SP)->a_track_pos_err,\
(SP)->c_track_pos_err,\
(SP)->radial_pos_err,\
(SP)->reserved2

struct platform {
char    reserved1[128];
char    num_data_points[4];
char    year_of_data_points[4];
char    month_of_data_points[4];
char    day_of_data_points[4];
char    day_of_data_points_in_year[4];
char    sec_of_day_of_data[22];
char    data_points_time_gap[22];
char    ref_coord_sys[64];
char    greenwhich_mean_hour_angle[22];
char    a_track_pos_err[16];
char    c_track_pos_err[16];
char    radial_pos_err[16];
char    reserved2[48];
};

#define POSITION_VECTOR_RCS "%22c%22c%22c%22c%22c%22c"

#define POSITION_VECTOR_RVL(SP)\
(SP)->pos_x,\
(SP)->pos_y,\
(SP)->pos_z,\
(SP)->vel_x,\
(SP)->vel_y,\
(SP)->vel_z

struct position_vector {
char pos_x[22] ;
char pos_y[22] ;
char pos_z[22] ;
char vel_x[22] ;
char vel_y[22] ;
char vel_z[22] ;
};
