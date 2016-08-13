#ifndef PRM_H
#define PRM_H
struct PRM {
	char input_file[256];
	char SLC_file[256];
	char out_amp_file[256];
	char out_data_file[256];
	char deskew[8];
	char iqflip[8];
	char offset_video[8];
	char srm[8];
	char ref_file[128];
	char led_file[128];
	char orbdir[8];	/* orbit direction A or D (ASCEND or DESCEND) - added by RJM*/
        char lookdir[8];/* look direction R or L (RIGHT or LEFT) */
	char dtype[8];  /* SLC data type a-SCOMPLEX integer complex, c-FCOMPLEX float complex */
	char date[16];  /* yymmdd format - skip first two digits of year - added by RJM*/

	int debug_flag;
	int bytes_per_line;
	int good_bytes;
	int first_line;
	int num_patches;
	int first_sample;
	int num_valid_az;
	int st_rng_bin;
	int num_rng_bins;
	int chirp_ext;
	int nlooks;
	int rshift;
	int ashift;
	int fdc_ystrt;
	int fdc_strt;
	int rec_start;
	int rec_stop;
	int SC_identity;	/* (1)-ERS1 (2)-ERS2 (3)-Radarsat (4)-Envisat (5)-ALOS (6)-  (7)-TSX (8)-CSK (9)-RS2 (10)-S1A */
	int ref_identity;	/* (1)-ERS1 (2)-ERS2 (3)-Radarsat (4)-Envisat (5)-ALOS (6)-  (7)-TSX (8)-CSK (9)-RS2 (10)-S1A */
	int nrows;
	int num_lines;
	int SLC_format;		/* 1 => complex ints (2 bytes)	2 => complex floats (4 bytes) */

	double SC_clock_start;	/* YYDDD.DDDD */
	double SC_clock_stop;	/* YYDDD.DDDD */
	double icu_start;	/* onboard clock counter */
	double clock_start; /* DDD.DDDDDDDD clock without year has more precision */
	double clock_stop;  /* DDD.DDDDDDDD clock without year has more precision */
	double caltone;
	double RE;			/*local earth eadius */
	double rc;			/* polar radius */
	double ra;			/* equatorial radius */
	double vel;			/* Equivalent SC velocity */
	double ht;			/* (SC_radius - RE) center */
	double ht_start;		/* (SC_radius - RE) start */
	double ht_end;			/* (SC_radius - RE) end */
	double near_range;
	double far_range;
	double prf;
	double xmi;
	double xmq;
	double az_res;
	double fs;
	double chirp_slope;
	double pulsedur;
	double lambda;
	double rhww;
	double pctbw;
	double pctbwaz;
	double fd1;
	double fdd1;
	double fddd1;
	double delr;			/* added RJM */
	double yaw;			/* added RJM 12/07*/
	double SLC_scale;

	double sub_int_r;
	double sub_int_a;
	double sub_double;
	double stretch_r;
	double stretch_a;
	double a_stretch_r;
	double a_stretch_a;
	double baseline_start;
	double baseline_center;
	double baseline_end;
	double alpha_start;
	double alpha_center;
	double alpha_end;
	double bpara;			/* parallel baseline - added by RJM */
	double bperp;			/* perpendicular baseline - added by RJM */
};
#endif /* PRM_H */
