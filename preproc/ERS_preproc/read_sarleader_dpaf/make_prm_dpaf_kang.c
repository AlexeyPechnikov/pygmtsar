/* write a PRM file */
/* modified by Paul F. Jamason, 3/19/98 */
/* added SC_identity read from file, R. Mellors. 10/1/98 */
/* changed SC_identity to sc_identity in fprintf line, K. Watson.  3/8/99 */
/* force I_mean Q_mean = 15.5 M.Wei. 5/24/06 (Some header is bad)*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/SARtape.h"
#include "../include/soi.h"
#define FACTOR 1000000

void calc_time(char time_string[], double *clock_time); 

void make_prm_dpaf(sar)
struct	SAR_info	sar;
{
double	radar_wavelength,pulse_dur,rng_samp_rate;
double	I_mean,Q_mean,sc_clock_start,chirp_slope;

int	SC_identity; /* ERS-1 = 1, ERS2 = 2 */

double SC_clock_start,SC_clock_stop;

/* translate from structure holding SAR info */
/* see SARtape.h and associated include files for more details */

	sscanf((sar.dpaf_dss->mission_identifier+3),"%d",&SC_identity);

	prf1 = atof(sar.dpaf_dss->nominal_prf);
	radar_wavelength = atof(sar.dpaf_dss->radar_wavelength);

/* convert into seconds from MHz */

	pulse_dur = (atof(sar.dpaf_dss->range_pulse_length)/FACTOR);
	rng_samp_rate = FACTOR*(atof(sar.dpaf_dss->sampling));

/* correct by factor of two */

	chirp_slope = 2.0*atof(sar.dpaf_dss->range_pulse_phase_quad);

/* pj: add 15.5 to these values (see ERS SAR.RAW documentation) */
	
/*	I_mean = 15.5+atof(sar.dpaf_dss->dc_bias_i);
	Q_mean = 15.5+atof(sar.dpaf_dss->dc_bias_q); */
	I_mean = 15.5;
	Q_mean = 15.5; 

/* we don't want the first two character; (of year)) */

	sc_clock_start = atof((sar.dpaf_dss->satelite_clock_time+2));

/* pj */	
	calc_time(sar.dpaf_dss->zero_dop_az_time_f_pixel,&SC_clock_start);
	calc_time(sar.dpaf_dss->zero_dop_az_time_l_pixel,&SC_clock_stop);

/*fprintf(stdout,"I_mean			= %lf\n",I_mean);
fprintf(stdout,"Q_mean			= %lf\n",Q_mean);
fprintf(stdout,"rng_samp_rate		= %lg\n",rng_samp_rate);
fprintf(stdout,"chirp_slope		= %lg\n",chirp_slope);
fprintf(stdout,"pulse_dur		= %lg\n",pulse_dur);
fprintf(stdout,"radar_wavelength		= %lg\n",radar_wavelength);
fprintf(stdout, "SC_identity          = %d\n",SC_identity);*/

/* pj: write these to PRM file now instead of in read_data_file.c */

fprintf(stdout, "SC_clock_start		= %16.10lf\n",SC_clock_start);
fprintf(stdout, "SC_clock_stop          = %16.10lf\n",SC_clock_stop);

}

/* pj */

void calc_time(char time_string[],double *clock_time)

  {
    int day_of_month,hour,minute,second,msec,jday,i,month,year,leap;
    int mon_nday[]={31,28,31,30,31,30,31,31,30,31,30,31};
    char c_month[4],tmp_time_string[24];
        char *string_tmp=(char*) malloc(3);
/* need a temporary time string else msec reads beginning of next
zero_dop string */

    strncpy(tmp_time_string,&time_string[0],24);
    
/* pj: calculate SC_clock_start/stop from sarleader for DPAF */

    day_of_month = atoi(&time_string[0]);

    strncpy(c_month,&time_string[3],3);
    if (!strncmp(c_month,"JAN",3)) {
      month = 1;
    }
    if (!strncmp(c_month,"FEB",3)) {

      month = 2;
    }
    if (!strncmp(c_month,"MAR",3)) {
      month = 3;
    }
    if (!strncmp(c_month,"APR",3)) {
      month = 4;
    }
    if (!strncmp(c_month,"MAY",3)) {
      month = 5;
    }
    if (!strncmp(c_month,"JUN",3)) {
      month = 6;
    }
    if (!strncmp(c_month,"JUL",3)) {
      month = 7;
    }
    if (!strncmp(c_month,"AUG",3)) {
      month = 8;
    }
    if (!strncmp(c_month,"SEP",3)) {
      month = 9;
    }
    if (!strncmp(c_month,"OCT",3)) {
      month = 10;
    }
    if (!strncmp(c_month,"NOV",3)) {
      month = 11;
    }
    if (!strncmp(c_month,"DEC",3)) {
      month = 12;
    }
    
    year = atoi(&time_string[9]);
    
    
    if(year % 4 == 0 && year % 100 !=0){
      leap=1;
    }
    else if(year == 00){
      leap=1;
    }
    else{
      leap=0;
    }
    
    if(leap){
      mon_nday[1]=29;
    }
    
    jday=0;
    
    if(month>1){
      for(i=0;i<month-1;i++) jday=jday+mon_nday[i];
    }
    
    jday=jday+day_of_month;
    
    hour = atoi(&time_string[12]);
    
    minute = atoi(&time_string[15]);
    
    second = atoi(&time_string[18]);
    
//    msec = atoi(&tmp_time_string[21]);
//   memcpy(string_tmp,&tmp_time_string[21],3);
   msec =atoi(strncpy(string_tmp,tmp_time_string+21,3));
//    printf("%3d\n",hour);
//   printf("%3d\n",minute);
//   printf("%3d\n",second);
//   printf("%3d\n",msec);
//   printf("%s\n",tmp_time_string);
/* add either 1900 or 2000 to the year */
    if(year < 85.) {
	year = year + 2000.;
    }
    else {
	year = year + 1900.;
    }

/*     calculate clock time; determine fraction of msecs in day */

    *clock_time=(double)(year*1000+jday)+
      ((double)(hour*3600*1000+minute*60*1000+second*1000+msec)/
       (double)(24*3600*1000));
  }
