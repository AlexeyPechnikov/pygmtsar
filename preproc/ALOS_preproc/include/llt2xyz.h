/*  @(#)physcon.h       1.1  98/12/08  */
/********1*********2*********3*********4*********5*********6*********7**
 * Name:        physcon
 * Version:     9602.17
 * Author:      M. Schenewerk
 * Purpose:     physical constants
 *
 * Global:
 * -----------
 * C_to_K            convert Celcius to Kelvin [degree K]
 * FIVE              (double)5.0
 * FOUR              (double)4.0
 * L1_frequency      GPS L1 frequency [Hz]
 * L1_wavelength     GPS L1 wavelength [m]
 * L2_frequency      GPS L2 frequency [Hz]
 * L2_wavelength     GPS L2 wavelength [m]
 * ONE               (double)1.0
 * THREE             (double)3.0
 * TWO               (double)2.0
 * ZERO              (double)0.0
 * cee               speed of light [m/s]
 * deg_to_rad        conversion for degrees to radians [rad/deg]
 * eflat             Earth flattening factor
 * emajor            Earth's semi-major axis [m]
 * eom               Earth / Moon mass ratio
 * erate             Earth's rotation rate [rad/sec]
 * pi                pi
 * pio2              pi over two
 * rad_to_deg        conversion for radians to degrees [deg/rad]
 * rad_to_hours      conversion for radians to hours [hours/rad]
 * st_to_ut          conversion for siderial to UT time [sid s/UT s]
 * soe               Sun / Earth mass ratio
 * soem              Sun / (Earth + Moon) mass ratio
 * tropical_year     tropical year [day]
 * twopi             two pi
 * ut_to_st          conversion factor for UT to siderial time [UT s/sid s]
 *
 * Notes:
 * -----------
 *
 * References:
 * -----------
 * "IERS Standards (1992)", IERS Technical Note 13, ed. D.D. McCarthy,
 * July 1992
 *
 ********1*********2*********3*********4*********5*********6*********7**
 * Modification History:
 * 9602.17, MSS, Converted to C.
 ********1*********2*********3*********4*********5*********6*********7*/

#ifndef physcon_h
#define physcon_h

#define ZERO    ((double)0.0)
#define ONE     ((double)1.0)
#define TWO     ((double)2.0)
#define THREE   ((double)3.0)
#define FOUR    ((double)4.0)
#define FIVE    ((double)5.0)

#define pi      ((double)3.14159265358979)
#define pio2    ((double)0.5*pi)
#define twopi   ((double)2.0*pi)

#define deg_to_rad      (twopi/(double)360.0)
#define rad_to_deg      ((double)360.0/twopi)
#define rad_to_hours    ((double)24.0/twopi)
#define m_to_km    ((double)1.0e-3)
#define km_to_m    ((double)1.0e+3)

#define cee     ((double)299792458.0)
#define C_to_K  ((double)273.15)

#define J2000   ((long)51545)

#define tropical_year   ((double)365.2421910)
#define ut_to_st        ((double)1.00273790934)
#define st_to_ut        ((double)0.9972695663399999)

#define emajor  ((double)6378137.0)
#define eflat   ((double)0.00335281068118 )
#define erate   ((double)7.292115855228083e-5)
#define soem    ((double)328900.550)
#define eom     ((double)81.3005870)
#define soe     (soem*((double)1.0 + (double)1.0/eom) )

#define L1_frequency    ((double)1575.420e+6)
#define L2_frequency    ((double)1227.600e+6)
#define L1_wavelength   ((double)cee/L1_frequency)
#define L2_wavelength   ((double)cee/L2_frequency)

#endif /* physcon_h */
