/* adapted from Remko Scharroo's geoxyz/geocen/polcar                           */
/* This subroutine converts geodetic latitude, longitude and height above the   */
/* GRS80 reference ellipsoid to Earth Centered Fixed coordinates X, Y and Z.    */
/*                                                                              */
/* LAT     (input) : Geodetic latitude (rad).                                   */
/* LON     (input) : Geodetic longitude (rad).                                  */
/* HEIGHT  (input) : Height above the reference ellipsoid (m).                  */
/* XYZ(3) (output) : Earth Centered Fixed coordinates X, Y, and Z (m).          */
/* R      (output) : Distance to geocenter (m).                                 */

#include <stdio.h>
#include <math.h>
#include "lib_functions.h"
#include "sfd_complex.h"

double	AE,FLAT,OMFSQ,AP,GM,AM,FFACT,AMSQ;
/*------------------------------------------------------------------------------*/
void geocen(double lat, double height, double *latc, double *r)
{
double rs, lats;

lats = atan(OMFSQ*tan(lat));
rs = AP / sqrt(1.0 + FFACT * cos(lats)*cos(lats));
*r = sqrt(height*height + rs*rs + 2*height*rs*cos(lat - lats));
*latc = lats + asin(height * (sin(lat - lats))/(*r));
}
/*------------------------------------------------------------------------------*/
void polcar(double lat, double lon, double *r, double *xyz)
{
xyz[0] = *r * cos(lat) * cos(lon);
xyz[1] = *r * cos(lat) * sin(lon);
xyz[2] = *r * sin(lat);
}
/*------------------------------------------------------------------------------*/
void geoxyz(double lat, double lon, double height, double *xyz, double *r)
{
double latc;

/* from GRS80.inc */
/* set parameters for GRS80	*/

AE = 6378137.0;
FLAT = 1.0/298.257;
GM = 398600.436;
OMFSQ = (1.0 - FLAT)*(1.0 - FLAT);
AP = AE*(1.0 - FLAT);
FFACT = OMFSQ - 1.0;
AM = AE*(1.0 - FLAT/3.0 - FLAT*FLAT/5.0);
AMSQ = AM*AM;

geocen(lat, height, &latc, r);
polcar(latc, lon, r, xyz);
}
