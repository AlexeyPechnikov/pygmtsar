/************************************************************************
* plh2xyz  converts elliptic lat, lon, hgt to geocentric X, Y, Z        *
*          the library found at http://www.ngs.noaa.gov/PC_PROD/XYZWIN/ *
************************************************************************/
/************************************************************************
* Creator: Clyde Goad           (National Geodetic Survey)              *
************************************************************************/
/************************************************************************
* Modification history:                                                 *
*:8301.00,  CG, Creation                                                *
*:9406.16, MSS, Conversion to C.                                        *
*:9602.20, MSS, Stripped plh to xyz convertion from tlate.              *
*:0708.07, DTS, changed comments slightly                               *
************************************************************************/

#include <math.h>
#include "../include/llt2xyz.h"

void plh2xyz( double *plh, double *xyz, double A, double FL )
/********1*********2*********3*********4*********5*********6*********7*********
 * input parameters
 * ----------------
 * A                semi-major axis of ellipsoid [units are of distance]
 * FL               flattening of ellipsoid [unitless]
 * plh[]            ellipsoidal coordinates of point, in geodetic latitude,
 *                  longitude east of Greenwich, and height [units for
 *                  latitude=plh[0] and longitude=plh[1] are in degrees;
 *                  height=plh[2] are distance and will be the same as
 *                  those of the input parameters]
 *
 * output parameters
 * -----------------
 * xyz[]            geocentric Cartesian coordinates [units are of distance]
 *
 *
 * ------------------------------
 * Escobal, "Methods of Orbit Determination", 1965, Wiley & Sons, Inc.,
 * pp. 27-29.
 *
 * see also:
 * ------------------------------
 * xyz2plh
 */

{
        double flatfn= (TWO - FL)*FL;
        double funsq= (ONE - FL)*(ONE - FL);
        double g1;
        double g2;
        double lat_rad= deg_to_rad * plh[0];
        double lon_rad= deg_to_rad * plh[1];
        double sin_lat;


        sin_lat= sin( lat_rad );

        g1= A / sqrt( ONE - flatfn*sin_lat*sin_lat );
        g2= g1*funsq + plh[2];
        g1= g1 + plh[2];

        xyz[0]= g1 * cos( lat_rad );
        xyz[1]= xyz[0] * sin( lon_rad );
        xyz[0]= xyz[0] * cos( lon_rad );
        xyz[2]= g2 * sin_lat;
}
