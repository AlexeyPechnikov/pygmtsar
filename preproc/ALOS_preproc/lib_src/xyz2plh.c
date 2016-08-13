/************************************************************************
* xyz2plh Converts XYZ geocentric coordinates to Phi (latitude),        *
*              Lambda (longitude), H (height) referred to an            *
*              ellipsoid of semi-major axis A and flattening FL.        *
************************************************************************/
/************************************************************************
* Creator: B. Archinal (USNO)                                           *
************************************************************************/
/************************************************************************
* Modification history:                                                 *
* 9007.20, BA,  Creation
* 9507,21, JR,  Modified for use with the page programs
* 9602.17, MSS, Converted to C.
************************************************************************/

#include <math.h>
#include "../include/llt2xyz.h"

void xyz2plh( double *xyz, double *plh, double A, double FL )
/********1*********2*********3*********4*********5*********6*********7**
 * input:
 * -----------
 * A                semi-major axis of ellipsoid [units are of distance]
 * FL               flattening of ellipsoid [unitless]
 * xyz[]            geocentric Cartesian coordinates [units are of distance]
 *
 * output:
 * -----------
 * plh[]            ellipsoidal coordinates of point, in geodetic latitude,
 *                  longitude east of Greenwich, and height [units for
 *                  latitude=plh[0] and longitude=plh[1] are in degrees;
 *                  height=plh[2] are distance and will be the same as
 *                  those of the input parameters]
 *
 * Local:
 * -----------
 * B                semi-minor axis of ellipsoid [same units as A]
 *
 * Notes:
 * -----------
 * This routine will fail for points on the Z axis, i.e. if X= Y= 0
 * (Phi = +/- 90 degrees).
 *
 * Units of input parameters `A' and `xyz' must be the same.
 *
 * References:
 * -----------
 * Borkowski, K. M. (1989).  "Accurate algorithms to transform geocentric
 * to geodetic coordinates", *Bulletin Geodesique*, v. 63, pp. 50-56.
 *
 * Borkowski, K. M. (1987).  "Transformation of geocentric to geodetic
 * coordinates without approximations", *Astrophysics and Space Science*,
 * v. 139, n. 1, pp. 1-4.  Correction in (1988), v. 146, n. 1, p. 201.
 *
 * An equivalent formulation is recommended in the IERS Standards
 * (1995), draft.
 *
 ********1*********2*********3*********4*********5*********6*********7*/
{
        double B;
        double d;
        double e;
        double f;
        double g;
        double p;
        double q;
        double r;
        double t;
        double v;
        double x= xyz[0];
        double y= xyz[1];
        double z= xyz[2];
        double zlong;
/*
 *   1.0 compute semi-minor axis and set sign to that of z in order
 *       to get sign of Phi correct
 */
        B= A * (ONE - FL);
        if( z < ZERO )
                B= -B;
/*
 *   2.0 compute intermediate values for latitude
 */
        r= sqrt( x*x + y*y );
        e= ( B*z - (A*A - B*B) ) / ( A*r );
        f= ( B*z + (A*A - B*B) ) / ( A*r );
/*
 *   3.0 find solution to:
 *       t^4 + 2*E*t^3 + 2*F*t - 1 = 0
 */
        p= (FOUR / THREE) * (e*f + ONE);
        q= TWO * (e*e - f*f);
        d= p*p*p + q*q;

        if( d >= ZERO ) {
                v= pow( (sqrt( d ) - q), (ONE / THREE) )
                 - pow( (sqrt( d ) + q), (ONE / THREE) );
        } else {
                v= TWO * sqrt( -p )
                 * cos( acos( q/(p * sqrt( -p )) ) / THREE );
        }
/*
 *   4.0 improve v
 *       NOTE: not really necessary unless point is near pole
 */
        if( v*v < fabs(p) ) {
                v= -(v*v*v + TWO*q) / (THREE*p);
        }
        g= (sqrt( e*e + v ) + e) / TWO;
        t = sqrt( g*g  + (f - v*g)/(TWO*g - e) ) - g;

        plh[0] = atan( (A*(ONE - t*t)) / (TWO*B*t) );
/*
 *   5.0 compute height above ellipsoid
 */
        plh[2]= (r - A*t)*cos( plh[0] ) + (z - B)*sin( plh[0] );
/*
 *   6.0 compute longitude east of Greenwich
 */
        zlong = atan2( y, x );
        if( zlong < ZERO )
                zlong= zlong + twopi;

        plh[1]= zlong;
/*
 *   7.0 convert latitude and longitude to degrees
 */
        plh[0] = plh[0] * rad_to_deg;
        plh[1] = plh[1] * rad_to_deg;

        return;
}
