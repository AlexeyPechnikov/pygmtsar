/* rewritten Dennis Milbert' solid_24d.f */
/* by
   Xiaohua Xu, 05/20/2018
   for
   GMTSAR
*/

char *USAGE = "\nsolid_tide yyyyddd.fffffff < lon_lat > lon_lat_dx_dy_dz \n\n"
              "    translated from Dennis Milbert's solid.for\n\n"
              "    compute solid earth tide given certain lon lat\n\n";

#include "lib_functions.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double pi, pi2, rad;
double a = 6378137.0;
double e2 = 6.69438002290341574957e-3;

int mjd0;

int day2date(double, double, double *, double *);
int isleapyr(int);

int compute_tide(double, double, double, double, double *, double *, double *);
int geo2xyz(double, double, double, double *, double *, double *);
int civmjd(double, double, double, double, double, double, double *, double *);
int mjdciv(double, double, double *, double *, double *, double *, double *, double *);
int sunxyz(double, double, double *);
int moonxyz(double, double, double *);
int detide(double *, double, double, double *, double *, double *);
int rge(double, double, double *, double *, double *, double, double, double);

double gps2utc(double);
double gpsleap(double);
double gps2tt(double);
int getghar(double, double, double *);
int rot3(double, double, double, double, double *);
int rot1(double, double, double, double, double *);
int sprod(double *, double *, double *, double *, double *);
double enorm8(double *);
int st1idiu(double *, double *, double *, double, double, double *);
int st1isem(double *, double *, double *, double, double, double *);
int st1l1(double *, double *, double *, double, double, double *);
int step2diu(double *, double, double, double *);
int step2lon(double *, double, double, double *);

int main(int argc, char **argv) {

	double yr, day, lon, lat, du, dv, dw;
	pi = 4.0 * atan(1.0);
	pi2 = pi * 2;
	rad = 180.0 / pi;

	if (argc != 2) {

		die("", USAGE);
	}

	yr = floor(atof(argv[1]) / 1000.0);
	day = atof(argv[1]) - yr * 1000.0;
	while (scanf(" %lf %lf ", &lon, &lat) == 2) {
		if (lon < 0.0)
			lon = lon + 360.0;
		// as quote from Dennis Milbert's website "The output file name is
		// solid.txt. It is plain ASCII text. After the header, solid earth tide
		// components are computed for 24 hours, at 1 minute intervals. Note: the
		// time stamps refer to UTC time. The solid earth tide components are NORTH,
		// EAST, UP in the local geodetic (ellipsoidal) horizon system."
		compute_tide(yr, day, lon, lat, &du, &dv, &dw);
		if (lon > 180.0)
			lon = lon - 360.0;
		fprintf(stdout, "%.9f %.9f %.12e %.12e %.12e \n", lon, lat, dv, du,
		        dw); // du = north, dv = east, dw = up
	}
	return (1);
}

int compute_tide(double yr, double day, double glod, double glad, double *du, double *dv, double *dw) {

	double gla0, glo0, eht0, x0, y0, z0, xsta[4], rsun[4], rmoon[4], etide[4];
	double fmjd, mjd, iyr, imo, idy, ihr, imn, sec;

	gla0 = glad / rad;
	glo0 = glod / rad;
	eht0 = 0.0;
	// convert yyyyddd.fffff to iyr,imo,idy,ihr,imn,sec
	day2date(yr, day, &imo, &idy);
	// printf("%lf %lf %lf %lf %lf %lf\n",yr,day,imo,idy,glad,glod);
	iyr = yr;
	ihr = floor((day - floor(day)) * 86400.0 / 3600.0);
	imn = floor(((day - floor(day)) * 86400.0 - ihr * 3600.0) / 60.0);
	sec = (day - floor(day)) * 86400.0 - ihr * 3600.0 - imn * 60.0;

	// sec = round(sec);   //for testing

	geo2xyz(gla0, glo0, eht0, &x0, &y0, &z0);
	xsta[1] = x0;
	xsta[2] = y0;
	xsta[3] = z0;

	//*** here comes the sun  (and the moon)  (go, tide!)
	// ihr=0.0;
	// imn=0.0;
	// sec=49920.0;        //GPS time system;

	civmjd(iyr, imo, idy, ihr, imn, sec, &mjd, &fmjd);
	mjd0 = mjd;
	// fprintf(stderr,"%.17lf %.17lf %lf %lf %lf %lf %lf
	// %.17lf\n",mjd,fmjd,iyr,imo,idy,ihr,imn,sec);

	// mjdciv(mjd,fmjd,&iyr,&imo,&idy,&ihr,&imn,&sec);
	// printf("%lf %lf %lf %lf %lf %lf %lf
	// %lf\n",mjd,fmjd,iyr,imo,idy,ihr,imn,sec);
	sunxyz(mjd, fmjd, rsun);
	// printf("%.17lf %.17lf %.17lf %.17lf
	// %.17lf\n",mjd,fmjd,rsun[1],rsun[2],rsun[3]);
	moonxyz(mjd, fmjd, rmoon);
	// printf("%.17lf %.17lf %.17lf %.17lf
	// %.17lf\n",mjd,fmjd,rmoon[1],rmoon[2],rmoon[3]);
	detide(xsta, mjd, fmjd, rsun, rmoon, etide);
	// printf("%.17lf %.17lf %.17lf %.17lf
	// %.17lf\n",mjd,fmjd,etide[1],etide[2],etide[3]);
	rge(gla0, glo0, du, dv, dw, etide[1], etide[2], etide[3]);
	// printf("%.17lf %.17lf %.17lf %.17lf %.17lf\n",mjd,fmjd,*du,*dv,*dw);
	// printf("hahahahaha\n");
	return (1);
}

int rge(double gla, double glo, double *u, double *v, double *w, double x, double y, double z) {

	//*** given a rectangular cartesian system (x,y,z)
	//*** compute a geodetic h cartesian sys   (u,v,w)

	// implicit double precision(a-h,o-z)

	double sb, cb, sl, cl;

	sb = sin(gla);
	cb = cos(gla);
	sl = sin(glo);
	cl = cos(glo);

	*u = -sb * cl * x - sb * sl * y + cb * z;
	*v = -sl * x + cl * y;
	*w = cb * cl * x + cb * sl * y + sb * z;

	return (1);
}

int detide(double *xsta, double mjd, double fmjd, double *xsun, double *xmon, double *dxtide) {

	//*** computation of tidal corrections of station displacements caused
	//***    by lunar and solar gravitational attraction

	//*** step 1 (here general degree 2 and 3 corrections +
	//***         call st1idiu + call st1isem + call st1l1)
	//***   + step 2 (call step2diu + call step2lon + call step2idiu)
	//*** it has been decided that the step 3 un-correction for permanent tide
	//*** would *not* be applied in order to avoid jump in the reference frame
	//*** (this step 3 must added in order to get the mean tide station position
	//*** and to be conformed with the iag resolution.)

	//*** inputs
	//***   xsta(i),i=1,2,3   -- geocentric position of the station (ITRF/ECEF)
	//***   xsun(i),i=1,2,3   -- geoc. position of the sun (ECEF)
	//***   xmon(i),i=1,2,3   -- geoc. position of the moon (ECEF)
	//***   mjd,fmjd          -- modified julian day (and fraction) (in GPS time)

	//****old calling
	// sequence*****************************************************
	//***   dmjd               -- time in mean julian date (including day
	// fraction)
	//***   fhr=hr+zmin/60.+sec/3600.   -- hr in the day

	//*** outputs
	//***   dxtide(i),i=1,2,3           -- displacement vector (ITRF)

	//*** author iers 1996 :  v. dehant, s. mathews and j. gipson
	//***    (test between two subroutines)
	//*** author iers 2000 :  v. dehant, c. bruyninx and s. mathews
	//***    (test in the bernese program by c. bruyninx)

	//*** created:  96/03/23 (see above)
	//*** modified from dehanttideinelMJD.f by Dennis Milbert 2006sep10
	//*** bug fix regarding fhr (changed calling sequence, too)
	//*** modified to reflect table 7.5a and b IERS Conventions 2003
	//*** modified to use TT time system to call step 2 functions
	//*** sign correction by V.Dehant to match eq.16b, p.81, Conventions
	//*** applied by Dennis Milbert 2007may05

	// implicit double precision(a-h,o-z)
	// double precision xsta(3),xsun(3),xmon(3),dxtide(3),xcorsta(3)
	// double precision h20,l20,h3,l3,h2,l2
	// double precision mass_ratio_sun,mass_ratio_moon

	double h20, l20, h3, l3, h2, l2;
	double rsta, rsun, rmon, scs, scm, scsun, scmon, cosphi;
	double tsecgps, tsectt, fmjdtt, dmjdtt, t, fhr;
	double p2sun, p2mon, p3sun, p3mon, x2sun, x2mon, x3sun, x3mon;
	double mass_ratio_sun, mass_ratio_moon, re, fac2sun, fac2mon, fac3sun, fac3mon, xcorsta[4];
	int i;

	//*** nominal second degree and third degree love numbers and shida numbers

	// data h20/0.6078d0/,l20/0.0847d0/,h3/0.292d0/,l3/0.015d0/
	h20 = 0.6078;
	l20 = 0.0847;
	h3 = 0.292;
	l3 = 0.015;

	//*** internal support for new calling sequence
	//*** also convert GPS time into TT time

	tsecgps = fmjd * 86400.0;  //!*** GPS time (sec of day)
	tsectt = gps2tt(tsecgps);  //!*** TT  time (sec of day)
	fmjdtt = tsectt / 86400.0; //!*** TT  time (fract. day)

	dmjdtt = mjd + fmjdtt; //!*** float MJD in TT
	//*** commented line was live code in dehanttideinelMJD.f
	//*** changed on the suggestion of Dr. Don Kim, UNB -- 09mar21
	//*** Julian date for 2000 January 1 00:00:00.0 UT is  JD 2451544.5
	//*** MJD         for 2000 January 1 00:00:00.0 UT is MJD   51544.0
	//***** t=(dmjdtt-51545.d0)/36525.d0                !*** days to centuries, TT
	t = (dmjdtt - 51544.0) / 36525.0;      //!*** days to centuries, TT
	fhr = (dmjdtt - floor(dmjdtt)) * 24.0; //!*** hours in the day, TT

	//*** scalar product of station vector with sun/moon vector

	sprod(xsta, xsun, &scs, &rsta, &rsun);
	sprod(xsta, xmon, &scm, &rsta, &rmon);
	scsun = scs / rsta / rsun;
	scmon = scm / rsta / rmon;

	//*** computation of new h2 and l2

	cosphi = sqrt(xsta[1] * xsta[1] + xsta[2] * xsta[2]) / rsta;
	h2 = h20 - 0.0006 * (1.0 - 3.0 / 2.0 * cosphi * cosphi);
	l2 = l20 + 0.0002 * (1.0 - 3.0 / 2.0 * cosphi * cosphi);

	//*** p2-term

	p2sun = 3.0 * (h2 / 2.0 - l2) * scsun * scsun - h2 / 2.0;
	p2mon = 3.0 * (h2 / 2.0 - l2) * scmon * scmon - h2 / 2.0;

	//*** p3-term

	p3sun = 5.0 / 2.0 * (h3 - 3.0 * l3) * scsun * scsun * scsun + 3.0 / 2.0 * (l3 - h3) * scsun;
	p3mon = 5.0 / 2.0 * (h3 - 3.0 * l3) * scmon * scmon * scmon + 3.0 / 2.0 * (l3 - h3) * scmon;

	//*** term in direction of sun/moon vector

	x2sun = 3.0 * l2 * scsun;
	x2mon = 3.0 * l2 * scmon;
	x3sun = 3.0 * l3 / 2.0 * (5.0 * scsun * scsun - 1.0);
	x3mon = 3.0 * l3 / 2.0 * (5.0 * scmon * scmon - 1.0);

	//*** factors for sun/moon

	mass_ratio_sun = 332945.943062;
	mass_ratio_moon = 0.012300034;
	re = 6378136.55;
	fac2sun = mass_ratio_sun * re * pow(re / rsun, 3);
	fac2mon = mass_ratio_moon * re * pow(re / rmon, 3);
	fac3sun = fac2sun * (re / rsun);
	fac3mon = fac2mon * (re / rmon);

	//*** total displacement

	for (i = 1; i <= 3; i++) {
		dxtide[i] = fac2sun * (x2sun * xsun[i] / rsun + p2sun * xsta[i] / rsta) +
		            fac2mon * (x2mon * xmon[i] / rmon + p2mon * xsta[i] / rsta) +
		            fac3sun * (x3sun * xsun[i] / rsun + p3sun * xsta[i] / rsta) +
		            fac3mon * (x3mon * xmon[i] / rmon + p3mon * xsta[i] / rsta);
	}
	// call zero_vec8(xcorsta)
	xcorsta[1] = 0;
	xcorsta[2] = 0;
	xcorsta[3] = 0;
	xcorsta[0] = 0;
	// printf("%.15lf %.15lf %.15lf %.15lf %.15lf %.15lf
	// %.15lf\n",rsta,xsta[1],xsta[2],xsta[3],dxtide[1],dxtide[2],dxtide[3]);
	//*** corrections for the out-of-phase part of love numbers
	//***     (part h_2^(0)i and l_2^(0)i )

	//*** first, for the diurnal band

	st1idiu(xsta, xsun, xmon, fac2sun, fac2mon, xcorsta);
	dxtide[1] = dxtide[1] + xcorsta[1];
	dxtide[2] = dxtide[2] + xcorsta[2];
	dxtide[3] = dxtide[3] + xcorsta[3];
	// printf("%.15lf %.15lf %.15lf\n",xcorsta[1],xcorsta[2],xcorsta[3]);

	//*** second, for the semi-diurnal band

	st1isem(xsta, xsun, xmon, fac2sun, fac2mon, xcorsta);
	dxtide[1] = dxtide[1] + xcorsta[1];
	dxtide[2] = dxtide[2] + xcorsta[2];
	dxtide[3] = dxtide[3] + xcorsta[3];
	// printf("%.15lf %.15lf %.15lf\n",xcorsta[1],xcorsta[2],xcorsta[3]);
	//*** corrections for the latitude dependence of love numbers (part l^(1) )

	st1l1(xsta, xsun, xmon, fac2sun, fac2mon, xcorsta);
	dxtide[1] = dxtide[1] + xcorsta[1];
	dxtide[2] = dxtide[2] + xcorsta[2];
	dxtide[3] = dxtide[3] + xcorsta[3];
	// printf("%.15lf %.15lf %.15lf\n",xcorsta[1],xcorsta[2],xcorsta[3]);

	//*** consider corrections for step 2
	//*** corrections for the diurnal band:

	//***  first, we need to know the date converted in julian centuries

	//***  this is now handled at top of code   (also convert to TT time system)
	//***** t=(dmjd-51545.)/36525.
	//***** fhr=dmjd-int(dmjd)             !*** this is/was a buggy line (day vs.
	// hr)

	//***  second, the diurnal band corrections,
	//***   (in-phase and out-of-phase frequency dependence):

	step2diu(xsta, fhr, t, xcorsta);
	dxtide[1] = dxtide[1] + xcorsta[1];
	dxtide[2] = dxtide[2] + xcorsta[2];
	dxtide[3] = dxtide[3] + xcorsta[3];

	//***  corrections for the long-period band,
	//***   (in-phase and out-of-phase frequency dependence):

	step2lon(xsta, fhr, t, xcorsta);
	dxtide[1] = dxtide[1] + xcorsta[1];
	dxtide[2] = dxtide[2] + xcorsta[2];
	dxtide[3] = dxtide[3] + xcorsta[3];

	/**** consider corrections for step 3
	 *-----------------------------------------------------------------------
	 * The code below is commented to prevent restoring deformation
	 * due to permanent tide.  All the code above removes
	 * total tidal deformation with conventional Love numbers.
	 * The code above realizes a conventional tide free crust (i.e. ITRF).
	 * This does NOT conform to Resolution 16 of the 18th General Assembly
	 * of the IAG (1983).  This resolution has not been implemented by
	 * the space geodesy community in general (c.f. IERS Conventions 2003).
	 *-----------------------------------------------------------------------*/

	//*** uncorrect for the permanent tide  (only if you want mean tide system)

	//***   pi=3.141592654
	//***   sinphi=xsta(3)/rsta
	//***   cosphi=dsqrt(xsta(1)**2+xsta(2)**2)/rsta
	//***   cosla=xsta(1)/cosphi/rsta
	//***   sinla=xsta(2)/cosphi/rsta
	//***   dr=-dsqrt(5./4./pi)*h2*0.31460*(3./2.*sinphi**2-0.5)
	//***   dn=-dsqrt(5./4./pi)*l2*0.31460*3.*cosphi*sinphi
	//***   dxtide(1)=dxtide(1)-dr*cosla*cosphi+dn*cosla*sinphi
	//***   dxtide(2)=dxtide(2)-dr*sinla*cosphi+dn*sinla*sinphi
	//***   dxtide(3)=dxtide(3)-dr*sinphi      -dn*cosphi

	return (1);
}

int step2lon(double *xsta, double fhr, double t, double *xcorsta) {

	// implicit double precision (a-h,o-z)
	// double precision deg2rad
	// double precision xsta(3),xcorsta(3),datdi(9,5)
	double deg2rad = 0.017453292519943295769;
	double s, pr, h, p, zns, ps, rsta, sinphi, cosphi, cosla, sinla, dr_tot, dn_tot, thetaf, dr, dn, de;
	int i, j;
	//*** cf. table 7.5b of IERS conventions 2003 (TN.32, pg.82)
	//*** columns are s,h,p,N',ps, dR(ip),dT(ip),dR(op),dT(op)
	//*** IERS cols.= s,h,p,N',ps, dR(ip),dR(op),dT(ip),dT(op)
	//*** units of mm
	double datdi[45] /* was [9][5] */ = {0.,   0.,   0.,   1.,   0.,   .47,  .23, .16, .07,  0.,   2.,   0.,   0.,   0.,   -.2,
	                                     -.12, -.11, -.05, 1.,   0.,   -1.,  0.,  0.,  -.11, -.08, -.09, -.04, 2.,   0.,   0.,
	                                     0.,   0.,   -.13, -.11, -.15, -.07, 2.,  0.,  0.,   1.,   0.,   -.05, -.05, -.06, -.03};
	// data ((datdi(i,j),i=1,9),j=1,5)/
	//*   0, 0, 0, 1, 0,   0.47, 0.23, 0.16, 0.07,
	//*   0, 2, 0, 0, 0,  -0.20,-0.12,-0.11,-0.05,
	//*   1, 0,-1, 0, 0,  -0.11,-0.08,-0.09,-0.04,
	//*   2, 0, 0, 0, 0,  -0.13,-0.11,-0.15,-0.07,
	//*   2, 0, 0, 1, 0,  -0.05,-0.05,-0.06,-0.03/

	s = 218.31664563 + 481267.88194 * t - 0.0014663889 * t * t + 0.00000185139 * t * t * t;
	pr = 1.396971278 * t + 0.000308889 * t * t + 0.000000021 * t * t * t + 0.000000007 * t * t * t * t;
	s = s + pr;
	h = 280.46645 + 36000.7697489 * t + 0.00030322222 * t * t + 0.000000020 * t * t * t - 0.00000000654 * t * t * t * t;
	p = 83.35324312 + 4069.01363525 * t - 0.01032172222 * t * t - 0.0000124991 * t * t * t + 0.00000005263 * t * t * t * t;
	zns = 234.95544499 + 1934.13626197 * t - 0.00207561111 * t * t - 0.00000213944 * t * t * t + 0.00000001650 * t * t * t * t;
	ps = 282.93734098 + 1.71945766667 * t + 0.00045688889 * t * t - 0.00000001778 * t * t * t - 0.00000000334 * t * t * t * t;

	rsta = sqrt(xsta[1] * xsta[1] + xsta[2] * xsta[2] + xsta[3] * xsta[3]);
	sinphi = xsta[3] / rsta;
	cosphi = sqrt(xsta[1] * xsta[1] + xsta[2] * xsta[2]) / rsta;
	cosla = xsta[1] / cosphi / rsta;
	sinla = xsta[2] / cosphi / rsta;

	//*** reduce angles to between 0 and 360

	s = fmod(s, 360.0);
	//***** tau=dmod(tau,360.d0)       !*** tau not used here--09jul28
	h = fmod(h, 360.0);
	p = fmod(p, 360.0);
	zns = fmod(zns, 360.0);
	ps = fmod(ps, 360.0);

	dr_tot = 0.0;
	dn_tot = 0.0;
	for (i = 1; i <= 3; i++) {
		xcorsta[i] = 0.0;
	}

	//***             1 2 3 4   5   6      7      8      9
	//*** columns are s,h,p,N',ps, dR(ip),dT(ip),dR(op),dT(op)

	for (j = 1; j <= 5; j++) {
		thetaf = (datdi[j * 9 - 9] * s + datdi[j * 9 - 8] * h + datdi[j * 9 - 7] * p + datdi[j * 9 - 6] * zns +
		          datdi[j * 9 - 5] * ps) *
		         deg2rad;
		dr = datdi[j * 9 - 4] * (3.0 * sinphi * sinphi - 1.0) / 2.0 * cos(thetaf) +
		     datdi[j * 9 - 2] * (3.0 * sinphi * sinphi - 1.0) / 2.0 * sin(thetaf);
		dn = datdi[j * 9 - 3] * (cosphi * sinphi * 2.0) * cos(thetaf) + datdi[j * 9 - 1] * (cosphi * sinphi * 2.0) * sin(thetaf);
		de = 0.0;
		dr_tot = dr_tot + dr;
		dn_tot = dn_tot + dn;
		xcorsta[1] = xcorsta[1] + dr * cosla * cosphi - de * sinla - dn * sinphi * cosla;
		xcorsta[2] = xcorsta[2] + dr * sinla * cosphi + de * cosla - dn * sinphi * sinla;
		xcorsta[3] = xcorsta[3] + dr * sinphi + dn * cosphi;
	}

	for (i = 1; i <= 3; i++) {
		xcorsta[i] = xcorsta[i] / 1000.0;
	}

	return (1);
}

int step2diu(double *xsta, double fhr, double t, double *xcorsta) {

	//*** last change:  vd   17 may 00   1:20 pm
	//*** these are the subroutines for the step2 of the tidal corrections.
	//*** they are called to account for the frequency dependence
	//*** of the love numbers.

	// implicit double precision (a-h,o-z)
	// double precision xsta(3),xcorsta(3),datdi(9,31)
	// double precision deg2rad
	// data deg2rad/0.017453292519943295769d0/
	double deg2rad = 0.017453292519943295769;
	double s, tau, pr, h, zns, p, ps, rsta, sinphi, cosphi, cosla, sinla, zla, thetaf, dr, dn, de;
	int i, j;
	//*** note, following table is derived from dehanttideinelMJD.f (2000oct30
	// 16:10)
	//*** has minor differences from that of dehanttideinel.f (2000apr17 14:10)
	//*** D.M. edited to strictly follow published table 7.5a (2006aug08 13:46)

	//*** cf. table 7.5a of IERS conventions 2003 (TN.32, pg.82)
	//*** columns are s,h,p,N',ps, dR(ip),dR(op),dT(ip),dT(op)
	//*** units of mm

	double datdi[279] /* was [9][31] */ = {
	    -3.,   0.,   2.,   0.,   0.,   -.01, -.01, 0.,   0.,   -3., 2.,   0.,  0.,   0.,  -.01, -.01, 0.,  0.,   -2.,  0.,
	    1.,    -1.,  0.,   -.02, -.01, 0.,   0.,   -2.,  0.,   1.,  0.,   0.,  -.08, 0.,  .01,  .01,  -2., 2.,   -1.,  0.,
	    0.,    -.02, -.01, 0.,   0.,   -1.,  0.,   0.,   -1.,  0.,  -.1,  0.,  0.,   0.,  -1.,  0.,   0.,  0.,   0.,   -.51,
	    0.,    -.02, .03,  -1.,  2.,   0.,   0.,   0.,   .01,  0.,  0.,   0.,  0.,   -2., 1.,   0.,   0.,  .01,  0.,   0.,
	    0.,    0.,   0.,   -1.,  0.,   0.,   .02,  .01,  0.,   0.,  0.,   0.,  1.,   0.,  0.,   .06,  0.,  0.,   0.,   0.,
	    0.,    1.,   1.,   0.,   .01,  0.,   0.,   0.,   0.,   2.,  -1.,  0.,  0.,   .01, 0.,   0.,   0.,  1.,   -3.,  0.,
	    0.,    1.,   -.06, 0.,   0.,   0.,   1.,   -2.,  0.,   1.,  0.,   .01, 0.,   0.,  0.,   1.,   -2., 0.,   0.,   0.,
	    -1.23, -.07, .06,  .01,  1.,   -1.,  0.,   0.,   -1.,  .02, 0.,   0.,  0.,   1.,  -1.,  0.,   0.,  1.,   .04,  0.,
	    0.,    0.,   1.,   0.,   0.,   -1.,  0.,   -.22, .01,  .01, 0.,   1.,  0.,   0.,  0.,   0.,   12., -.78, -.67, -.03,
	    1.,    0.,   0.,   1.,   0.,   1.73, -.12, -.1,  0.,   1.,  0.,   0.,  2.,   0.,  -.04, 0.,   0.,  0.,   1.,   1.,
	    0.,    0.,   -1.,  -.5,  -.01, .03,  0.,   1.,   1.,   0.,  0.,   1.,  .01,  0.,  0.,   0.,   1.,  1.,   0.,   1.,
	    -1.,   -.01, 0.,   0.,   0.,   1.,   2.,   -2.,  0.,   0.,  -.01, 0.,  0.,   0.,  1.,   2.,   0.,  0.,   0.,   -.11,
	    .01,   .01,  0.,   2.,   -2.,  1.,   0.,   0.,   -.01, 0.,  0.,   0.,  2.,   0.,  -1.,  0.,   0.,  -.02, .02,  0.,
	    .01,   3.,   0.,   0.,   0.,   0.,   0.,   .01,  0.,   .01, 3.,   0.,  0.,   1.,  0.,   0.,   .01, 0.,   0.};
	/*
	data ((datdi(i,j),i=1,9),j=1,31)/
	* -3., 0., 2., 0., 0.,-0.01,-0.01, 0.0 , 0.0,
	* -3., 2., 0., 0., 0.,-0.01,-0.01, 0.0 , 0.0,
	* -2., 0., 1.,-1., 0.,-0.02,-0.01, 0.0 , 0.0,
	*****-----------------------------------------------------------------------
	****** -2., 0., 1., 0., 0.,-0.08,-0.05, 0.01,-0.02,      !*** original entry
	* -2., 0., 1., 0., 0.,-0.08, 0.00, 0.01, 0.01,      !*** table 7.5a
	*****-----------------------------------------------------------------------
	* -2., 2.,-1., 0., 0.,-0.02,-0.01, 0.0 , 0.0,
	*****-----------------------------------------------------------------------
	****** -1., 0., 0.,-1., 0.,-0.10,-0.05, 0.0 ,-0.02,      !*** original entry
	* -1., 0., 0.,-1., 0.,-0.10, 0.00, 0.00, 0.00,      !*** table 7.5a
	*****-----------------------------------------------------------------------
	****** -1., 0., 0., 0., 0.,-0.51,-0.26,-0.02,-0.12,      !*** original entry
	* -1., 0., 0., 0., 0.,-0.51, 0.00,-0.02, 0.03,      !*** table 7.5a
	*****-----------------------------------------------------------------------
	* -1., 2., 0., 0., 0., 0.01, 0.0 , 0.0 , 0.0,
	*  0.,-2., 1., 0., 0., 0.01, 0.0 , 0.0 , 0.0,
	*  0., 0.,-1., 0., 0., 0.02, 0.01, 0.0 , 0.0,
	*****-----------------------------------------------------------------------
	******  0., 0., 1., 0., 0., 0.06, 0.02, 0.0 , 0.01,      !*** original entry
	*  0., 0., 1., 0., 0., 0.06, 0.00, 0.00, 0.00,      !*** table 7.5a
	*****-----------------------------------------------------------------------
	*  0., 0., 1., 1., 0., 0.01, 0.0 , 0.0 , 0.0,
	*  0., 2.,-1., 0., 0., 0.01, 0.0 , 0.0 , 0.0,
	*  1.,-3., 0., 0., 1.,-0.06, 0.00, 0.00, 0.00,      !*** table 7.5a
	*  1.,-2., 0., 1., 0., 0.01, 0.0 , 0.0 , 0.0,
	*****-----------------------------------------------------------------------
	******  1.,-2., 0., 0., 0.,-1.23,-0.05, 0.06,-0.06,      !*** original entry
	*  1.,-2., 0., 0., 0.,-1.23,-0.07, 0.06, 0.01,      !*** table 7.5a
	*****-----------------------------------------------------------------------
	*  1.,-1., 0., 0.,-1., 0.02, 0.0 , 0.0 , 0.0,
	*  1.,-1., 0., 0., 1., 0.04, 0.0 , 0.0 , 0.0,
	*  1., 0., 0.,-1., 0.,-0.22, 0.01, 0.01, 0.00,      !*** table 7.5a
	*****-----------------------------------------------------------------------
	******  1., 0., 0., 0., 0.,12.02,-0.45,-0.66, 0.17,      !*** original entry
	*  1., 0., 0., 0., 0.,12.00,-0.78,-0.67,-0.03,      !*** table 7.5a
	*****-----------------------------------------------------------------------
	******  1., 0., 0., 1., 0., 1.73,-0.07,-0.10, 0.02,      !*** original entry
	*  1., 0., 0., 1., 0., 1.73,-0.12,-0.10, 0.00,      !*** table 7.5a
	*****-----------------------------------------------------------------------
	*  1., 0., 0., 2., 0.,-0.04, 0.0 , 0.0 , 0.0,
	*****-----------------------------------------------------------------------
	******  1., 1., 0., 0.,-1.,-0.50, 0.0 , 0.03, 0.0,       !*** original entry
	*  1., 1., 0., 0.,-1.,-0.50,-0.01, 0.03, 0.00,      !*** table 7.5a
	*****-----------------------------------------------------------------------
	*  1., 1., 0., 0., 1., 0.01, 0.0 , 0.0 , 0.0,
	*****-----------------------------------------------------------------------
	******  0., 1., 0., 1.,-1.,-0.01, 0.0 , 0.0 , 0.0,       !*** original entry
	*  1., 1., 0., 1.,-1.,-0.01, 0.0 , 0.0 , 0.0,       !*** v.dehant 2007
	*****-----------------------------------------------------------------------
	*  1., 2.,-2., 0., 0.,-0.01, 0.0 , 0.0 , 0.0,
	*****-----------------------------------------------------------------------
	******  1., 2., 0., 0., 0.,-0.12, 0.01, 0.01, 0.0,       !*** original entry
	*  1., 2., 0., 0., 0.,-0.11, 0.01, 0.01, 0.00,      !*** table 7.5a
	*****-----------------------------------------------------------------------
	*  2.,-2., 1., 0., 0.,-0.01, 0.0 , 0.0 , 0.0,
	*  2., 0.,-1., 0., 0.,-0.02, 0.02, 0.0 , 0.01,
	*  3., 0., 0., 0., 0., 0.0 , 0.01, 0.0 , 0.01,
	*  3., 0., 0., 1., 0., 0.0 , 0.01, 0.0 , 0.0/
	*/

	s = 218.31664563 + 481267.88194 * t - 0.0014663889 * t * t + 0.00000185139 * t * t * t;
	tau = fhr * 15.0 + 280.46061840 + 36000.77005360 * t + 0.000387930 * t * t - 0.0000000258 * t * t * t - s;
	pr = 1.396971278 * t + 0.000308889 * t * t + 0.000000021 * t * t * t + 0.000000007 * t * t * t * t;
	s = s + pr;
	h = 280.46645 + 36000.7697489 * t + 0.00030322222 * t * t + 0.000000020 * t * t * t - 0.00000000654 * t * t * t * t;
	p = 83.35324312 + 4069.01363525 * t - 0.01032172222 * t * t - 0.0000124991 * t * t * t + 0.00000005263 * t * t * t * t;
	zns = 234.95544499 + 1934.13626197 * t - 0.00207561111 * t * t - 0.00000213944 * t * t * t + 0.00000001650 * t * t * t * t;
	ps = 282.93734098 + 1.71945766667 * t + 0.00045688889 * t * t - 0.00000001778 * t * t * t - 0.00000000334 * t * t * t * t;

	//*** reduce angles to between 0 and 360

	s = fmod(s, 360.0);
	tau = fmod(tau, 360.0);
	h = fmod(h, 360.0);
	p = fmod(p, 360.0);
	zns = fmod(zns, 360.0);
	ps = fmod(ps, 360.0);

	rsta = sqrt(xsta[1] * xsta[1] + xsta[2] * xsta[2] + xsta[3] * xsta[3]);
	sinphi = xsta[3] / rsta;
	cosphi = sqrt(xsta[1] * xsta[1] + xsta[2] * xsta[2]) / rsta;

	cosla = xsta[1] / cosphi / rsta;
	sinla = xsta[2] / cosphi / rsta;
	zla = atan2(xsta[2], xsta[1]);
	for (i = 1; i <= 3; i++) {
		xcorsta[i] = 0.0;
	}
	for (j = 1; j <= 31; j++) {
		thetaf = (tau + datdi[j * 9 - 9] * s + datdi[j * 9 - 8] * h + datdi[j * 9 - 7] * p + datdi[j * 9 - 6] * zns +
		          datdi[j * 9 - 5] * ps) *
		         deg2rad;
		dr = datdi[j * 9 - 4] * 2.0 * sinphi * cosphi * sin(thetaf + zla) +
		     datdi[j * 9 - 3] * 2.0 * sinphi * cosphi * cos(thetaf + zla);
		dn = datdi[j * 9 - 2] * (cosphi * cosphi - sinphi * sinphi) * sin(thetaf + zla) +
		     datdi[j * 9 - 1] * (cosphi * cosphi - sinphi * sinphi) * cos(thetaf + zla);
		//***** following correction by V.Dehant to match eq.16b, p.81, 2003
		// Conventions
		//*****   de=datdi(8,j)*sinphi*cos(thetaf+zla)+
		de = datdi[j * 9 - 2] * sinphi * cos(thetaf + zla) - datdi[j * 9 - 1] * sinphi * sin(thetaf + zla);
		xcorsta[1] = xcorsta[1] + dr * cosla * cosphi - de * sinla - dn * sinphi * cosla;
		xcorsta[2] = xcorsta[2] + dr * sinla * cosphi + de * cosla - dn * sinphi * sinla;
		xcorsta[3] = xcorsta[3] + dr * sinphi + dn * cosphi;
	}

	for (i = 1; i <= 3; i++) {
		xcorsta[i] = xcorsta[i] / 1000.0;
	}

	return (1);
}

int st1l1(double *xsta, double *xsun, double *xmon, double fac2sun, double fac2mon, double *xcorsta) {

	//*** this subroutine gives the corrections induced by the latitude dependence
	//*** given by l^(1) in mahtews et al (1991)

	//***  input: xsta,xsun,xmon,fac3sun,fac3mon
	//*** output: xcorsta

	// implicit double precision (a-h,o-z)
	// dimension xsta(3),xsun(3),xmon(3),xcorsta(3)
	double l1, l1d, l1sd, rsta, sinphi, cosphi, sinla, cosla, costwola, sintwola, rmon, rsun, dnsun, dnmon, desun, demon, dn, de;
	;
	l1d = 0.0012;
	l1sd = 0.0024;

	rsta = enorm8(xsta);
	sinphi = xsta[3] / rsta;
	cosphi = sqrt(xsta[1] * xsta[1] + xsta[2] * xsta[2]) / rsta;
	sinla = xsta[2] / cosphi / rsta;
	cosla = xsta[1] / cosphi / rsta;
	rmon = enorm8(xmon);
	rsun = enorm8(xsun);

	//*** for the diurnal band

	l1 = l1d;
	dnsun = -l1 * sinphi * sinphi * fac2sun * xsun[3] * (xsun[1] * cosla + xsun[2] * sinla) / rsun / rsun;
	dnmon = -l1 * sinphi * sinphi * fac2mon * xmon[3] * (xmon[1] * cosla + xmon[2] * sinla) / rmon / rmon;
	desun =
	    l1 * sinphi * (cosphi * cosphi - sinphi * sinphi) * fac2sun * xsun[3] * (xsun[1] * sinla - xsun[2] * cosla) / rsun / rsun;
	demon =
	    l1 * sinphi * (cosphi * cosphi - sinphi * sinphi) * fac2mon * xmon[3] * (xmon[1] * sinla - xmon[2] * cosla) / rmon / rmon;
	de = 3.0 * (desun + demon);
	dn = 3.0 * (dnsun + dnmon);
	xcorsta[1] = -de * sinla - dn * sinphi * cosla;
	xcorsta[2] = de * cosla - dn * sinphi * sinla;
	xcorsta[3] = dn * cosphi;

	//*** for the semi-diurnal band

	l1 = l1sd;
	costwola = cosla * cosla - sinla * sinla;
	sintwola = 2.0 * cosla * sinla;
	dnsun = -l1 / 2.0 * sinphi * cosphi * fac2sun *
	        ((xsun[1] * xsun[1] - xsun[2] * xsun[2]) * costwola + 2.0 * xsun[1] * xsun[2] * sintwola) / rsun / rsun;
	dnmon = -l1 / 2.0 * sinphi * cosphi * fac2mon *
	        ((xmon[1] * xmon[1] - xmon[2] * xmon[2]) * costwola + 2.0 * xmon[1] * xmon[2] * sintwola) / rmon / rmon;
	desun = -l1 / 2.0 * sinphi * sinphi * cosphi * fac2sun *
	        ((xsun[1] * xsun[1] - xsun[2] * xsun[2]) * sintwola - 2.0 * xsun[1] * xsun[2] * costwola) / rsun / rsun;
	demon = -l1 / 2.0 * sinphi * sinphi * cosphi * fac2mon *
	        ((xmon[1] * xmon[1] - xmon[2] * xmon[2]) * sintwola - 2.0 * xmon[1] * xmon[2] * costwola) / rmon / rmon;
	de = 3.0 * (desun + demon);
	dn = 3.0 * (dnsun + dnmon);
	xcorsta[1] = xcorsta[1] - de * sinla - dn * sinphi * cosla;
	xcorsta[2] = xcorsta[2] + de * cosla - dn * sinphi * sinla;
	xcorsta[3] = xcorsta[3] + dn * cosphi;

	return (1);
}

int st1isem(double *xsta, double *xsun, double *xmon, double fac2sun, double fac2mon, double *xcorsta) {

	//*** this subroutine gives the out-of-phase corrections induced by
	//*** mantle inelasticity in the diurnal band

	//***  input: xsta,xsun,xmon,fac2sun,fac2mon
	//*** output: xcorsta

	// implicit double precision (a-h,o-z)
	// dimension xsta(3),xsun(3),xmon(3),xcorsta(3)
	// data dhi/-0.0022d0/,dli/-0.0007d0/

	double dhi, dli, rsta, sinphi, cosphi, sinla, cosla, costwola, sintwola, rmon, rsun, drsun, drmon, dnsun, dnmon, desun, demon,
	    dr, dn, de;

	dhi = -0.0022;
	dli = -0.0007;

	rsta = enorm8(xsta);
	sinphi = xsta[3] / rsta;
	cosphi = sqrt(xsta[1] * xsta[1] + xsta[2] * xsta[2]) / rsta;
	sinla = xsta[2] / cosphi / rsta;
	cosla = xsta[1] / cosphi / rsta;
	costwola = cosla * cosla - sinla * sinla;
	sintwola = 2.0 * cosla * sinla;
	rmon = enorm8(xmon);
	rsun = enorm8(xsun);
	drsun = -3.0 / 4.0 * dhi * cosphi * cosphi * fac2sun *
	        ((xsun[1] * xsun[1] - xsun[2] * xsun[2]) * sintwola - 2.0 * xsun[1] * xsun[2] * costwola) / rsun / rsun;
	drmon = -3.0 / 4.0 * dhi * cosphi * cosphi * fac2mon *
	        ((xmon[1] * xmon[1] - xmon[2] * xmon[2]) * sintwola - 2.0 * xmon[1] * xmon[2] * costwola) / rmon / rmon;
	dnsun = 1.50 * dli * sinphi * cosphi * fac2sun *
	        ((xsun[1] * xsun[1] - xsun[2] * xsun[2]) * sintwola - 2.0 * xsun[1] * xsun[2] * costwola) / rsun / rsun;
	dnmon = 1.50 * dli * sinphi * cosphi * fac2mon *
	        ((xmon[1] * xmon[1] - xmon[2] * xmon[2]) * sintwola - 2.0 * xmon[1] * xmon[2] * costwola) / rmon / rmon;
	desun = -3.0 / 2.0 * dli * cosphi * fac2sun *
	        ((xsun[1] * xsun[1] - xsun[2] * xsun[2]) * costwola + 2.0 * xsun[1] * xsun[2] * sintwola) / rsun / rsun;
	demon = -3.0 / 2.0 * dli * cosphi * fac2mon *
	        ((xmon[1] * xmon[1] - xmon[2] * xmon[2]) * costwola + 2.0 * xmon[1] * xmon[2] * sintwola) / rmon / rmon;
	dr = drsun + drmon;
	dn = dnsun + dnmon;
	de = desun + demon;
	xcorsta[1] = dr * cosla * cosphi - de * sinla - dn * sinphi * cosla;
	xcorsta[2] = dr * sinla * cosphi + de * cosla - dn * sinphi * sinla;
	xcorsta[3] = dr * sinphi + dn * cosphi;

	return (1);
}

int st1idiu(double *xsta, double *xsun, double *xmon, double fac2sun, double fac2mon, double *xcorsta) {

	//*** this subroutine gives the out-of-phase corrections induced by
	//*** mantle inelasticity in the diurnal band

	//***  input: xsta,xsun,xmon,fac2sun,fac2mon
	//*** output: xcorsta

	// implicit double precision (a-h,o-z)
	// dimension xsta(3),xsun(3),xmon(3),xcorsta(3)
	// data dhi/-0.0025d0/,dli/-0.0007d0/
	double dhi, dli, rsta, sinphi, cosphi, cos2phi, sinla, cosla, rmon, rsun, drsun, drmon, dnsun, dnmon, desun, demon, dr, dn,
	    de;

	dhi = -0.0025;
	dli = -0.0007;

	rsta = enorm8(xsta);
	sinphi = xsta[3] / rsta;
	cosphi = sqrt(xsta[1] * xsta[1] + xsta[2] * xsta[2]) / rsta;
	cos2phi = cosphi * cosphi - sinphi * sinphi;
	sinla = xsta[2] / cosphi / rsta;
	cosla = xsta[1] / cosphi / rsta;
	rmon = enorm8(xmon);
	rsun = enorm8(xsun);
	drsun = -3.0 * dhi * sinphi * cosphi * fac2sun * xsun[3] * (xsun[1] * sinla - xsun[2] * cosla) / rsun / rsun;
	drmon = -3.0 * dhi * sinphi * cosphi * fac2mon * xmon[3] * (xmon[1] * sinla - xmon[2] * cosla) / rmon / rmon;
	dnsun = -3.0 * dli * cos2phi * fac2sun * xsun[3] * (xsun[1] * sinla - xsun[2] * cosla) / rsun / rsun;
	dnmon = -3.0 * dli * cos2phi * fac2mon * xmon[3] * (xmon[1] * sinla - xmon[2] * cosla) / rmon / rmon;
	desun = -3.0 * dli * sinphi * fac2sun * xsun[3] * (xsun[1] * cosla + xsun[2] * sinla) / rsun / rsun;
	demon = -3.0 * dli * sinphi * fac2mon * xmon[3] * (xmon[1] * cosla + xmon[2] * sinla) / rmon / rmon;
	dr = drsun + drmon;
	dn = dnsun + dnmon;
	de = desun + demon;
	xcorsta[1] = dr * cosla * cosphi - de * sinla - dn * sinphi * cosla;
	xcorsta[2] = dr * sinla * cosphi + de * cosla - dn * sinphi * sinla;
	xcorsta[3] = dr * sinphi + dn * cosphi;

	return (1);
}

double enorm8(double *a) {

	//*** compute euclidian norm of a vector (of length 3)

	// double precision a(3)
	double enorm8_return;
	enorm8_return = sqrt(a[1] * a[1] + a[2] * a[2] + a[3] * a[3]);

	return (enorm8_return);
}

int sprod(double *x, double *y, double *scal, double *r1, double *r2) {

	//***  computation of the scalar-product of two vectors and their norms

	//***  input:   x(i),i=1,2,3  -- components of vector x
	//***           y(i),i=1,2,3  -- components of vector y
	//***  output:  scal          -- scalar product of x and y
	//***           r1,r2         -- lengths of the two vectors x and y

	// implicit double precision (a-h,o-z)
	// double  x(4),y(4);

	*r1 = sqrt(x[1] * x[1] + x[2] * x[2] + x[3] * x[3]);
	*r2 = sqrt(y[1] * y[1] + y[2] * y[2] + y[3] * y[3]);
	*scal = x[1] * y[1] + x[2] * y[2] + x[3] * y[3];

	return (1);
}

int moonxyz(double mjd, double fmjd, double *rm) {

	//*** get low-precision, geocentric coordinates for moon (ECEF)

	//*** input:  mjd/fmjd, is Modified Julian Date (and fractional) in GPS time
	//*** output: rm, is geocentric lunar position vector [m] in ECEF
	//*** 1."satellite orbits: models, methods, applications" montenbruck &
	// gill(2000)
	//*** section 3.3.2, pg. 72-73
	//*** 2."astronomy on the personal computer, 4th ed." montenbruck & pfleger
	//(2005)
	//*** section 3.2, pg. 38-39  routine MiniMoon

	// implicit double precision(a-h,o-z)
	// dimension rm(3)
	// common/stuff/rad,pi,pi2

	//*** use TT for lunar ephemerides

	double tsecgps, tsectt, fmjdtt, tjdtt, t, el0, el, elp, f, d, selond, q, selatd, rse, oblir, sselat, cselat, sselon, cselon,
	    t1, t2, t3;
	double ghar;

	tsecgps = fmjd * 86400.0;  //!*** GPS time (sec of day)
	tsectt = gps2tt(tsecgps);  //!*** TT  time (sec of day)
	fmjdtt = tsectt / 86400.0; //!*** TT  time (fract. day)

	//*** julian centuries since 1.5 january 2000 (J2000)
	//***   (note: also low precision use of mjd --> tjd)

	tjdtt = round(mjd) + fmjdtt + 2400000.50; //!*** Julian Date, TT
	t = (tjdtt - 2451545.0) / 36525.0;        //!*** julian centuries, TT

	//*** el0 -- mean longitude of Moon (deg)
	//*** el  -- mean anomaly of Moon (deg)
	//*** elp -- mean anomaly of Sun  (deg)
	//*** f   -- mean angular distance of Moon from ascending node (deg)
	//*** d   -- difference between mean longitudes of Sun and Moon (deg)

	//*** equations 3.47, p.72

	el0 = 218.316170 + (481267.880880 - 1.3972) * t;
	el = 134.962920 + 477198.867530 * t;
	elp = 357.525430 + 35999.049440 * t;
	f = 93.272830 + 483202.018730 * t;
	d = 297.850270 + 445267.111350 * t;

	//*** longitude w.r.t. equinox and ecliptic of year 2000

	selond = el0 //!*** eq 3.48, p.72
	         + 22640.0 / 3600.0 * sin((el) / rad) + 769.0 / 3600.0 * sin((el + el) / rad) -
	         4586.0 / 3600.0 * sin((el - d - d) / rad) + 2370.0 / 3600.0 * sin((d + d) / rad) -
	         668.0 / 3600.0 * sin((elp) / rad) - 412.0 / 3600.0 * sin((f + f) / rad) -
	         212.0 / 3600.0 * sin((el + el - d - d) / rad) - 206.0 / 3600.0 * sin((el + elp - d - d) / rad) +
	         192.0 / 3600.0 * sin((el + d + d) / rad) - 165.0 / 3600.0 * sin((elp - d - d) / rad) +
	         148.0 / 3600.0 * sin((el - elp) / rad) - 125.0 / 3600.0 * sin((d) / rad) - 110.0 / 3600.0 * sin((el + elp) / rad) -
	         55.0 / 3600.0 * sin((f + f - d - d) / rad);

	//*** latitude w.r.t. equinox and ecliptic of year 2000

	q = 412.0 / 3600.0 * sin((f + f) / rad) + 541.0 / 3600.0 * sin((elp) / rad);

	selatd = //!*** eq 3.49, p.72
	    +18520.0 / 3600.0 * sin((f + selond - el0 + q) / rad) - 526.0 / 3600.0 * sin((f - d - d) / rad) +
	    44.0 / 3600.0 * sin((el + f - d - d) / rad) - 31.0 / 3600.0 * sin((-el + f - d - d) / rad) -
	    25.0 / 3600.0 * sin((-el - el + f) / rad) - 23.0 / 3600.0 * sin((elp + f - d - d) / rad) +
	    21.0 / 3600.0 * sin((-el + f) / rad) + 11.0 / 3600.0 * sin((-elp + f - d - d) / rad);

	//*** distance from Earth center to Moon (m)

	rse = 385000.0 * 1000.0 //!*** eq 3.50, p.72
	      - 20905.0 * 1000.0 * cos((el) / rad) - 3699.0 * 1000.0 * cos((d + d - el) / rad) -
	      2956.0 * 1000.0 * cos((d + d) / rad) - 570.0 * 1000.0 * cos((el + el) / rad) +
	      246.0 * 1000.0 * cos((el + el - d - d) / rad) - 205.0 * 1000.0 * cos((elp - d - d) / rad) -
	      171.0 * 1000.0 * cos((el + d + d) / rad) - 152.0 * 1000.0 * cos((el + elp - d - d) / rad);

	// printf("%.15lf %.15lf %.15lf %.15lf
	// %.15lf\n",selond,selatd,rse,oblir,selond);

	//*** convert spherical ecliptic coordinates to equatorial cartesian

	//*** precession of equinox wrt. J2000   (p.71)

	selond = selond + 1.39720 * t; //!*** degrees

	//*** position vector of moon (mean equinox & ecliptic of J2000) (EME2000,
	// ICRF)
	//***                         (plus long. advance due to precession -- eq.
	// above)

	oblir = 23.439291110 / rad; //!*** obliquity of the J2000 ecliptic

	sselat = sin(selatd / rad);
	cselat = cos(selatd / rad);
	sselon = sin(selond / rad);
	cselon = cos(selond / rad);

	t1 = rse * cselon * cselat; //!*** meters          !*** eq. 3.51, p.72
	t2 = rse * sselon * cselat; //!*** meters          !*** eq. 3.51, p.72
	t3 = rse * sselat;          //!*** meters          !*** eq. 3.51, p.72

	rot1(-oblir, t1, t2, t3, rm); //!*** eq. 3.51, p.72

	//*** convert position vector of moon to ECEF  (ignore polar motion/LOD)

	getghar(mjd, fmjd, &ghar);           //!*** sec 2.3.1,p.33
	rot3(ghar, rm[1], rm[2], rm[3], rm); //!*** eq. 2.89, p.37

	return (1);
}

int rot1(double theta, double x, double y, double z, double *rm) {

	//*** rotate coordinate axes about 1 axis by angle of theta radians
	//*** x,y,z transformed into u,v,w

	// implicit double precision(a-h,o-z)

	double s, c;
	s = sin(theta);
	c = cos(theta);

	rm[1] = x;
	rm[2] = c * y + s * z;
	rm[3] = c * z - s * y;

	return (1);
}

int sunxyz(double mjd, double fmjd, double *rs) {

	//*** get low-precision, geocentric coordinates for sun (ECEF)

	//*** input, mjd/fmjd, is Modified Julian Date (and fractional) in GPS time
	//*** output, rs, is geocentric solar position vector [m] in ECEF
	//*** 1."satellite orbits: models, methods, applications" montenbruck &
	// gill(2000)
	//*** section 3.3.2, pg. 70-71
	//*** 2."astronomy on the personal computer, 4th ed." montenbruck & pfleger
	//(2005)
	//*** section 3.2, pg. 39  routine MiniSun

	double obe, sobe, cobe, opod, tsecgps, tsectt, fmjdtt, tjdtt, t, emdeg, em, em2, r, slond, slon, sslon, cslon;
	double rs1, rs2, rs3, ghar;

	//*** mean elements for year 2000, sun ecliptic orbit wrt. Earth

	obe = 23.43929111 / rad; //!*** obliquity of the J2000 ecliptic
	sobe = sin(obe);
	cobe = cos(obe);
	opod = 282.9400; //!*** RAAN + arg.peri.  (deg.)

	//*** use TT for solar ephemerides

	tsecgps = fmjd * 86400.0;  //!*** GPS time (sec of day)
	tsectt = gps2tt(tsecgps);  //!*** TT  time (sec of day)
	fmjdtt = tsectt / 86400.0; //!*** TT  time (fract. day)

	//*** julian centuries since 1.5 january 2000 (J2000)
	//***   (note: also low precision use of mjd --> tjd)

	tjdtt = round(mjd) + fmjdtt + 2400000.5; //!*** Julian Date, TT
	t = (tjdtt - 2451545.0) / 36525.0;       //!*** julian centuries, TT
	emdeg = 357.5256 + 35999.049 * t;        //!*** degrees
	em = emdeg / rad;                        //!*** radians
	em2 = em + em;                           //!*** radians

	//*** series expansions in mean anomaly, em   (eq. 3.43, p.71)

	r = (149.619 - 2.499 * cos(em) - 0.021 * cos(em2)) * 1.0e9; //!*** m.
	slond = opod + emdeg + (6892.0 * sin(em) + 72.0 * sin(em2)) / 3600.0;

	//*** precession of equinox wrt. J2000   (p.71)

	slond = slond + 1.39720 * t; //!*** degrees
	// fprintf(stderr,"%.12lf\n",em);

	//*** position vector of sun (mean equinox & ecliptic of J2000) (EME2000,
	// ICRF)
	//***                        (plus long. advance due to precession -- eq.
	// above)

	slon = slond / rad; //!*** radians
	sslon = sin(slon);
	cslon = cos(slon);

	rs1 = r * cslon;        //!*** meters             !*** eq. 3.46, p.71
	rs2 = r * sslon * cobe; //!*** meters             !*** eq. 3.46, p.71
	rs3 = r * sslon * sobe; //!*** meters             !*** eq. 3.46, p.71

	// printf("%lf %lf %lf\n",rs1,rs2,rs3);
	//*** convert position vector of sun to ECEF  (ignore polar motion/LOD)

	getghar(mjd, fmjd, &ghar);     //!*** sec 2.3.1,p.33
	rot3(ghar, rs1, rs2, rs3, rs); //!*** eq. 2.89, p.37

	return (1);
}

int rot3(double theta, double x, double y, double z, double *rs) {

	//*** rotate coordinate axes about 3 axis by angle of theta radians
	//*** x,y,z transformed into u,v,w

	// implicit double precision(a-h,o-z)

	double s, c;

	s = sin(theta);
	c = cos(theta);

	rs[1] = c * x + s * y;
	rs[2] = c * y - s * x;
	rs[3] = z;

	return (1);
}

double gps2tt(double tsec) {

	//*** convert tsec in GPS to tsec in TT

	double gps2tt_return;
	gps2tt_return = tsec + 51.1840; //!*** fixed offset

	return (gps2tt_return);
}

int getghar(double mjd, double fmjd, double *ghar) {

	//*** convert mjd/fmjd in GPS time to Greenwich hour angle (in radians)
	//*** "satellite orbits: models, methods, applications" montenbruck &
	// gill(2000)
	//*** section 2.3.1, pg. 33

	//*** need UT to get sidereal time  ("astronomy on the personal computer", 4th
	// ed)
	//***                               (pg.43, montenbruck & pfleger, springer,
	// 2005)
	double tsecgps, tsecutc, fmjdutc, d, ghad, i;

	tsecgps = fmjd * 86400.0;    //!*** GPS time (sec of day)
	tsecutc = gps2utc(tsecgps);  //!*** UTC time (sec of day)
	fmjdutc = tsecutc / 86400.0; //!*** UTC time (fract. day)

	//***** d = MJD - 51544.5d0                               !*** footnote
	d = (round(mjd) - 51544) + (fmjdutc - 0.5); //!*** days since J2000

	//*** greenwich hour angle for J2000  (12:00:00 on 1 Jan 2000)

	//***** ghad = 100.46061837504d0 + 360.9856473662862d0*d  !*** eq. 2.85
	//(+digits)
	ghad = 280.460618375040 + 360.98564736628620 * d; //!*** corrn.   (+digits)

	//**** normalize to 0-360 and convert to radians

	i = floor(ghad / 360.0);
	*ghar = (ghad - i * 360.0) / rad;
	while (*ghar > pi2) {
		*ghar = *ghar - pi2;
	}
	while (*ghar < 0.0) {
		*ghar = *ghar + pi2;
	}

	return (1);
}

double gps2utc(double tsec) {
	//*** convert tsec in GPS to tsec in UTC

	//*** GPS is ahead of UTC  (c.f. USNO)
	//*** UTC is behind GPS
	//*** gpsleap() is (so far) positive (and increasing)
	//*** so, must subtract gpsleap from GPS to get UTC
	double gps2utc_return;
	gps2utc_return = tsec - gpsleap(tsec);

	return (gps2utc_return);
}

double gpsleap(double tsec) {

	//*** return total leap seconds since GPS epoch 1980jan06

	//*** note: does **NOT** return the full TAI-UTC delta
	//*** input time is GPS seconds -- initialized by setjd0()
	//*** Y2K -- only functional between 1980jan06-00:00:00  (GPS time start)
	//***                            and hard-coded date

	// implicit double precision(a-h,o-z)

	//***** "Julian Date Converter"
	//***** http://aa.usno.navy.mil/data/docs/JulianDate.php
	//***** "Bulletin C"
	//***** http://hpiers.obspm.fr/eoppc/bul/bulc/bulletinc.dat
	//***** parameter(mjdhard=55196)            !*** cut-off date 2009dec31
	//***** parameter(mjdhard=55377)            !*** cut-off date 2010jun30
	//***** parameter(mjdhard=55561)            !*** cut-off date 2010dec31
	//***** parameter(mjdhard=55742)            !*** cut-off date 2011jun30
	//***** parameter(mjdhard=55926)            !*** cut-off date 2011dec31
	//***** parameter(mjdhard=56108)            !*** cut-off date 2012jun30
	//***** parameter(mjdhard=56292)            !*** cut-off date 2012dec31
	//***** parameter(mjdhard=56473)            !*** cut-off date 2013jun30
	//***** parameter(mjdhard=56657)            !*** cut-off date 2013dec31
	//***** parameter(mjdhard=56838)            !*** cut-off date 2014jun30
	//***** parameter(mjdhard=57022)            !*** cut-off date 2014dec31
	//***** parameter(mjdhard=57203)            !*** cut-off date 2015jun30
	//***** parameter(mjdhard=57387)            !*** cut-off date 2015dec31
	//***** parameter(mjdhard=57569)            !*** cut-off date 2016jun30
	//***** parameter(mjdhard=57753)            !*** cut-off date 2016dec31
	//***** parameter(mjdhard=57934)            !*** cut-off date 2017jun30
	//***** parameter(mjdhard=58118)            !*** cut-off date 2017dec31
	// parameter(mjdhard=58299)            !*** cut-off date 2018jun30
	// parameter(mjdhard=58664)            !*** cut-off date 2019jun30
	double mjdhard = 59030.0;                  // temporary new cut-off date until 2020jun30
	double ttsec, mjd0t, tai_utc, gpsleap_return;

	// save  /mjdoff/
	// common/mjdoff/mjd0
	//*** clone for tests (and do any rollover)

	ttsec = tsec;
	mjd0t = mjd0;

	while (ttsec >= 86400.0) {
		ttsec = ttsec - 86400.0;
		mjd0t = mjd0t + 1;
	}

	while (ttsec < 0.0) {
		ttsec = ttsec + 86400.0;
		mjd0t = mjd0t - 1;
	}

	//*** test date limits

	if (mjd0t > mjdhard) {
		fprintf(stderr, "FATAL ERROR --\n");
		fprintf(stderr, "exceeded cut-off date in gpsleap()\n");
		die("", "");
	}

	if (mjd0t < 44244.0) { // then             !*** 1980jan06
		die("FATAL ERROR --\n", "cut-off date underflow in gpsleap()\n");
	}

	//*** http://maia.usno.navy.mil/ser7/tai-utc.dat
	//*** 1980 JAN  1 =JD 2444239.5  TAI-UTC=  19.0s
	//*** 1981 JUL  1 =JD 2444786.5  TAI-UTC=  20.0s
	//*** 1982 JUL  1 =JD 2445151.5  TAI-UTC=  21.0s
	//*** 1983 JUL  1 =JD 2445516.5  TAI-UTC=  22.0s
	//*** 1985 JUL  1 =JD 2446247.5  TAI-UTC=  23.0s
	//*** 1988 JAN  1 =JD 2447161.5  TAI-UTC=  24.0s
	//*** 1990 JAN  1 =JD 2447892.5  TAI-UTC=  25.0s
	//*** 1991 JAN  1 =JD 2448257.5  TAI-UTC=  26.0s
	//*** 1992 JUL  1 =JD 2448804.5  TAI-UTC=  27.0s
	//*** 1993 JUL  1 =JD 2449169.5  TAI-UTC=  28.0s
	//*** 1994 JUL  1 =JD 2449534.5  TAI-UTC=  29.0s
	//*** 1996 JAN  1 =JD 2450083.5  TAI-UTC=  30.0s
	//*** 1997 JUL  1 =JD 2450630.5  TAI-UTC=  31.0s
	//*** 1999 JAN  1 =JD 2451179.5  TAI-UTC=  32.0s
	//*** 2006 JAN  1 =JD 2453736.5  TAI-UTC=  33.0s
	//*** 2009 JAN  1 =JD 2454832.5  TAI-UTC=  34.0s
	//*** 2012 JUL  1 =JD 2456109.5  TAI-UTC=  35.0s
	//*** 2015 JUL  1 =JD 2457204.5  TAI-UTC=  36.0s
	//*** 2017 JAN  1 =JD 2457754.5  TAI-UTC=  37.0s

	//*** test against newest leaps first

	tai_utc = 0.0;

	if (mjd0t >= 57754) { // then       !*** 2017 JAN 1 = 57754
		tai_utc = 37.0;
	}
	else if (mjd0t >= 57204) { // then       !*** 2015 JUL 1 = 57204
		tai_utc = 36.0;
	}
	else if (mjd0t >= 56109) { // then       !*** 2012 JUL 1 = 56109
		tai_utc = 35.0;
	}
	else if (mjd0t >= 54832) { // then       !*** 2009 JAN 1 = 54832
		tai_utc = 34.0;
	}
	else if (mjd0t >= 53736) { // then       !*** 2006 JAN 1 = 53736
		tai_utc = 33.0;
	}
	else if (mjd0t >= 51179) { // then       !*** 1999 JAN 1 = 51179
		tai_utc = 32.0;
	}
	else if (mjd0t >= 50630) { // then       !*** 1997 JUL 1 = 50630
		tai_utc = 31.0;
	}
	else if (mjd0t >= 50083) { // then       !*** 1996 JAN 1 = 50083
		tai_utc = 30.0;
	}
	else if (mjd0t >= 49534) { // then       !*** 1994 JUL 1 = 49534
		tai_utc = 29.0;
	}
	else if (mjd0t >= 49169) { // then       !*** 1993 JUL 1 = 49169
		tai_utc = 28.0;
	}
	else if (mjd0t >= 48804) { // then       !*** 1992 JUL 1 = 48804
		tai_utc = 27.0;
	}
	else if (mjd0t >= 48257) { // then       !*** 1991 JAN 1 = 48257
		tai_utc = 26.0;
	}
	else if (mjd0t >= 47892) { // then       !*** 1990 JAN 1 = 47892
		tai_utc = 25.0;
	}
	else if (mjd0t >= 47161) { // then       !*** 1988 JAN 1 = 47161
		tai_utc = 24.0;
	}
	else if (mjd0t >= 46247) { // then       !*** 1985 JUL 1 = 46247
		tai_utc = 23.0;
	}
	else if (mjd0t >= 45516) { // then       !*** 1983 JUL 1 = 45516
		tai_utc = 22.0;
	}
	else if (mjd0t >= 45151) { // then       !*** 1982 JUL 1 = 45151
		tai_utc = 21.0;
	}
	else if (mjd0t >= 44786) { // then       !*** 1981 JUL 1 = 44786
		tai_utc = 20.0;
	}
	else if (mjd0t >= 44239) { // then       !*** 1980 JAN 1 = 44239
		tai_utc = 19.0;
	}

	//*** should never get here

	else {
		die("FATAL ERROR --\n", "fell thru tests in gpsleap()\n");
	}

	//*** convert TAI-UTC into GPS leap seconds

	gpsleap_return = tai_utc - 19.0;

	return (gpsleap_return);
}

int civmjd(double iyr, double imo, double idy, double ihr, double imn, double sec, double *mjd, double *fmjd) {

	//*** convert civil date to modified julian date
	//*** only valid in range mar-1900 thru feb-2100     (leap year protocols)
	//*** ref: hofmann-wellenhof, 2nd ed., pg 34-35
	//*** operation confirmed against table 3.3 values on pg.34
	double m, it1, it2;
	int y;
	if (iyr < 1900)
		die("year should not be smaller than 1900", ""); // stop 34588

	if (imo <= 2) {
		y = iyr - 1.0;
		m = imo + 12.0;
	}
	else {
		y = iyr;
		m = imo;
	}
	it1 = floor(365.25 * y);
	it2 = floor(30.6001 * (m + 1.0));
	*mjd = it1 + it2 + idy - 679019.0;

	*fmjd = (3600.0 * ihr + 60.0 * imn + sec) / 86400.0;
	return (1);
}

int mjdciv(double mjd, double fmjd, double *iyr, double *imo, double *idy, double *ihr, double *imn, double *sec) {

	//*** convert modified julian date to civil date

	//*** imo in range 1-12, idy in range 1-31
	//*** only valid in range mar-1900 thru feb-2100
	//*** ref: hofmann-wellenhof, 2nd ed., pg 34-35
	//*** operation confirmed for leap years (incl. year 2000)

	double rjd, ia, ib, ic, id, ie, it1, it2, it3, tmp;

	rjd = mjd + fmjd + 2400000.5;
	ia = floor(rjd + 0.5);
	ib = floor(ia + 1537.0);
	ic = floor((ib - 122.1) / 365.25);
	id = floor(365.25 * ic);
	ie = floor((ib - id) / 30.6001);

	//*** the fractional part of a julian day is fractional mjd + 0.5
	//*** therefore, fractional part of julian day + 0.5 is fractional mjd

	it1 = floor(ie * 30.6001);
	*idy = floor(ib - id - it1 + fmjd);
	it2 = floor(ie / 14.0);
	*imo = ie - 1.0 - 12 * it2;
	it3 = floor((7.0 + *imo) / 10.0);
	*iyr = ic - 4715.0 - it3;

	tmp = fmjd * 24.0;
	*ihr = floor(tmp);
	tmp = (tmp - *ihr) * 60.0;
	*imn = floor(tmp);
	*sec = (tmp - *imn) * 60.0;

	return (1);
}

int geo2xyz(double gla, double glo, double eht, double *x, double *y, double *z) {
	// convert geodetic lat, long, ellip ht. to x,y,z
	double sla, cla, w, w2, en;

	sla = sin(gla);
	cla = cos(gla);
	w2 = 1.0 - e2 * sla * sla;
	w = sqrt(w2);
	en = a / w;

	*x = (en + eht) * cla * cos(glo);
	*y = (en + eht) * cla * sin(glo);
	*z = (en * (1.0 - e2) + eht) * sla;

	return (1);
}

/* below are functions by Kang Wang */

int isleapyr(int yr) {
	if ((yr % 4) != 0) {
		return 0;
	}
	else if ((yr % 100) != 0) {
		return 1;
	}
	else if ((yr % 400) != 0) {
		return 0;
	}
	else {
		return 1;
	}
}

int day2date(double year, double day, double *mon, double *id) {
	int yr, iday;
	int leapyr;
	// char tmp_str[200];

	int month = 0, day_of_month = 0;
	int day_of_yr[13];
	// char month_str[10];
	// char yr_str[10],yr_short[10];
	// char date_str[10];
	// char str_out[10];

	yr = (int)year;
	iday = (int)ceil(day);
	//  printf("yr = %d, day_of_yr = %d\n",yr,iday);

	leapyr = isleapyr(yr);
	if (leapyr == 1) {
		day_of_yr[0] = 0;
		day_of_yr[1] = 31;
		day_of_yr[2] = 60;
		day_of_yr[3] = 91;
		day_of_yr[4] = 121;
		day_of_yr[5] = 152;
		day_of_yr[6] = 182;
		day_of_yr[7] = 213;
		day_of_yr[8] = 244;
		day_of_yr[9] = 274;
		day_of_yr[10] = 305;
		day_of_yr[11] = 335;
		day_of_yr[12] = 366;
	}
	else {
		day_of_yr[0] = 0;
		day_of_yr[1] = 31;
		day_of_yr[2] = 59;
		day_of_yr[3] = 90;
		day_of_yr[4] = 120;
		day_of_yr[5] = 151;
		day_of_yr[6] = 181;
		day_of_yr[7] = 212;
		day_of_yr[8] = 243;
		day_of_yr[9] = 273;
		day_of_yr[10] = 304;
		day_of_yr[11] = 334;
		day_of_yr[12] = 365;
	}

	if (iday > day_of_yr[0] && iday <= day_of_yr[1]) {
		month = 1;
		day_of_month = iday;
	}
	else if (iday > day_of_yr[1] && iday <= day_of_yr[2]) {
		month = 2;
		day_of_month = iday - day_of_yr[1];
	}
	else if (iday > day_of_yr[2] && iday <= day_of_yr[3]) {
		month = 3;
		day_of_month = iday - day_of_yr[2];
	}
	else if (iday > day_of_yr[3] && iday <= day_of_yr[4]) {
		month = 4;
		day_of_month = iday - day_of_yr[3];
	}
	else if (iday > day_of_yr[4] && iday <= day_of_yr[5]) {
		month = 5;
		day_of_month = iday - day_of_yr[4];
	}
	else if (iday > day_of_yr[5] && iday <= day_of_yr[6]) {
		month = 6;
		day_of_month = iday - day_of_yr[5];
	}
	else if (iday > day_of_yr[6] && iday <= day_of_yr[7]) {
		month = 7;
		day_of_month = iday - day_of_yr[6];
	}

	else if (iday > day_of_yr[7] && iday <= day_of_yr[8]) {
		month = 8;
		day_of_month = iday - day_of_yr[7];
	}

	else if (iday > day_of_yr[8] && iday <= day_of_yr[9]) {
		month = 9;
		day_of_month = iday - day_of_yr[8];
	}
	else if (iday > day_of_yr[9] && iday <= day_of_yr[10]) {
		month = 10;
		day_of_month = iday - day_of_yr[9];
	}
	else if (iday > day_of_yr[10] && iday <= day_of_yr[11]) {
		month = 11;
		day_of_month = iday - day_of_yr[10];
	}
	else if (iday > day_of_yr[11] && iday <= day_of_yr[12]) {
		month = 12;
		day_of_month = iday - day_of_yr[11];
	}

	//  sprintf(yr_str,"%d",yr);

	//  memset(yr_short, '\0', sizeof(yr_short));

	//  yr_str[4]='\0';
	//  strcpy(yr_str,yr);
	// sprintf(month_str,"%02d",month);
	// sprintf(date_str,"%02d",day_of_month);
	//  strncpy(yr_short,yr_str+2,2);
	//  strcat(str_out,yr_short);
	//  strcat(str_out,yr_str);
	//  strcat(str_out,month_str);
	//  strcat(str_out,month_str);
	//  strcat(str_out,date_str);
	//  fprintf(stdout,"%s\n",str_out);
	// printf("%04d%02d%02d\n",yr,month,day_of_month);

	*mon = month;
	*id = day_of_month;
	return (1);
}
