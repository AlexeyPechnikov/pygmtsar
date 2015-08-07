/*
 *  Copyright (C) 2012 Walter M. Szeliga
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include "Doppler.hh"
#include "metadata/Orbit.hh"
#include "metadata/StateVector.hh"
#include "metadata/Instrument.hh"
#include "metadata/Scene.hh"
#include "metadata/Planet.hh"
#include "boost/date_time/posix_time/posix_time.hpp" 

/**
 * Calculate latitude,longitude,height from cartesian coordinates and
 * ellipsoidal parameters
 * 
 * @param xyz cartesian coordinates
 * @param e2 the eccentricity squared
 * @param a the semi-major axis
 * @return the latitude,longitude,height
 */
std::vector<double> 
xyz_to_llh(double *xyz,double e2,double a)
{
    std::vector<double> llh(3,0.0);

    double q2 = 1.0/(1.0-e2);
    double q = sqrt(q2);
    double q3 = q2 - 1.0;
    double b = a*sqrt(1.0 - e2);
    
    llh[1] = atan2(xyz[1],xyz[0]);

    double p = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
    double tant = (xyz[2]/p)*q;
    double theta = atan(tant);
    tant = (xyz[2] + q3*b*pow(sin(theta),3.0))/(p - e2*a*pow(cos(theta),3.0));
    llh[0] = atan(tant);
    double re = a/sqrt(1.0 - e2*sin(llh[0])*sin(llh[0]));
    llh[2] = p/cos(llh[0]) - re;

    llh[0] = llh[0]*180.0/M_PI;
    llh[1] = llh[1]*180.0/M_PI;

    return llh;
}

/**
 * Calculate the absolute velocity of a satellite at a given time.
 *
 * @param orbit the orbit of the satellite
 * @param time the time at which to calculate the velocity
 * @return the total velocity of the satellite
 */
double
calculateVelocity(Orbit *orbit,boost::posix_time::ptime time)
{
	StateVector sv = orbit->interpolate(time,"linear");
	double *vel = sv.getVelocity();
	return sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
}

/**
 * Calculate the height at mid-swath
 *
 * @param orbit the satellite orbit
 * @param time the time at which to calculate the height
 * @return the height at mid-swath
 */
double
calculateHeight(Orbit *orbit,boost::posix_time::ptime time)
{

    StateVector sv0 = orbit->interpolate(time,"linear");
    double a = Planet::Earth.getSemiMajorAxis();
    double e = Planet::Earth.getEccentricity();
    double e2 = e*e;

    std::vector<double> LLH = xyz_to_llh(sv0.getPosition(),e2,a);
    
    double height = LLH[2];

    return height;
}

/**
 * Calculate the absolute change in satellite height above the ellipsoid
 * over some time interval.
 *
 * @note the satellite is assumed to be orbiting around a WGS84 Earth
 *
 * @param orbit the satellite orbit
 * @param startTime the beginning of the time interval
 * @param midTime the end of the time interval
 * @param the change in satellite height above the ellipsoid
 */
double
calculateHeightDt(Orbit *orbit,boost::posix_time::ptime startTime,boost::posix_time::ptime midTime)
{
    StateVector sv0 = orbit->interpolate(startTime,"linear");
    StateVector sv1 = orbit->interpolate(midTime,"linear");

    double a = Planet::Earth.getSemiMajorAxis();
    double e = Planet::Earth.getEccentricity();
    double e2 = e*e;

    std::vector<double> startLLH = xyz_to_llh(sv0.getPosition(),e2,a);
    std::vector<double> midLLH = xyz_to_llh(sv1.getPosition(),e2,a);

    double startHeight = startLLH[2];
    double midHeight = midLLH[2];

    double heightDt = (midHeight - startHeight)/((midTime-startTime).total_microseconds()/1e6);
    return heightDt;
}

/**
 * Calculate the squint angle of the satellite in degrees
 *
 */
double
calculateSquint(Doppler *doppler,Scene *scene, Instrument *instrument, double velocity)
{
	double height = scene->getSatelliteHeight();
	double startingRange = scene->getStartingRange();
	double wavelength = instrument->getWavelength();
	double prf = instrument->getPulseRepetitionFrequency();

	if (height > startingRange)
	{
		std::cerr << "Spacecraft height too large (" << height << ">" << startingRange << ")" << std::endl;
		return 0.0;
	}
	double sinTheta = sqrt(1.0 - (height/startingRange)*(height/startingRange));
	double fd = doppler->getConstant()*prf;
	double sinSquint = fd/(2.0*velocity*sinTheta)*wavelength;
	if (sinSquint*sinSquint > 1.0)
	{
		std::cerr << "Error in squint calculation" << std::endl;
		return 0.0;
	}
	double squint = atan2(sinSquint,sqrt(1.0-sinSquint*sinSquint));
	squint = squint*180.0/M_PI;

	return squint;
}

/**
 * Calculate the radius of the Earth at a given point specified
 * by latitude, longitude, and height using the formula from 
 * Bomford, "Geodesy" 1971, eqs A.53 and A.55 pg 565.
 *
 * @param llh
 * @param heading
 * @param a
 * @param e2
 * @return the radius of curvature of the Earth 
 */
double
radiusOfCurvature(std::vector<double> llh, double heading, double e2, double a)
{
    double lat = llh[0]*M_PI/180.0;
    heading = heading*M_PI/180.0;

    double east = a/sqrt(1.0 - e2*sin(lat)*sin(lat));
    double north = (a*(1.0 - e2))/pow(sqrt(1.0 - e2*sin(lat)*sin(lat)),3.0);

    double dir = (east*north)/(east*cos(heading)*cos(heading) + north*sin(heading)*sin(heading));

    return dir;
}

/**
 * A simplified wrapper for calculating the radius of curvature
 * at the center of the scene, for ROI_PAC.  
 *
 * @param orbit
 * @param midTime
 * @return the radius of curvature
 */
double
calculateEarthRadius(Orbit *orbit, boost::posix_time::ptime midTime)
{
    StateVector sv1 = orbit->interpolate(midTime,"linear");

    double a = Planet::Earth.getSemiMajorAxis();
    double e = Planet::Earth.getEccentricity();
    double e2 = e*e;

    std::vector<double> midLLH = xyz_to_llh(sv1.getPosition(),e2,a);
    double rcurv = radiusOfCurvature(midLLH,0.0,e2,a);

    return rcurv;
}

boost::posix_time::ptime
getSensingMid(const Scene *scene)
{
	boost::posix_time::ptime startTime = scene->getSensingStart();
	boost::posix_time::ptime stopTime = scene->getSensingStop();
	long duration = (stopTime-startTime).total_microseconds();
	boost::posix_time::ptime midTime = startTime + boost::posix_time::microseconds((long)duration/2.0);

	return midTime;
}
