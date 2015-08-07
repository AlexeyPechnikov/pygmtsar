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
#ifndef UTILITIES_HH
#define UTILITIES_HH 1
#include <vector>
#include "Doppler.hh"
#include "metadata/Orbit.hh"
#include "boost/date_time/posix_time/posix_time.hpp" 

std::vector<double> xyz_to_llh(double *xyz,double e2,double a);
double calculateVelocity(Orbit *orbit,boost::posix_time::ptime time);
double calculateHeight(Orbit *orbit,boost::posix_time::ptime time);
double calculateHeightDt(Orbit *orbit,boost::posix_time::ptime startTime,boost::posix_time::ptime midTime);
double calculateSquint(Doppler *doppler,Scene *scene, Instrument *instrument, double velocity);
double calculateEarthRadius(Orbit *orbit, boost::posix_time::ptime midTime);
double radiusOfCurvature(double *llh, double heading, double e2, double a);
boost::posix_time::ptime getSensingMid(const Scene *scene);

#endif
