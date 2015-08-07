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
#ifndef STATEVECTOR_HH
#define STATEVECTOR_HH 1

#include <string>
#include "boost/date_time/posix_time/posix_time.hpp" 

/**
 * A StateVector holds position and velocity information at a given epoch.
 *
 * @author Walter Szeliga
 */
class StateVector
{
	private:
		boost::posix_time::ptime time;
		double position[3];
		double velocity[3];
	public:
		StateVector();
		~StateVector();
		boost::posix_time::ptime getTime();
		double *getPosition();
		double *getVelocity();
		void setTime(boost::posix_time::ptime time);
		void setPosition(double *position);
		void setVelocity(double *velocity);
                std::string toString();
};

#endif
