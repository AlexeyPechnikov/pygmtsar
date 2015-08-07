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
#ifndef ORBIT_HH
#define ORBIT_HH 1

#include <string>
#include <vector>
#include "StateVector.hh"
#include "boost/date_time/posix_time/posix_time.hpp" 


/**
 * An Orbit holds StateVectors containing the position and velocity 
 * of a satellite over time.
 *
 * @author Walter Szeliga
 */
class Orbit
{
	private:
		boost::posix_time::ptime minTime;
		boost::posix_time::ptime maxTime;
		std::vector<StateVector> *stateVectors; //!< A vector of StateVector objects
		std::string referenceFrame; //!< A text string indicating the orbital reference frame
		std::string quality; //!< A text string indicating the orbit quality
		std::string source; //!< A text string to indicate the source of the orbital data
                StateVector linear(boost::posix_time::ptime time);
	public:
		Orbit();
		~Orbit();
		void setQuality(std::string quality);
		void setReferenceFrame(std::string referenceFrame);
		void setSource(std::string source);
		void addStateVector(StateVector sv);
		std::string getQuality() const;
		std::string getReferenceFrame() const;
		std::string getSource() const;
		std::vector<StateVector> *getStateVectors(); // Ugh, this is hideous, I should write a proper iterator at some point
		Orbit trim(boost::posix_time::ptime start, boost::posix_time::ptime stop);
                StateVector interpolate(boost::posix_time::ptime time, std::string method);
};

#endif
