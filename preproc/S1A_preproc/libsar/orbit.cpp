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
#include <iostream>
#include "metadata/Orbit.hh"

Orbit::Orbit() 
{
	this->stateVectors = new std::vector<StateVector>();
	this->minTime = boost::posix_time::ptime(boost::posix_time::max_date_time);
	this->maxTime = boost::posix_time::ptime(boost::posix_time::min_date_time);
}

Orbit::~Orbit()
{
	delete this->stateVectors;
}

void
Orbit::setQuality(std::string quality)
{
	this->quality = quality;
}

void
Orbit::setReferenceFrame(std::string referenceFrame)
{
	this->referenceFrame = referenceFrame;
}

void
Orbit::setSource(std::string source)
{
	this->source = source;
}

void
Orbit::addStateVector(StateVector sv)
{
	this->stateVectors->push_back(sv);
	if (sv.getTime() < this->minTime)
	{
		this->minTime = sv.getTime();
	}
	if (sv.getTime() > this->maxTime)
	{
		this->maxTime = sv.getTime();
	}
}

std::vector<StateVector> *
Orbit::getStateVectors()
{
	return this->stateVectors;
}

std::string
Orbit::getQuality() const
{
	return this->quality;
}

std::string
Orbit::getSource() const
{
	return this->source;
}

std::string
Orbit::getReferenceFrame() const
{
	return this->referenceFrame;
}

Orbit
Orbit::trim(boost::posix_time::ptime start, boost::posix_time::ptime stop)
{
	Orbit newOrbit = Orbit();
	newOrbit.setSource(this->source);
	newOrbit.setReferenceFrame(this->referenceFrame);
	newOrbit.setQuality(this->quality);
	for (std::vector<StateVector>::iterator itr = this->stateVectors->begin(); itr!=this->stateVectors->end(); ++itr)
	{
		StateVector sv = *itr;
		if ((sv.getTime() > start) && (sv.getTime() < stop))
		{
			newOrbit.addStateVector(sv);
		}

	}

	return newOrbit;
}

StateVector
Orbit::interpolate(boost::posix_time::ptime time, std::string method)
{
 StateVector sv;

 if (method.compare("linear") == 0)
 {
   sv = linear(time);
 }
 else
 {
   std::cerr << "Unsupported interpolation method " << method << std::endl;
 }

 return sv;
}

StateVector
Orbit::linear(boost::posix_time::ptime time)
{
  if (this->stateVectors->size() < 2)
  {
    std::cerr << "Fewer than 2 state vectors present in orbit, cannot interpolate" << std::endl;
    return StateVector();
  };

  double position [] = {0.0, 0.0, 0.0};
  double velocity [] = {0.0, 0.0, 0.0};

  for (std::vector<StateVector>::iterator itr1 = this->stateVectors->begin();itr1!=this->stateVectors->end();++itr1)
  {
      StateVector sv1 = *itr1;
      double tmp=1.0;
      for (std::vector<StateVector>::iterator itr2 = this->stateVectors->begin();itr2!=this->stateVectors->end();++itr2)
      {
          StateVector sv2 = *itr2;
          if (sv1.getTime() == sv2.getTime()) {continue;}
          double numerator = (double)((sv2.getTime()-time).total_microseconds()/1e6);
          double denominator = (double)((sv2.getTime() - sv1.getTime()).total_microseconds()/1e6);
          tmp = tmp*(numerator)/(denominator);
      }
      for (int i=0;i<3;i++)
      {
          position[i] = position[i] + sv1.getPosition()[i]*tmp;
          velocity[i] = velocity[i] + sv1.getVelocity()[i]*tmp;
      }
  }

  StateVector newSV = StateVector();
  newSV.setTime(time);
  newSV.setPosition(position);
  newSV.setVelocity(velocity);
  return newSV;
}
