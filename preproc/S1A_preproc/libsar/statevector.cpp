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
#include "metadata/StateVector.hh"

StateVector::StateVector() {}
StateVector::~StateVector() {}

void
StateVector::setTime(boost::posix_time::ptime time)
{
	this->time = time;
}

void
StateVector::setPosition(double pos[])
{
	for (int i=0;i<3;i++)
	{
		this->position[i] = pos[i];
	}
}

void
StateVector::setVelocity(double vel[])
{
	for (int i=0;i<3;i++)
	{
		this->velocity[i] = vel[i];
	}
}

boost::posix_time::ptime
StateVector::getTime()
{
	return this->time;
}

double *
StateVector::getPosition()
{
	return this->position;
}

double *
StateVector::getVelocity()
{
	return this->velocity;
}

std::string
StateVector::toString()
{
    std::stringstream ss;
    ss << time << " ";
    ss << position[0] << " " << position[1] << " " << position[2] << " ";
    ss << velocity[0] << " " << velocity[1] << " " << velocity[2];

    return ss.str();
}
