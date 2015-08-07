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
#include "metadata/Planet.hh"


const double
Planet::getMass() const
{
	return this->mass;
}

const double
Planet::getSpinRate() const
{
	return this->spinRate;
}

const double 
Planet::getSemiMajorAxis() const
{
	return this->semiMajorAxis;
}

const double 
Planet::getSemiMinorAxis() const
{
	return this->semiMinorAxis;
}

const double
Planet::getEccentricity() const
{
	return this->eccentricity;
}

const double
Planet::getFlattening() const
{
	return this->flattening;
}

const double
Planet::getGM() const
{
	// We should just store GM, since it is known more precisely than either G or M
	double G = 6.673e-11; // Gravitational Constant m^3 kg^-1 s^-2
	return this->mass*G;
}

const double
Planet::getEccentricitySquared() const
{
	return this->eccentricity*this->eccentricity;
}

// WGS84 Earth parameters
const Planet Planet::Earth = Planet(5.9736e24,7.29211573052e-5,6378137,6356752.314245,0.081819190842621,(1.0/298.25722356));
