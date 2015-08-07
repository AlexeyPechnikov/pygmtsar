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
#include "metadata/Platform.hh"

Platform::Platform() {}

Platform::~Platform() {}

void
Platform::setMission(const std::string mission)
{
	this->mission = mission;
}

std::string
Platform::getMission() const
{
	return this->mission;
}

/*
void
Platform::setPlanet(const Planet planet)
{
	this->planet = planet;
}

const Planet 
Platform::getPlanet()
{
	return this->planet;
}*/
