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
#include "metadata/Scene.hh"

Scene::Scene() {}

Scene::~Scene() {}

void
Scene::setStartingRange(double startingRange)
{
    this->startingRange = startingRange;
}

void
Scene::setIBias(double iBias)
{
    this->iBias = iBias;
}

void
Scene::setQBias(double qBias)
{
    this->qBias = qBias;
}

void
Scene::setPassDirection(std::string passDirection)
{
    this->passDirection = passDirection;
}

void
Scene::setOrbitNumber(int orbitNumber)
{
    this->orbitNumber = orbitNumber;
}

void
Scene::setProcessingFacility(std::string processingFacility)
{
    this->processingFacility = processingFacility;
}

void
Scene::setProcessingLevel(std::string processingLevel)
{
    this->processingLevel = processingLevel;
}

void
Scene::setPolarization(std::string polarization)
{
    this->polarization = polarization;
}

void
Scene::setBeam(std::string beam)
{
    this->beam = beam;
}

void
Scene::setSatelliteHeight(double satelliteHeight)
{
    this->satelliteHeight = satelliteHeight;
}

/*void
Scene::setNumberOfLines(int numberOfLines)
{
    this->numberOfLines = numberOfLines;
}

void
Scene::setNumberOfSamples(int numberOfSamples)
{
    this->numberOfSamples = numberOfSamples;
}*/

void
Scene::setSensingStart(boost::posix_time::ptime sensingStart)
{
	this->sensingStart = sensingStart;
}

void
Scene::setSensingStop(boost::posix_time::ptime sensingStop)
{
	this->sensingStop = sensingStop;
}

// Getters

boost::posix_time::ptime
Scene::getSensingStart() const
{
	return this->sensingStart;
}

boost::posix_time::ptime
Scene::getSensingStop() const
{
	return this->sensingStop;
}

double
Scene::getStartingRange() const
{
    return this->startingRange;
}

double
Scene::getSatelliteHeight() const
{
    return this->satelliteHeight;
}

double
Scene::getIBias() const
{
    return this->iBias;
}

double
Scene::getQBias() const
{
    return this->qBias;
}

std::string
Scene::getProcessingFacility() const
{
    return this->processingFacility;
}

std::string
Scene::getProcessingLevel() const
{
    return this->processingLevel;
}

std::string
Scene::getPolarization() const
{
    return this->polarization;
}

std::string
Scene::getBeam() const
{
    return this->beam;
}

int
Scene::getOrbitNumber() const
{
    return this->orbitNumber;
}

/*int
Scene::getNumberOfLines()
{
    return this->numberOfLines;
}

int
Scene::getNumberOfSamples()
{
    return this->numberOfSamples;
}*/

std::string
Scene::getPassDirection() const
{
	return this->passDirection;
}

std::string
Scene::toString()
{
	std::stringstream ss;
	ss << "Starting Range " << startingRange << "\n";
	ss << "Satellite Height " << satelliteHeight << "\n";
	ss << "Orbit Number " << orbitNumber << "\n";
	ss << "Sensing Start " << sensingStart << "\n";
	ss << "Sensing Stop " << sensingStop;

	return ss.str();
}
