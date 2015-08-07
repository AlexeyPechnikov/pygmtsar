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
#include "metadata/Instrument.hh"
#include <iostream>
#include <sstream>

Instrument::Instrument() : azimuthPixelSize(0.0)
{
	this->platform = new Platform();
}

Instrument::~Instrument() {}

void
Instrument::setPlatform(Platform *platform)
{
	this->platform = platform;
}

void
Instrument::setWavelength(double wavelength)
{
	this->wavelength = wavelength;
}

void
Instrument::setIncidenceAngle(double incidenceAngle)
{
	this->incidenceAngle = incidenceAngle;
}

void
Instrument::setPulseRepetitionFrequency(double pulseRepetitionFrequency)
{
	this->pulseRepetitionFrequency = pulseRepetitionFrequency;
}

void
Instrument::setRangePixelSize(double rangePixelSize)
{
	this->rangePixelSize = rangePixelSize;
}

void
Instrument::setAzimuthPixelSize(double azimuthPixelSize)
{
	this->azimuthPixelSize = azimuthPixelSize;
}

void
Instrument::setPulseLength(double pulseLength)
{
	this->pulseLength = pulseLength;
}

void
Instrument::setChirpSlope(double chirpSlope)
{
	this->chirpSlope = chirpSlope;
}

void
Instrument::setAntennaLength(double antennaLength)
{
	this->antennaLength = antennaLength;
}

void
Instrument::setAntennaSide(int antennaSide)
{
	this->antennaSide = antennaSide;
}

void
Instrument::setRangeSamplingFrequency(double rangeSamplingFrequency)
{
	this->rangeSamplingFrequency = rangeSamplingFrequency;
}

Platform *
Instrument::getPlatform()
{
	return this->platform;
}

double
Instrument::getWavelength() const
{
	return this->wavelength;
}

double
Instrument::getIncidenceAngle() const
{
	return this->incidenceAngle;
}

double
Instrument::getPulseRepetitionFrequency() const
{
	return this->pulseRepetitionFrequency;
}

double
Instrument::getRangePixelSize() const
{
	return this->rangePixelSize;
}

double
Instrument::getAzimuthPixelSize() const
{
	return this->azimuthPixelSize;
}

double
Instrument::getPulseLength() const
{
	return this->pulseLength;
}

double
Instrument::getChirpSlope() const
{
	return this->chirpSlope;
}

double
Instrument::getAntennaLength() const
{
	return this->antennaLength;
}

int
Instrument::getAntennaSide() const
{
	return this->antennaSide;
}

double
Instrument::getRangeSamplingFrequency() const
{
	return this->rangeSamplingFrequency;
}

std::string
Instrument::toString()
{
	std::stringstream ss;
	ss << "Wavelength " << wavelength << "\n";
	ss << "PRF " << pulseRepetitionFrequency << "\n";
	ss << "Pulse Length " << pulseLength << "\n";
	ss << "Range Pixel Size " << rangePixelSize << "\n";
	ss << "Antenna Length " << antennaLength << "\n";
	ss << "Chirp Slope " << chirpSlope << "\n";
	ss << "Antenna Side " << antennaSide;
	return ss.str();
}
