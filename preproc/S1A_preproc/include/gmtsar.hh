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
#ifndef ROIPAC_HH
#define ROIPAC_HH 1

#include <string>
#include <complex>
#include "metadata/Orbit.hh"
#include "metadata/Instrument.hh"
#include "metadata/Scene.hh"
#include "image/SingleBandImage.hh"

namespace GMTSAR {

void writeEphemeris(const std::string filename, Orbit *orbit, Scene *scene);
template <typename T>
void writeResourceFile(const std::string filename,Scene *scene,
		Instrument *instrument,Orbit *orbit, SingleBandImage<T> *image);


};

#endif
