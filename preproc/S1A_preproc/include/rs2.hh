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

#ifndef RS2_HH
#define RS2_HH 1

#include <complex>
#include <string>
#include "metadata/Orbit.hh"
#include "metadata/Instrument.hh"
#include "metadata/Scene.hh"
#include "image/SingleBandImage.hh"
#include "tinyxml.h"

class RS2
{
	private:
		std::string directory;
		TiXmlDocument *document;
		Instrument *instrument;
		Scene *scene;
		Orbit *orbit;
		void populatePlatform();
		void populateInstrument();
		void populateScene();
		void populateOrbit();
		std::string readElement(const char *xpath);
	public:
		RS2(const std::string filename);
		~RS2();
		SingleBandImage<std::complex<short> > *extractSlc_short_Image(const std::string outFile);
		SingleBandImage<std::complex<float> > *extractSlcImage(const std::string outFile);
		Instrument *getInstrument();
		Scene *getScene();
		Orbit *getOrbit();
};

#endif
