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

/*
 * Modification History:
 *
 * DATE: 2014.06.16             by: Xiaohua Xu
 * MODIFICATION:
 * converted from rs2.hh for purpose of processing sentinel1-a data
 *
 */

#ifndef S1A_HH
#define S1A_HH 1

#include <complex>
#include <string>
#include "metadata/Orbit.hh"
#include "metadata/Instrument.hh"
#include "metadata/Scene.hh"
#include "image/SingleBandImage.hh"
#include "tinyxml.h"

class S1A
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
		void populateBurst();
		std::string readElement(const char *xpath);
	public:
		S1A(const std::string filename);
		~S1A();
		SingleBandImage<std::complex<short> > *extractSlc_short_Image(const std::string outFile,const std::string inFile);
		SingleBandImage<std::complex<float> > *extractSlcImage(const std::string outFile);
		Instrument *getInstrument();
		Scene *getScene();
		Orbit *getOrbit();
};

#endif
