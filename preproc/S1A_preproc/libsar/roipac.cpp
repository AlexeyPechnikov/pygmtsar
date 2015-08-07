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
#include <sstream>
#include <fstream>
#include "boost/date_time/posix_time/posix_time.hpp" 
#include "boost/algorithm/string.hpp"
#include "metadata/StateVector.hh"
#include "metadata/Planet.hh"
#include "Doppler.hh"
#include "roipac.hh"
#include "utilities.hh"

void ROIPAC::writeEphemeris(const std::string filename,Orbit *orbit)
{
    std::stringstream ss;

    std::vector<StateVector> *stateVectors = orbit->getStateVectors();
    for (std::vector<StateVector>::iterator itr = stateVectors->begin(); itr != stateVectors->end(); ++itr)
    {
        StateVector sv = *itr;
	double t = (double)(sv.getTime().time_of_day().total_microseconds()/1e6);
	double *pos = sv.getPosition();
	double *vel = sv.getVelocity();
	ss << t << " " << std::setprecision(12) << 
		pos[0] << " " <<
		pos[1] << " " <<
		pos[2] << " " <<
		vel[0] << " " <<
		vel[1] << " " <<
		vel[2] << "\n";
    }

    std::ofstream ephFile (filename.c_str(),std::ofstream::binary);
    ephFile.write(ss.str().c_str(), ss.str().length());
    ephFile.close();
}

template <typename T>
void ROIPAC::writeResourceFile(const std::string filename,Scene *scene,
		Instrument *instrument,Orbit *orbit,Doppler *doppler,
		SingleBandImage<T> *image)
{
	std::stringstream ss;

	// Calculate the scene mid-time
	boost::posix_time::ptime midTime = getSensingMid(scene);

	double velocity = calculateVelocity(orbit,midTime);
	double height = calculateHeight(orbit,midTime);
	double heightDt = calculateHeightDt(orbit,scene->getSensingStart(),midTime);
	double squint = calculateSquint(doppler,scene,instrument,velocity);
	double earthRadius = calculateEarthRadius(orbit,midTime);
	double azimuthPixelSize = instrument->getAzimuthPixelSize();
	if ( azimuthPixelSize == 0.0) {
		azimuthPixelSize = velocity/instrument->getPulseRepetitionFrequency();
	}

	// Set up the formatting for the date times
	boost::posix_time::time_facet *facet = new boost::posix_time::time_facet();
	ss.imbue(std::locale(ss.getloc(),facet));

	// Calculate the center line
        int centerLine = (int)(image->getHeight()/2);
	double firstLineUTC = (double)(scene->getSensingStart().time_of_day().total_microseconds()/1e6);
	double lastLineUTC = (double)(scene->getSensingStop().time_of_day().total_microseconds()/1e6);
	double centerLineUTC = (lastLineUTC+firstLineUTC)/2.0;


	//ss << "FIRST_FRAME\n";
	facet->format("%Y%m%d%H%M%S%F");
	ss << "FIRST_FRAME_SCENE_CENTER_TIME\t" << scene->getSensingStart() << "\n"; // Need to switch this to center line
	ss << "FIRST_FRAME_SCENE_CENTER_LINE\t" << centerLine << "\n";
	facet->format("%Y%m%d");
	ss << "DATE\t\t\t\t" << scene->getSensingStart() << "\n";
	ss << "FIRST_LINE_YEAR\t\t\t" << scene->getSensingStart().date().year() << "\n";
	ss << "FIRST_LINE_MONTH_OF_YEAR\t" << scene->getSensingStart().date().month() << "\n";
	ss << "FIRST_LINE_DAY_OF_MONTH\t\t" << scene->getSensingStart().date().day() << "\n";
	ss << "FIRST_CENTER_HOUR_OF_DAY\t" << scene->getSensingStart().time_of_day().hours() << "\n";
	ss << "FIRST_CENTER_MN_OF_HOUR\t\t" << scene->getSensingStart().time_of_day().minutes() << "\n";
	ss << "FIRST_CENTER_S_OF_MN\t\t" << scene->getSensingStart().time_of_day().seconds()<< "\n";
	ss << "FIRST_CENTER_MS_OF_S\t\t" << (int)(scene->getSensingStart().time_of_day().fractional_seconds()/1e3) << "\n";
	ss << "PROCESSING_FACILITY\t\t" << scene->getProcessingFacility() << "\n";
	ss << "PROCESSING_VERSION\t\t" << scene->getProcessingLevel() << "\n";
	ss << "PLATFORM\t\t\t" << instrument->getPlatform()->getMission() << "\n";
	ss << "BEAM\t\t\t\t" << scene->getBeam() << "\n";
	ss << "ORBIT_NUMBER\t\t\t" << scene->getOrbitNumber() << "\n";
	ss << "STARTING_RANGE\t\t\t" << std::setprecision(9) << scene->getStartingRange() << "\n";
	ss << "RANGE_PIXEL_SIZE\t\t" << instrument->getRangePixelSize() << "\n";
	ss << "PRF\t\t\t\t" << instrument->getPulseRepetitionFrequency() << "\n";
	ss << "ANTENNA_SIDE\t\t\t" << instrument->getAntennaSide() << "\n";
	ss << "ANTENNA_LENGTH\t\t\t" << instrument->getAntennaLength() << "\n";
	ss << "FILE_START\t\t\t1\n";
	ss << "FILE_LENGTH\t\t\t" << image->getHeight() << "\n";
	ss << "XMIN\t\t\t\t0\n";
	// Kind of a hack to account for making slc files
	if (boost::algorithm::ends_with(filename,".slc.rsc")) {
	    ss << "XMAX\t\t\t\t" << (image->getWidth()) << "\n";
	    ss << "WIDTH\t\t\t\t" << image->getWidth() << "\n";
	} else {
	    ss << "XMAX\t\t\t\t" << 2*(image->getWidth()) << "\n";
	    ss << "WIDTH\t\t\t\t" << 2*(image->getWidth()) << "\n";
	}
	ss << "YMIN\t\t\t\t0\n";
	ss << "YMAX\t\t\t\t" << image->getHeight() << "\n";
	ss << "RANGE_SAMPLING_FREQUENCY\t" << instrument->getRangeSamplingFrequency() << "\n";
	ss << "PLANET_GM\t\t\t" << Planet::Earth.getGM() << "\n";
	ss << "PLANET_SPINRATE\t\t\t" << Planet::Earth.getSpinRate() << "\n";
	ss << "FIRST_LINE_UTC\t\t\t" << std::setprecision(9) <<  firstLineUTC << "\n";
	ss << "CENTER_LINE_UTC\t\t\t" << std::setprecision(9) << centerLineUTC << "\n";
	ss << "LAST_LINE_UTC\t\t\t" << std::setprecision(9) << lastLineUTC << "\n";
	ss << "HEIGHT\t\t\t\t" << height << "\n";
	ss << "HEIGHT_DT\t\t\t" << heightDt << "\n";
	ss << "VELOCITY\t\t\t" << velocity << "\n";
	ss << "EQUATORIAL_RADIUS\t\t" << Planet::Earth.getSemiMajorAxis() << "\n";
	ss << "ECCENTRICITY_SQUARED\t\t" << Planet::Earth.getEccentricitySquared() << "\n";
	ss << "EARTH_RADIUS\t\t\t" << std::setprecision(12) << earthRadius << "\n";
	ss << "ORBIT_DIRECTION\t\t\t" << scene->getPassDirection() << "\n";
	ss << "WAVELENGTH\t\t\t" << instrument->getWavelength() << "\n";
	ss << "PULSE_LENGTH\t\t\t" << instrument->getPulseLength() << "\n";
	ss << "CHIRP_SLOPE\t\t\t" << instrument->getChirpSlope() << "\n";
	ss << "I_BIAS\t\t\t\t" << scene->getIBias() << "\n";
	ss << "Q_BIAS\t\t\t\t" << scene->getQBias() << "\n";
	ss << "DOPPLER_RANGE0\t\t\t" << doppler->getConstant() << "\n";
	ss << "DOPPLER_RANGE1\t\t\t" << doppler->getLinear() << "\n";
	ss << "DOPPLER_RANGE2\t\t\t" << doppler->getQuadratic() << "\n";
	ss << "DOPPLER_RANGE3\t\t\t" << 0.0 << "\n";
	ss << "SQUINT\t\t\t\t" << squint << "\n";
	ss << "AZIMUTH_PIXEL_SIZE\t\t" << azimuthPixelSize << "\n";
	ss << "DELTA_LINE_UTC\t\t\t" << 1.0/instrument->getPulseRepetitionFrequency() << "\n";
	ss << "ROI_PAC_VERSION\t\t\t3\n";

	std::ofstream rscFile (filename.c_str(),std::ofstream::binary);
	rscFile.write(ss.str().c_str(), ss.str().length());
	rscFile.close();
}

template void ROIPAC::writeResourceFile<std::complex<unsigned char> >(const std::string filename,Scene *scene,Instrument *instrument,Orbit *orbit,Doppler *doppler,SingleBandImage<std::complex<unsigned char> > *image);
template void ROIPAC::writeResourceFile<std::complex<float> >(const std::string filename,Scene *scene,Instrument *instrument,Orbit *orbit,Doppler *doppler,SingleBandImage<std::complex<float> > *image);
