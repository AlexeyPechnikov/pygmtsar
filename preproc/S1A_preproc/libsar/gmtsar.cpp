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
#include <cmath>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include "boost/date_time/posix_time/posix_time.hpp" 
#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/algorithm/string.hpp"
#include "metadata/StateVector.hh"
#include "metadata/Planet.hh"
#include "Doppler.hh"
#include "gmtsar.hh"
#include "utilities.hh"

void GMTSAR::writeEphemeris(const std::string filename,Orbit *orbit,Scene *scene)
{
    std::stringstream ss;

    using namespace boost::gregorian;

    std::vector<StateVector> *stateVectors = orbit->getStateVectors();

/*  do the loop once to get the number of vectors and their time spacing */

    int jj = 0;
    double dt = 0.;
    double t00 = 0.;
    double t0 = 0.;
    double t1 = 0.;
    for (std::vector<StateVector>::iterator itr = stateVectors->begin(); itr != stateVectors->end(); ++itr){
        StateVector sv = *itr;
        t0 = t1;
        t1 = (double)(sv.getTime().time_of_day().total_microseconds()/1e6);
        if(jj == 0) t00 = t1;
        jj++;
    }
    dt = t1 - t0;

/* compute the julian day */

    int year = scene->getSensingStart().date().year();
    int month = scene->getSensingStart().date().month();
    int day = scene->getSensingStart().date().day();
    date tacq = date(year,month,day);
    date tjan1 = date(year,1,1);
    days julday = tacq -  tjan1;

    ss << jj << " "<< year << " "<< julday << std::setprecision(12) << " " << t00 << " "<< dt << " "<< "\n";
    for (std::vector<StateVector>::iterator itr = stateVectors->begin(); itr != stateVectors->end(); ++itr)
    {
        StateVector sv = *itr;
	double t = (double)(sv.getTime().time_of_day().total_microseconds()/1e6);
	double *pos = sv.getPosition();
	double *vel = sv.getVelocity();
	ss << " "<< year << " "<< julday << std::setprecision(12) << " "<< t << " " << 
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
void GMTSAR::writeResourceFile(const std::string filename,Scene *scene,
		Instrument *instrument,Orbit *orbit, SingleBandImage<T> *image)
{
	std::stringstream ss;
        using namespace boost::gregorian;

	// Calculate the azimuth spacing
	boost::posix_time::ptime midTime = getSensingMid(scene);

	double velocity = calculateVelocity(orbit,midTime);
	double azimuthPixelSize = instrument->getAzimuthPixelSize();
	if ( azimuthPixelSize == 0.0) {
		azimuthPixelSize = velocity/instrument->getPulseRepetitionFrequency();
	}

	// Set up the formatting for the date times
	boost::posix_time::time_facet *facet = new boost::posix_time::time_facet();
	ss.imbue(std::locale(ss.getloc(),facet));

	// Calculate the center line
        //int centerLine = (int)(image->getHeight()/2);
	double firstLineUTC = (double)(scene->getSensingStart().time_of_day().total_microseconds()/1e6);
	//double lastLineUTC = (double)(scene->getSensingStop().time_of_day().total_microseconds()/1e6);
	//double centerLineUTC = (lastLineUTC+firstLineUTC)/2.0;

	// make sure the number of lines is divisable by 4 
	int nlines = image->getHeight() - image->getHeight()%4;

	// compute the julian day 
	int year = scene->getSensingStart().date().year();
	int month = scene->getSensingStart().date().month();
	int day = scene->getSensingStart().date().day();
	date tacq = date(year,month,day);
	date tjan = date(year,1,1);
	days julday = tacq - tjan;

	//* compute the start and end time 
        double SC_clock_start = (double)1000.*year + julday.days() +  firstLineUTC/86400.;
        double prf = instrument->getPulseRepetitionFrequency();
        double SC_clock_stop =  SC_clock_start + nlines/(prf*86400.);
        //double SC_clock_stop2 = (double)1000.*year + julday.days() + lastLineUTC/86400.;

	// get the stem name
	std::string stemname = filename.substr(0,filename.size()-3);
        std::string dir = scene->getPassDirection().substr(0,1);

	// determine the spacecraft identity and set the number
	int SC_identity = 9999;
	std::string SCname = instrument->getPlatform()->getMission().substr(0,1);
	std::string TSX ("T");
	if(SCname.compare(TSX) == 0) SC_identity = 7;
	std::string CSK ("C");
	if(SCname.compare(CSK) == 0) SC_identity = 8;
	std::string RS2 ("R");
	if(SCname.compare(RS2) == 0) SC_identity = 9;
        std::string S1A ("S");
	if(SCname.compare(S1A) == 0) SC_identity = 10; 

	// parameters needed for both raw and SLC
	ss << "first_line\t\t= 1\n";
	ss << "st_rng_bin\t\t= 1\n";
	ss << "nlooks\t\t\t= 1\n";
	ss << "rshift\t\t\t= 0\n";
	ss << "ashift\t\t\t= 0\n";
	ss << "sub_int_r\t\t= 0.0\n";
	ss << "sub_int_a\t\t= 0.0\n";
	ss << "stretch_r\t\t= 0.0\n";
	ss << "stretch_a\t\t= 0.0\n";
	ss << "a_stretch_r\t\t= 0.0\n";
	ss << "a_stretch_a\t\t= 0.0\n";
	ss << "first_sample\t\t= 0\n";
	ss << "dtype\t\t\t= a\n";
	ss << "rng_samp_rate\t\t= " << std::setprecision(7) << instrument->getRangeSamplingFrequency() << "\n";
	ss << "SC_identity\t\t= " << SC_identity << "\n";
	ss << "radar_wavelength\t= " << std::setprecision(9) << instrument->getWavelength() << "\n";
	ss << "chirp_slope\t\t= " << std::setprecision(9) << instrument->getChirpSlope() << "\n";
	ss << "pulse_dur\t\t= " << instrument->getPulseLength() << "\n";
	ss << "I_mean\t\t\t= "<< scene->getIBias() << "\n";
	ss << "Q_mean\t\t\t= "<< scene->getQBias() << "\n";
	ss << "PRF\t\t\t= " << std::setprecision(9) << instrument->getPulseRepetitionFrequency() << "\n";
	ss << "near_range\t\t= " << std::setprecision(9) << scene->getStartingRange() << "\n";
	ss << "equatorial_radius\t= " << Planet::Earth.getSemiMajorAxis() << "\n";
	ss << "polar_radius\t\t= " << Planet::Earth.getSemiMinorAxis() << "\n";
	ss << "orbdir\t\t\t= " << dir << "\n";
	ss << "input_file\t\t= " << stemname << "raw\n";
	ss << "led_file\t\t= " << stemname << "LED\n";
	ss << "SLC_file\t\t= " << stemname << "SLC\n";
	ss << "SLC_scale\t\t= 1.0\n";
	ss << "SC_clock_start\t\t= " << std::setprecision(17) << SC_clock_start << "\n";
	ss << "SC_clock_stop\t\t= " << std::setprecision(17) << SC_clock_stop  << "\n";

	// parameters needed for the raw data
	ss << "Flip_iq\t\t\t= n\n";
	ss << "deskew\t\t\t= n\n";
	ss << "offset_video\t\t= n\n";
	ss << "caltone\t\t\t= 0.0\n";
	ss << "rm_az_band\t\t= 0.0\n";
	ss << "rm_rng_band\t\t= 0.2\n";
	ss << "rng_spec_wgt\t\t= 1.0\n";
	ss << "scnd_rng_mig\t\t= n\n";
	//ss << "az_res\t\t\t= " << azimuthPixelSize << "\n";
	ss << "az_res\t\t\t= 3.0\n";
	ss << "antenna_side\t\t= " << instrument->getAntennaSide() << "\n";
	ss << "fdd1\t\t\t= 0.0\n";
	ss << "fddd1\t\t\t= 0.0\n";

	std::string procLevel = scene->getProcessingLevel();
	std::string SLC ("SLC");
	if(procLevel.compare(SLC) == 0) {
	// parameters only for SLC
	  ss << "bytes_per_line\t\t= " << 4*image->getWidth() << "\n";
          ss << "good_bytes_per_line\t= " << 4*image->getWidth() << "\n";
	  ss << "num_lines\t\t= " << nlines << "\n";
	  ss << "nrows\t\t\t= " << nlines << "\n";
	  ss << "num_valid_az\t\t= " << nlines << "\n";
	  ss << "num_patches\t\t= 1\n";
	  ss << "num_rng_bins\t\t= " << image->getWidth() << "\n";
	  ss << "chirp_ext\t\t= 0\n";
	}
	else {
	// parameters only for raw
	  ss << "bytes_per_line\t\t= " << 2*image->getWidth() << "\n";
          ss << "good_bytes_per_line\t= " << 2*image->getWidth() << "\n";
	  ss << "num_lines\t\t= " << nlines << "\n";
	  ss << "nrows\t\t\t= 8192\n";
	  ss << "num_valid_az\t\t= 6400\n";
	  int num_patches = (int)((float)nlines/6400. + 0.5);
	  ss << "num_patches\t\t= " << num_patches << "\n";
	  // add some padding to the number of range bins in increments of 1200 and make the chirp extension 2400
	  ss << "chirp_ext\t\t= 2400\n";
	  int num_rng_bins = image->getWidth();
	  int nrng = 2 + (int)(num_rng_bins/1200.);
	  num_rng_bins = 1200 * nrng;
	  ss << "num_rng_bins\t\t= " << num_rng_bins << "\n";
	}
	
	std::ofstream rscFile (filename.c_str(),std::ofstream::binary);
	rscFile.write(ss.str().c_str(), ss.str().length());
	rscFile.close();
}

template void GMTSAR::writeResourceFile<std::complex<unsigned char> >(const std::string filename,Scene *scene,Instrument *instrument,Orbit *orbit,SingleBandImage<std::complex<unsigned char> > *image);
template void GMTSAR::writeResourceFile<std::complex<short> >(const std::string filename,Scene *scene,Instrument *instrument,Orbit *orbit,SingleBandImage<std::complex<short> > *image);
