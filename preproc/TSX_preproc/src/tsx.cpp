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
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>
#include <libgen.h>
#include "xpath_static.h"
#include "tsx.hh"

TSX::TSX(const std::string filename)
{
  boost::filesystem::path pathname(filename);
  directory = pathname.parent_path().string();
  // Empty paths should be converted to '.'
  if (pathname.parent_path().empty()) {
	  directory = ".";
  }
  document = new TiXmlDocument();
  bool load_ok = this->document->LoadFile(filename.c_str());
  if (!load_ok) {
    fprintf(stderr,"\n ****** xml file not loaded \n \n");
    exit(-1);
  }

  this->instrument = new Instrument();
  this->orbit = new Orbit();
  this->populatePlatform();
  this->populateInstrument();
  this->populateScene();
  this->populateOrbit();
}

TSX::~TSX() 
{
	delete this->document;
}

void TSX::populatePlatform()
{
   std::string satelliteID;

   // Xpath for the satelliteID is /productInfo/missionInfo/mission
   TIXML_STRING mission = TinyXPath::S_xpath_string(this->document->RootElement(),"/level1Product/productInfo/missionInfo/mission/text()");
   satelliteID = std::string(mission.c_str());

   // Extract the satellite name
   this->instrument->getPlatform()->setMission(satelliteID);
}

void TSX::populateInstrument()
{
    double c = 299792458.0;
    int antennaSide;
    TIXML_STRING tmp = TinyXPath::S_xpath_string(this->document->RootElement(),
	    "/level1Product/productInfo/acquisitionInfo/lookDirection/text()");
    std::string pointingDirection(tmp.c_str());
    if (pointingDirection.compare("RIGHT") != 0)
    {
        antennaSide = 1;
    }
    else
    {
        antennaSide = -1;
    }

    double frequency = TinyXPath::d_xpath_double(this->document->RootElement(),
            "/level1Product/instrument/radarParameters/centerFrequency/text()");
    double wavelength = c/frequency;
    double prf = TinyXPath::d_xpath_double(this->document->RootElement(),
            "/level1Product/productSpecific/complexImageInfo/commonPRF/text()");
    double pulseLength = TinyXPath::d_xpath_double(this->document->RootElement(),
	    "/level1Product/processing/processingParameter/rangeCompression/chirps/referenceChirp/pulseLength/text()");
    pulseLength = pulseLength/1e9; // picoseconds to seconds, the pulse length appears to be encoded as picoseconds
    double chirpPulseBandwidth = TinyXPath::d_xpath_double(this->document->RootElement(),
	    "/level1Product/processing/processingParameter/rangeCompression/chirps/referenceChirp/pulseBandwidth/text()");
    double chirpSlope = chirpPulseBandwidth/pulseLength;
    double antennaLength = 4.8;
    double rowSpacing = TinyXPath::d_xpath_double(this->document->RootElement(),
            "/level1Product/productInfo/imageDataInfo/imageRaster/rowSpacing/text()");
    //double azimuthPixelSize = TinyXPath::d_xpath_double(this->document->RootElement(),
    //        "/level1Product/productInfo/imageDataInfo/imageRaster/azimuthResolution/text()");
    // Range Sampling Frequency is stored in /level1Product/instrument/settings/RSF
    //double rangeSamplingFrequency = 1/(2.0*rowSpacing);
    double rangeSamplingFrequency = 1/rowSpacing;
    double rangePixelSize = (c*rowSpacing/2.0);
    tmp = TinyXPath::S_xpath_string(this->document->RootElement(),
		    "/level1Product/processing/processingParameter/rangeCompression/chirps/referenceChirp/chirpSlope/text()");
    std::string chirpDirection(tmp.c_str());
    if (chirpDirection.compare("DOWN") == 0)
    {
        chirpSlope = -1.0*chirpSlope;
    }

    this->instrument->setAntennaSide(antennaSide);
    this->instrument->setWavelength(wavelength);
    this->instrument->setPulseRepetitionFrequency(prf);
    this->instrument->setPulseLength(pulseLength);
    this->instrument->setChirpSlope(chirpSlope);
    this->instrument->setAntennaLength(antennaLength);
    this->instrument->setRangeSamplingFrequency(rangeSamplingFrequency);
    this->instrument->setRangePixelSize(rangePixelSize);
    //this->instrument->setAzimuthPixelSize(azimuthPixelSize);

}

void TSX::populateScene()
{
    double c = 299792458.0;
    // Sampling Window Start Time is stored in /level1Product/instrument/settings/settingRecord/echoWindowPosition
    // Starting Range
    double startingRange = TinyXPath::d_xpath_double(this->document->RootElement(),
		    "/level1Product/productInfo/sceneInfo/rangeTime/firstPixel/text()");
    startingRange = startingRange*c/2.0;
    // Satellite Height
    double satelliteHeight = 0.0;
    // Processing facility
    TIXML_STRING tmp = TinyXPath::S_xpath_string(this->document->RootElement(),
		"/level1Product/productInfo/generationInfo/level1ProcessingFacility/text()");
    std::string processingFacility(tmp.c_str());
    std::string processingLevel = "SLC";
    // Orbit number
    int orbitNumber = TinyXPath::i_xpath_int(this->document->RootElement(),
		    "/level1Product/productInfo/missionInfo/absOrbit/text()");
    // Orbit direction
    tmp = TinyXPath::S_xpath_string(this->document->RootElement(),
		    "/level1Product/productInfo/missionInfo/orbitDirection/text()");
    std::string orbitDirection(tmp.c_str());
    // Polarization
    tmp = TinyXPath::S_xpath_string(this->document->RootElement(),
		    "/level1Product/instrument/settings/polLayer/text()");
    std::string polarization(tmp.c_str());

    // Set up the date time formatter
    boost::posix_time::time_input_facet *facet = new boost::posix_time::time_input_facet("%Y-%m-%dT%H:%M:%S%F%Z");

    // Sensing Start
    tmp = TinyXPath::S_xpath_string(this->document->RootElement(),
		    "/level1Product/productInfo/sceneInfo/start/timeUTC/text()");
    std::stringstream sensingStartString(tmp.c_str());
    sensingStartString.imbue(std::locale(sensingStartString.getloc(),facet));
    boost::posix_time::ptime sensingStart;
    sensingStartString >> sensingStart;
    // Sensing Stop
    tmp = TinyXPath::S_xpath_string(this->document->RootElement(),
		    "/level1Product/productInfo/sceneInfo/stop/timeUTC/text()");
    std::stringstream sensingStopString(tmp.c_str());
    sensingStopString.imbue(std::locale(sensingStopString.getloc(),facet));
    boost::posix_time::ptime sensingStop;
    sensingStopString >> sensingStop;

    this->scene = new Scene();
    this->scene->setIBias(0.0);
    this->scene->setQBias(0.0);
    this->scene->setStartingRange(startingRange);
    this->scene->setSatelliteHeight(satelliteHeight);
    this->scene->setProcessingFacility(processingFacility);
    this->scene->setProcessingLevel(processingLevel);
    this->scene->setPolarization(polarization);
    this->scene->setOrbitNumber(orbitNumber);
    this->scene->setPassDirection(orbitDirection);
    this->scene->setSensingStart(sensingStart);
    this->scene->setSensingStop(sensingStop);
}

void TSX::populateOrbit()
{
    double pos[3];
    double vel[3];

    // Get the orbit quality
    TIXML_STRING tmp = TinyXPath::S_xpath_string(this->document->RootElement(),
		    "/level1Product/platform/orbit/orbitHeader/accuracy/text()");
    std::string quality(tmp.c_str());
    this->orbit->setQuality(quality);

    // Get the number of state vectors
    TinyXPath::xpath_processor xp_proc(this->document->RootElement(),"/level1Product/platform/orbit/stateVec");
    unsigned int count = xp_proc.u_compute_xpath_node_set();

    // Set up the date time formatter
    boost::posix_time::time_input_facet *facet = new boost::posix_time::time_input_facet("%Y-%m-%dT%H:%M:%S%F");
    std::stringstream epochString;
    epochString.imbue(std::locale(epochString.getloc(),facet));

    // Get the state vectors
    for (unsigned int i=0;i<count;i++) {
	// Parse XML
        TiXmlNode *stateVector = xp_proc.XNp_get_xpath_node(i);
	TIXML_STRING tmp = TinyXPath::S_xpath_string(stateVector,"//timeUTC/text()");
	pos[0] = TinyXPath::d_xpath_double(stateVector,"//posX/text()");
	pos[1] = TinyXPath::d_xpath_double(stateVector,"//posY/text()");
	pos[2] = TinyXPath::d_xpath_double(stateVector,"//posZ/text()");
	vel[0] = TinyXPath::d_xpath_double(stateVector,"//velX/text()");
	vel[1] = TinyXPath::d_xpath_double(stateVector,"//velY/text()");
	vel[2] = TinyXPath::d_xpath_double(stateVector,"//velZ/text()");

        // Translate time stamp
	epochString.str(tmp.c_str());
        boost::posix_time::ptime epoch;
        epochString >> epoch;

	// Create State Vector
	StateVector sv = StateVector();
	sv.setTime(epoch);
	sv.setPosition(pos);
	sv.setVelocity(vel);
	this->orbit->addStateVector(sv);
    }
    
    
}

bool
TSX::checkImagingMode()
{
	TIXML_STRING tmp = TinyXPath::S_xpath_string(this->document->RootElement(),
			"/level1Product/productInfo/acquisitionInfo/imagingMode/text()");
	std::string imagingMode(tmp.c_str());
	if (imagingMode.compare("SM") != 0) {
		std::cerr << "This does not appear to be a strip-map mode acquisition" << std::endl;
		return false;
	}
        this->scene->setBeam(imagingMode);
	return true;
}

SingleBandImage<std::complex<short> > *TSX::extractSlc_short_Image(const std::string output)
{
  // Check that this is a strip-map product
  if (!checkImagingMode()) {
	  return NULL;
  }

  // Get the image filename
  TIXML_STRING tmp = TinyXPath::S_xpath_string(this->document->RootElement(),
		  "/level1Product/productComponents/imageData/file/location/path/text()");
  std::string relativePath(tmp.c_str());
  tmp = TinyXPath::S_xpath_string(this->document->RootElement(),
		  "/level1Product/productComponents/imageData/file/location/filename/text()");
  std::string filename(tmp.c_str());

  std::string input = this->directory + "/" + relativePath + "/" + filename;
  std::cout << "Image File location " << input << std::endl;
  try {
      cosar = new Cosar(input,output);
      cosar->parse();
  } catch (const char *ex) {
	  std::cerr << ex << std::endl;
  }

  int rows = TinyXPath::i_xpath_int(this->document->RootElement(),
		  "/level1Product/productInfo/imageDataInfo/imageRaster/numberOfRows/text()");
  int cols = TinyXPath::i_xpath_int(this->document->RootElement(),
		  "/level1Product/productInfo/imageDataInfo/imageRaster/numberOfColumns/text()");
  SingleBandImage<std::complex<short> > *image = 
	  new SingleBandImage<std::complex<short> >(output.c_str(),"r",cols,rows);

  return image;
}

Instrument *
TSX::getInstrument()
{
    return this->instrument;
}

Scene *
TSX::getScene()
{
    return this->scene;
}

Orbit *
TSX::getOrbit()
{
    return this->orbit;
}
