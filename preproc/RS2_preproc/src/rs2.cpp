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
#include "ezETAProgressBar.hpp"
#include "xpath_static.h"
#include "tiffio.h"
#include "rs2.hh"

RS2::RS2(const std::string filename)
{
    boost::filesystem::path pathname(filename);
    directory = pathname.parent_path().string();
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

RS2::~RS2()
{
    delete this->document;
}

void RS2::populatePlatform()
{
  TIXML_STRING tmp = TinyXPath::S_xpath_string(this->document->RootElement(),
                  "/product/sourceAttributes/satellite/text()");
  std::string satelliteID(tmp.c_str());

  this->instrument->getPlatform()->setMission(satelliteID);
}

void RS2::populateInstrument()
{
    double c = 299792458.0;
    double frequency = TinyXPath::d_xpath_double(this->document->RootElement(), 
                    "/product/sourceAttributes/radarParameters/radarCenterFrequency/text()");
    double wavelength = c/frequency;

    double prf = TinyXPath::d_xpath_double(this->document->RootElement(),
		    "/product/sourceAttributes/radarParameters/pulseRepetitionFrequency/text()");

//* need to divided the prf by 2 if this is quadpole and multiply by 2 if this is fine beam
     prf = prf/2.;

    double pulseLength = TinyXPath::d_xpath_double(this->document->RootElement(),
		    "/product/sourceAttributes/radarParameters/pulseLength[1]/text()");
    double pulseBandwidth = TinyXPath::d_xpath_double(this->document->RootElement(),
		    "/product/sourceAttributes/radarParameters/pulseBandwidth[1]/text()");
    double chirpSlope = pulseBandwidth/pulseLength;
    double rangePixelSize = TinyXPath::d_xpath_double(this->document->RootElement(),
		    "/product/imageAttributes/rasterAttributes/sampledPixelSpacing/text()");
/*  this may need to be divided by 4 instead of 2. */
    double rangeSamplingFrequency = c/(2.0*rangePixelSize);

    this->instrument->setAntennaSide(-1);
    this->instrument->setWavelength(wavelength);
    this->instrument->setPulseRepetitionFrequency(prf);
    this->instrument->setPulseLength(pulseLength);
    this->instrument->setChirpSlope(chirpSlope);
    this->instrument->setAntennaLength(15.0);
    this->instrument->setRangeSamplingFrequency(rangeSamplingFrequency);
    this->instrument->setRangePixelSize(rangePixelSize);
}

void RS2::populateScene()
{
    double c = 299792458.0;
    double height = TinyXPath::d_xpath_double(this->document->RootElement(),
		    "/product/imageGenerationParameters/sarProcessingInformation/satelliteHeight/text()");
    double firstSampleTime = TinyXPath::d_xpath_double(this->document->RootElement(),
		    "/product/imageGenerationParameters/slantRangeToGroundRange/slantRangeTimeToFirstRangeSample/text()");
    double startingRange = firstSampleTime*c/2.0;
    std::string facility = this->readElement("/product/imageGenerationParameters/generalProcessingInformation/processingFacility/text()");
    //std::string version = this->readElement("/product/imageGenerationParameters/generalProcessingInformation/softwareVersion/text()");
    std::string version = "SLC";
    std::string passDirection = this->readElement("/product/sourceAttributes/orbitAndAttitude/orbitInformation/passDirection/text()");

    // get orbit direction
    std::string DIR = passDirection.substr(0,1);
    std::string ASC ("A");
    // if this is an ascending pass then use the last time as the start time 
      std::string dataStartString = this->readElement("/product/imageGenerationParameters/sarProcessingInformation/zeroDopplerTimeFirstLine/text()");
    if(DIR.compare(ASC) == 0){
      std::string dataStartString = this->readElement("/product/imageGenerationParameters/sarProcessingInformation/zeroDopplerTimeLastLine/text()");
    } 

    // Set up the date time formatter
    boost::posix_time::time_input_facet *facet = new boost::posix_time::time_input_facet("%Y-%m-%dT%H:%M:%S%FZ");
    std::stringstream sensingStartString(dataStartString.c_str());
    sensingStartString.imbue(std::locale(sensingStartString.getloc(),facet));
    boost::posix_time::ptime sensingStart;
    sensingStartString >> sensingStart;

    double prf = TinyXPath::d_xpath_double(this->document->RootElement(),
		    "/product/sourceAttributes/radarParameters/pulseRepetitionFrequency/text()");
    int lines = TinyXPath::i_xpath_int(this->document->RootElement(),
		    "/product/imageAttributes/rasterAttributes/numberOfLines/text()");
    double msec = 1e6* double(lines)/prf;
    boost::posix_time::ptime sensingStop = sensingStart + boost::posix_time::microseconds(msec);
    std::string polarization = this->readElement("/product/sourceAttributes/radarParameters/polarizations/text()");
    std::string beam = this->readElement("/product/sourceAttributes/beamModeMnemonic/text()");

    this->scene = new Scene();
    this->scene->setIBias(0.0);
    this->scene->setQBias(0.0);
    this->scene->setStartingRange(startingRange);
    this->scene->setSatelliteHeight(height);
    this->scene->setProcessingFacility(facility);
    this->scene->setProcessingLevel(version);
    this->scene->setPolarization(polarization);
    this->scene->setBeam(beam);
    this->scene->setOrbitNumber(0); // No orbit number information is present in the XML data file
    this->scene->setPassDirection(passDirection);
    this->scene->setSensingStart(sensingStart);
    this->scene->setSensingStop(sensingStop);
}

void RS2::populateOrbit()
{
    double pos[3];
    double vel[3];

    // Get the number of state vectors
    TinyXPath::xpath_processor xp_proc(this->document->RootElement(),
		    "/product/sourceAttributes/orbitAndAttitude/orbitInformation/stateVector");
    unsigned int count = xp_proc.u_compute_xpath_node_set();

    // Set up the date time formatter
    boost::posix_time::time_input_facet *facet = new boost::posix_time::time_input_facet("%Y-%m-%dT%H:%M:%S%F");
    std::stringstream epochString;
    epochString.imbue(std::locale(epochString.getloc(),facet));

    for (unsigned int i=0;i<count;i++) {
        // Parse XML
	TiXmlNode *stateVector = xp_proc.XNp_get_xpath_node(i);
	TIXML_STRING tmp = TinyXPath::S_xpath_string(stateVector,"//timeStamp/text()");
	pos[0] = TinyXPath::d_xpath_double(stateVector,"//xPosition/text()");
	pos[1] = TinyXPath::d_xpath_double(stateVector,"//yPosition/text()");
	pos[2] = TinyXPath::d_xpath_double(stateVector,"//zPosition/text()");
	vel[0] = TinyXPath::d_xpath_double(stateVector,"//xVelocity/text()");
	vel[1] = TinyXPath::d_xpath_double(stateVector,"//yVelocity/text()");
	vel[2] = TinyXPath::d_xpath_double(stateVector,"//zVelocity/text()");

	// Translate the time stamp
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

SingleBandImage<std::complex<short> > *RS2::extractSlc_short_Image(const std::string filename)
{
    int i,j;
    //float *float_buf;
    uint32 width,height,widthi;
    FILE *out;
    SingleBandImage<std::complex<short> > *image; 

    // Get the image filename
    TIXML_STRING tmp = TinyXPath::S_xpath_string(this->document->RootElement(),
		    "/product/imageAttributes/fullResolutionImageData/text()");
    std::string input(tmp.c_str());
    std::string imageFile = this->directory + "/" + input;

    TIFFSetWarningHandler(NULL);
    TIFF *tif = TIFFOpen(imageFile.c_str(),"r");

    if (!tif) {return NULL;}

    uint16 *buf;
    uint32 row;

    //out = fopen(filename.c_str(),"wb");
    TIFFGetField(tif,TIFFTAG_IMAGEWIDTH,&widthi);
    TIFFGetField(tif,TIFFTAG_IMAGELENGTH,&height);
    // make sure the width is divisible by 4
    width = widthi - widthi%4;
    try {
        image = new SingleBandImage<std::complex<short> >(filename.c_str(),"w",width,height);
    } catch (std::exception &e) {
	    std::cerr << e.what() << std::endl;
    }

    // get orbit direction
    std::string DIR = scene->getPassDirection().substr(0,1);
    std::string ASC ("A");

    buf = (uint16 *)_TIFFmalloc(TIFFScanlineSize(tif));
    std::complex<short> *cpx_buf = new std::complex<short>[width];

    ez::ezETAProgressBar pg((unsigned int)height/100);
    pg.start();

// if this is ascending then flip the image top to bottom
    if(DIR.compare(ASC) == 0){
      for (row = 0; row < height; row++) {
        if ((row%100) == 0) {++pg;}
        TIFFReadScanline(tif, buf, row);
        for (i=0,j=0;i<2*width;i+=2,j++) {
            short real = (short)buf[i];
            short imag = (short)buf[i+1];
            std::complex<short> cpx_val((short)real,(short)imag);
            cpx_buf[j] = cpx_val;
          }
        image->setRow((height-row-1),cpx_buf);
      }
    } else {

// this must be descending so reverse the lines
      for (row = 0; row < height; row++) {
        if ((row%100) == 0) {++pg;}
        TIFFReadScanline(tif, buf, row);
          for (i=0,j=0;i<2*width;i+=2,j++) {
            short real = (short)buf[i];
            short imag = (short)buf[i+1];
            std::complex<short> cpx_val((short)real,(short)imag);
            cpx_buf[width-j-1] = cpx_val;
           }
         image->setRow((row),cpx_buf);
       }
    }

    delete [] cpx_buf;
    TIFFClose(tif);
    //fclose(out);

    // Flush the counter
    std::cout << std::endl;

    return image;
}

SingleBandImage<std::complex<float> > *RS2::extractSlcImage(const std::string filename)
{
    int i,j;
    //float *float_buf;
    uint32 width,height;
    FILE *out;
    SingleBandImage<std::complex<float> > *image; 

    // Get the image filename
    TIXML_STRING tmp = TinyXPath::S_xpath_string(this->document->RootElement(),
		    "/product/imageAttributes/fullResolutionImageData/text()");
    std::string input(tmp.c_str());
    std::string imageFile = this->directory + "/" + input;

    TIFFSetWarningHandler(NULL);
    TIFF *tif = TIFFOpen(imageFile.c_str(),"r");

    if (!tif) {return NULL;}

    uint16 *buf;
    uint32 row;

    //out = fopen(filename.c_str(),"wb");
    TIFFGetField(tif,TIFFTAG_IMAGEWIDTH,&width);
    TIFFGetField(tif,TIFFTAG_IMAGELENGTH,&height);
    try {
        image = new SingleBandImage<std::complex<float> >(filename.c_str(),"w",width,height);
    } catch (std::exception &e) {
	    std::cerr << e.what() << std::endl;
    }

    buf = (uint16 *)_TIFFmalloc(TIFFScanlineSize(tif));
    std::complex<float> *cpx_buf = new std::complex<float>[width];

    ez::ezETAProgressBar pg((unsigned int)height/100);
    pg.start();
    for (row = 0; row < height; row++) {
	if ((row%100) == 0) {++pg;}
        TIFFReadScanline(tif, buf, row);
	for (i=0,j=0;i<2*width;i+=2,j++) {
	    short real = (short)buf[i];
	    short imag = (short)buf[i+1];
	    std::complex<float> cpx_val((float)real,(float)imag);
	    cpx_buf[j] = cpx_val;
	    // These images are stored from late to early in azimuth
            //image->setValue(j,(height-row-1),cpx_val);
	}
	image->setRow((height-row-1),cpx_buf);
	//fwrite(float_buf,sizeof(float),2*width,out);
    }
 
    delete [] cpx_buf;
    TIFFClose(tif);
    //fclose(out);

    // Flush the counter
    std::cout << std::endl;

    return image;
}

/**
 * Convenience function for querying an XML element by an XPATH
 *
 * @param xpath the xpath for the element
 * @return the string value of that element
 */
std::string
RS2::readElement(const char *xpath)
{
    TIXML_STRING tmp = TinyXPath::S_xpath_string(this->document->RootElement(),xpath);
    std::string value(tmp.c_str());

    return value;
}

Instrument *
RS2::getInstrument()
{
    return this->instrument;
}

Scene *
RS2::getScene()
{
    return this->scene;
}

Orbit *
RS2::getOrbit()
{
    return this->orbit;
}
