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
 * converted from rs2.cpp for purpose of processing sentinel1-a data
 *
 */

#include<stdio.h>
#include<stdlib.h>
#include <iostream>
#include <string>
#include <deque>
#include <sstream>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>
#include <libgen.h>
#include "ezETAProgressBar.hpp"
#include "xpath_static.h"
#include "tiffio.h"
#include "s1a.hh"

int ki[100000],ko[100000];
S1A::S1A(const std::string filename)
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
    this->populateBurst();
}

S1A::~S1A()
{
    delete this->document;
}

void S1A::populatePlatform()
{
  TIXML_STRING tmp = TinyXPath::S_xpath_string(this->document->RootElement(),
                  "/product/adsHeader/missionId/text()");
  std::string satelliteID(tmp.c_str());

  this->instrument->getPlatform()->setMission(satelliteID);
}

void S1A::populateInstrument()
{
    double c = 299792458.0;
    double frequency = TinyXPath::d_xpath_double(this->document->RootElement(), 
                    "/product/generalAnnotation/productInformation/radarFrequency/text()");
    double wavelength = c/frequency;
    double prf = TinyXPath::d_xpath_double(this->document->RootElement(),
		    "/product/imageAnnotation/imageInformation/azimuthFrequency/text()");
    double pulseLength = TinyXPath::d_xpath_double(this->document->RootElement(),
		    "/product/generalAnnotation/downlinkInformationList/downlinkInformation/downlinkValues/txPulseLength/text()");
    double pulseBandwidth = TinyXPath::d_xpath_double(this->document->RootElement(),
		    "/product/imageAnnotation/processingInformation/swathProcParamsList/swathProcParams/rangeProcessing/lookBandwidth/text()");
    double chirpSlope = pulseBandwidth/pulseLength;
    double rangePixelSize = TinyXPath::d_xpath_double(this->document->RootElement(),
		    "/product/imageAnnotation/imageInformation/rangePixelSpacing/text()");
    double rangeSamplingFrequency = TinyXPath::d_xpath_double(this->document->RootElement(),
            "/product/generalAnnotation/productInformation/rangeSamplingRate/text()");
    
    this->instrument->setAntennaSide(-1);
    this->instrument->setWavelength(wavelength);
    this->instrument->setPulseRepetitionFrequency(prf);
    this->instrument->setPulseLength(pulseLength);
    this->instrument->setChirpSlope(chirpSlope);
    this->instrument->setAntennaLength(12.3);//?
    this->instrument->setRangeSamplingFrequency(rangeSamplingFrequency);
    this->instrument->setRangePixelSize(rangePixelSize);
}

void S1A::populateScene()
{
    double c = 299792458.0;
    double height = 693000;
    double firstSampleTime = TinyXPath::d_xpath_double(this->document->RootElement(),
		    "/product/imageAnnotation/imageInformation/slantRangeTime/text()");
    double startingRange = firstSampleTime*c/2.0;
    std::string facility = this->readElement("/product/imageGenerationParameters/generalProcessingInformation/processingFacility/text()"); //?
    std::string version = "SLC";
    std::string passDirection = this->readElement("/product/generalAnnotation/productInformation/pass/text()");

    // get orbit direction
    std::string DIR = passDirection.substr(0,1);
    std::string ASC ("A");
    std::string dataStartString = this->readElement("/product/adsHeader/startTime/text()");

    // Set up the date time formatter
    boost::posix_time::time_input_facet *facet = new boost::posix_time::time_input_facet("%Y-%m-%dT%H:%M:%S%F");
    std::stringstream sensingStartString(dataStartString.c_str());
    sensingStartString.imbue(std::locale(sensingStartString.getloc(),facet));
    boost::posix_time::ptime sensingStart;
    sensingStartString >> sensingStart;

    double prf = TinyXPath::d_xpath_double(this->document->RootElement(),
		    "/product/imageAnnotation/imageInformation/azimuthFrequency/text()");
    int lines = TinyXPath::i_xpath_int(this->document->RootElement(),
		    "/product/imageAnnotation/imageInformation/numberOfLines/text()");
    double msec = 1e6* double(lines)/prf;
    boost::posix_time::ptime sensingStop = sensingStart + boost::posix_time::microseconds(msec);
    std::string polarization = this->readElement("/product/adsHeader/polarisation/text()");
    std::string beam = this->readElement("/product/adsHeader/mode/text()"); //?

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

void S1A::populateBurst()
{

//don't set the hash table if this not is an IW beam
    if(scene->getBeam() != "IW") {
      //fprintf(stdout," must not be an IW ");
      for (unsigned int ii=0;ii<100000;ii++) {
        ki[ii]=0;
        ko[ii]=0;
        //fprintf(stdout," %d %d\n",ki[ii],ko[ii]);
      }
    }

//set the hash table
    else {
    TinyXPath::xpath_processor xp_proc(this->document->RootElement(),
		    "/product/swathTiming/burstList/burst");
    unsigned int count = xp_proc.u_compute_xpath_node_set();
    int lpb = TinyXPath::i_xpath_int(this->document->RootElement(),
		"/product/swathTiming/linesPerBurst/text()");

    double prf = TinyXPath::d_xpath_double(this->document->RootElement(),
		    "/product/imageAnnotation/imageInformation/azimuthFrequency/text()");
    double time = 0.,time0 = -1.;
    int k=0,ktm=0;
    int k_start = -1, k_tot = -1;
    char *cflag, *cflag0;
    int flag[lpb];

    cflag = (char*)malloc(sizeof(char) * 10*( lpb + 1 ) );

    // loop over the bursts
    for (unsigned int i=0;i<count;i++) {

	TiXmlNode *stateVector = xp_proc.XNp_get_xpath_node(i);
	double btime  = TinyXPath::d_xpath_double(stateVector,"//azimuthAnxTime/text()");

	TIXML_STRING tmp = TinyXPath::S_xpath_string(stateVector,"//firstValidSample/text()");

        cflag0 = cflag;
	strcpy(cflag0,tmp.c_str());

    // loop over the lines in each burst and extract the flag information
        for (unsigned int j=0;j<lpb;j++) {
          flag[j] = (int)strtol(cflag0,&cflag0,10);
          time = btime + double(j)/prf;

    // don't use the flagged data
          ktm = -1;
    	  if(flag[j] != -1) {
            if(time0 < 0.) {
              time0 = time;
              k_start = k;
            }

    // map the echo sequence number into the original line number.  add a bit so it won't round down.
            ktm=(int)(((time-time0)*prf)+0.1);
          }
	  ki[k]=k;
          ko[k]=ktm;
          if(k_tot < ktm) k_tot = ktm;
          //fprintf(stdout," %d %d %d\n",ki[k],j,ko[k]);
          k++;
        }
        //fprintf(stdout," %d %d %lf \n",k_start,k_tot,prf);

   // now reset the start time
    }
	boost::posix_time::ptime sensingStart;
	sensingStart = this->scene->getSensingStart();
        double msec = 1e6* double(k_start)/prf;
        sensingStart = sensingStart + boost::posix_time::microseconds(msec);
        this->scene->setSensingStart(sensingStart);
    free(cflag);
  }
}

void S1A::populateOrbit()
{
    double pos[3];
    double vel[3];

    // Get the number of state vectors
    TinyXPath::xpath_processor xp_proc(this->document->RootElement(),
		    "/product/generalAnnotation/orbitList/orbit");
    unsigned int count = xp_proc.u_compute_xpath_node_set();
    // fprintf(stderr,"count = %d \n",count);
     
    // Set up the date time formatter
    boost::posix_time::time_input_facet *facet = new boost::posix_time::time_input_facet("%Y-%m-%dT%H:%M:%S%F");
    std::stringstream epochString;
    epochString.imbue(std::locale(epochString.getloc(),facet));
    
    for (unsigned int i=0;i<count;i++) {
        // Parse XML
	TiXmlNode *stateVector = xp_proc.XNp_get_xpath_node(i);
	TIXML_STRING tmp = TinyXPath::S_xpath_string(stateVector,"//time/text()");
	pos[0] = TinyXPath::d_xpath_double(stateVector,"//position/x/text()");
	pos[1] = TinyXPath::d_xpath_double(stateVector,"//position/y/text()");
	pos[2] = TinyXPath::d_xpath_double(stateVector,"//position/z/text()");
	vel[0] = TinyXPath::d_xpath_double(stateVector,"//velocity/x/text()");
	vel[1] = TinyXPath::d_xpath_double(stateVector,"//velocity/y/text()");
	vel[2] = TinyXPath::d_xpath_double(stateVector,"//velocity/z/text()");

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

SingleBandImage<std::complex<short> > *S1A::extractSlc_short_Image(const std::string filename,const std::string inputfile)
{
    int i,j,iheight=-1;
    uint32 width,height0,height,widthi;
    //FILE *out;
    SingleBandImage<std::complex<short> > *image; 

    // Get the image filename
    //std::cout << filename.c_str() << std::endl;
    std::string tmp = inputfile.c_str();
    unsigned pos = tmp.find(".xml");
    std::string input = tmp.substr(0,pos) + ".tiff";
    std::string imageFile = this->directory + "/" + input;
    
    std::cout << imageFile << std::endl;
    
    TIFFSetWarningHandler(NULL);
    TIFF *tif = TIFFOpen(imageFile.c_str(),"r");
    if (!tif) {return NULL;}
    //std::cout << "sp1" << std::endl;
    uint16 *buf;
    uint32 row,row_out;

    //out = fopen(filename.c_str(),"wb");
    TIFFGetField(tif,TIFFTAG_IMAGEWIDTH,&widthi);
    TIFFGetField(tif,TIFFTAG_IMAGELENGTH,&height0);
    
    // make sure the width is divisible by 4
    width = widthi - widthi%4;
    if(scene->getBeam() != "IW") {
      height = height0;
    } 
    else {
    // examine the ko[] array to reset the height
      for (row = 0; row < height0; row++) {
        if(iheight < ko[row]) iheight = ko[row];
      }
      height = iheight +1;
    }
    try {
        image = new SingleBandImage<std::complex<short> >(filename.c_str(),"w",width,height);
    } catch (std::exception &e) {
	    std::cerr << e.what() << std::endl;
    }

    buf = (uint16 *)_TIFFmalloc(TIFFScanlineSize(tif));
    std::complex<short> *cpx_buf = new std::complex<short>[width];

    ez::ezETAProgressBar pg((unsigned int)height/100);
    pg.start();

    // extract image
    // first case is standard swath mode data 
    if(scene->getBeam() != "IW") {
      for (row = 0; row < height0; row++) {
        if ((row%100) == 0) {++pg;}
        TIFFReadScanline(tif, buf, row);
         for (i=0,j=0;i<2*width;i+=2,j++) {
            short real = (short)buf[i];
            short imag = (short)buf[i+1];
            std::complex<short> cpx_val((short)real,(short)imag);
            cpx_buf[j] = cpx_val;
	   }
         image->setRow((row),cpx_buf);
       }
     }
   //  must be burst mode data so only write non zero lines
     else {
      for (row = 0; row < height0; row++) {
        if ((row%100) == 0) {++pg;}
        TIFFReadScanline(tif, buf, row);
        if(ko[row] > 0) {
         for (i=0,j=0;i<2*width;i+=2,j++) {
            short real = (short)buf[i];
            short imag = (short)buf[i+1];
            std::complex<short> cpx_val((short)real,(short)imag);
            cpx_buf[j] = cpx_val;
	   }
           row_out = ko[row];
         image->setRow(row_out,cpx_buf);
        }
       }
     }

    delete [] cpx_buf;
    TIFFClose(tif);
    //fclose(out);

    // Flush the counter
    std::cout << std::endl;

    return image;
}

SingleBandImage<std::complex<float> > *S1A::extractSlcImage(const std::string filename)
{
    int i,j;
    //float *float_buf;
    uint32 width,height;
    //FILE *out;
    SingleBandImage<std::complex<float> > *image; 

    // Get the image filename
    std::string tmp = filename.c_str();
    unsigned pos = tmp.find(".xml");
    std::string input = tmp.substr(0,pos-1) + ".tiff";
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
S1A::readElement(const char *xpath)
{
    TIXML_STRING tmp = TinyXPath::S_xpath_string(this->document->RootElement(),xpath);
    std::string value(tmp.c_str());

    return value;
}

Instrument *
S1A::getInstrument()
{
    return this->instrument;
}

Scene *
S1A::getScene()
{
    return this->scene;
}

Orbit *
S1A::getOrbit()
{
    return this->orbit;
}
