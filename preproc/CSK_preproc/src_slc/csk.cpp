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
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cfloat>
#include <sys/types.h>
#include <sys/fcntl.h>
#include <sys/mman.h>
#include <boost/date_time/posix_time/posix_time_io.hpp>
#include "ezETAProgressBar.hpp"
#include "H5Cpp.h"
#include "csk.hh"

CSK::CSK(const std::string filename)
{
  /* Open the file and get the dataset */ 
  try {
    // The chunk size of most CSK datasets is 128x128x2 entries of 2 bytes each
    // Set the chunk cache size to be ~100 Mb instead of the default 1 Mb
    H5::FileAccPropList propList = H5::FileAccPropList::DEFAULT;
    int wdc = 0;
    size_t rdcc_nelems = 512;
    size_t rdcc_nbytes = 512*1024*1024;
    double w0 = 0.75;
    propList.setCache(wdc,rdcc_nelems,rdcc_nbytes,w0);

    this->file = new H5::H5File(filename, H5F_ACC_RDONLY, H5::FileCreatPropList::DEFAULT,propList);
    propList = this->file->getAccessPlist();
    propList.getCache(wdc,rdcc_nelems,rdcc_nbytes,w0);
    std::cout << "Metadata Cache " << wdc << " (elements)"
	      << " Data Chunk Cache " << rdcc_nelems << " (chunks)"
	      << " Data Chunk Cache " << rdcc_nbytes << " (bytes)"
	      << " Preemption policy " << w0 << std::endl;
  } catch (H5::FileIException error) {
    error.printError();
  }

  this->instrument = new Instrument();
  this->scene = new Scene();
  this->orbit = new Orbit();
  // Get the Product Type and acquisition mode
  this->productType = getProductType();
  this->acquisitionMode = getAcquisitionMode();
  // Now populate the platform and instrument objects
  this->populatePlatform();
  this->populateInstrument();
  this->populateScene();
  this->populateOrbit();
}

CSK::~CSK()
{
 this->file->close();
}

void CSK::populatePlatform()
{
  std::string satelliteID;

  // Extract the satellite name
  satelliteID = this->readAttribute<std::string>("/","Satellite ID");

  //this->platform->setPlanet(Planet::Earth);
  this->instrument->getPlatform()->setMission(satelliteID);
}

void CSK::populateInstrument()
{
	int antennaSide;
	std::string antennaString;
	double wavelength;
	double prf;
	double pulseLength;
	double chirpSlope;
	double antennaLength;
	double rangeSamplingFrequency;
	H5::Attribute attr;
	H5::DataType type;

	double C = 299792458.0;

	// Antenna side
	antennaString = this->readAttribute<std::string>("/","Look Side");
	if (antennaString.compare("RIGHT") != 0)
	{
		antennaSide = 1;
	}
	else
	{
		antennaSide = -1;
	}

	// Antenna wavelength
	wavelength = this->readAttribute<double>("/","Radar Wavelength");

	// PRF
	prf = this->readAttribute<double>("/S01","PRF");
	// In spotlight acquisition mode, the focused image has a different effective PRF from the 
	// nominal PRF listed in the instrument data
	if (this->acquisitionMode == SPOTLIGHT) {
		double lineTimeInterval = this->readAttribute<double>("/S01","SBI","Line Time Interval");
		prf = 1.0/lineTimeInterval;
	}

	// Pulse Length
	pulseLength = this->readAttribute<double>("/S01","Range Chirp Length");

	// Chirp Slope	
	chirpSlope = this->readAttribute<double>("/S01","Range Chirp Rate");
	
	// Antenna Length
	antennaLength = this->readAttribute<double>("/","Antenna Length");

	// Range Sampling Frequency
	rangeSamplingFrequency = this->readAttribute<double>("/S01","Sampling Rate");

	// Range Pixel Size
	double rangePixelSize = C/(2.0*rangeSamplingFrequency);

        if (this->productType == SCS_B) {
		// Range Pixel Size
		rangePixelSize = this->readAttribute<double>("/S01","SBI","Column Spacing");

		// Azimuth Pixel Size
		double azimuthPixelSize = this->readAttribute<double>("/S01","SBI","Line Spacing");
		this->instrument->setAzimuthPixelSize(azimuthPixelSize);
	}

	this->instrument->setAntennaSide(antennaSide);
	this->instrument->setWavelength(wavelength);
	this->instrument->setPulseRepetitionFrequency(prf);
	this->instrument->setPulseLength(pulseLength);
	this->instrument->setChirpSlope(chirpSlope);
	this->instrument->setAntennaLength(antennaLength);
	this->instrument->setRangeSamplingFrequency(rangeSamplingFrequency);
	this->instrument->setRangePixelSize(rangePixelSize);
}

void CSK::populateScene()
{
    int orbitNumber;
    double rangeFirstTime;
    double slantRange;
    double satelliteHeight;
    std::string processingFacility;
    std::string processingLevel;
    std::string polarization;
    std::string orbitDirection;
    
    H5::Attribute attr;
    H5::DataType type;
    H5::DataSpace space;
   
    // Range First Times
    if (this->productType == RAW_B) {
	rangeFirstTime = this->readAttribute<double>("/S01","B001","Range First Times");
	processingLevel = "RAW";
    } else if (this->productType == SCS_B) {
	rangeFirstTime = this->readAttribute<double>("/S01","SBI","Zero Doppler Range First Time");
	processingLevel = "SLC";
    } else {
        std::cerr << "Unable to calculate slant range" << std::endl;
	return;
    }
    
    // Satellite Height
    satelliteHeight = this->readAttribute<double>("/","Satellite Height");
    
    // Processing Facility
    processingFacility = this->readAttribute<std::string>("/","Processing Centre");

    // Processing Level
    //processingLevel = this->readAttribute<std::string>("/","L0 Software Version");

    // Polarization
    polarization = this->readAttribute<std::string>("/S01","Polarisation");

    // Orbit Number
    orbitNumber = this->readAttribute<int>("/","Orbit Number");

    // Orbit Direction
    orbitDirection = this->readAttribute<std::string>("/","Orbit Direction");

    boost::posix_time::ptime sensingStart;
    boost::posix_time::ptime sensingStop;

    if (this->productType == SCS_B)
    {
        // Reference UTC time
	std::string referenceString = this->readAttribute<std::string>("/","Reference UTC");
        boost::posix_time::ptime referenceUTC (boost::posix_time::time_from_string(referenceString));

        // Scene starting time for focused data is the zero doppler time
	double zeroDopplerStartTime = this->readAttribute<double>("/S01","SBI","Zero Doppler Azimuth First Time");

        // Scene stopping time for focused data is the zero doppler time
	double zeroDopplerStopTime = this->readAttribute<double>("/S01","SBI","Zero Doppler Azimuth Last Time");

	boost::posix_time::time_duration dt = boost::posix_time::microseconds(zeroDopplerStartTime*1e6);
	sensingStart = referenceUTC + dt;
	dt = boost::posix_time::microseconds(zeroDopplerStopTime*1e6);
	sensingStop  = referenceUTC + dt;
    }
    else 
    {
        // Scene starting time
	std::string sensingStartString = this->readAttribute<std::string>("/","Scene Sensing Start UTC");
	boost::posix_time::ptime tmpStart (boost::posix_time::time_from_string(sensingStartString));
	sensingStart = tmpStart;

        // Scene stopping time
	std::string sensingStopString = this->readAttribute<std::string>("/","Scene Sensing Stop UTC");
	boost::posix_time::ptime tmpStop (boost::posix_time::time_from_string(sensingStopString));
	sensingStop = tmpStop;
    }

    double C = 299792458.0;
    slantRange = rangeFirstTime*C/2.0;
    
    this->scene->setIBias(127.5);
    this->scene->setQBias(127.5);
    this->scene->setStartingRange(slantRange);
    this->scene->setSatelliteHeight(satelliteHeight);
    this->scene->setProcessingFacility(processingFacility);
    this->scene->setProcessingLevel(processingLevel);
    this->scene->setPolarization(polarization);
    this->scene->setOrbitNumber(orbitNumber);
    this->scene->setPassDirection(orbitDirection);
    this->scene->setSensingStart(sensingStart);
    this->scene->setSensingStop(sensingStop);
    
}

void CSK::populateOrbit()
{
    int ndims;
    double *t,*pos,*vel;
    hsize_t dims[1];
    std::string referenceTime;
    H5::Attribute attr;
    H5::DataType type;
    H5::DataSpace space;

    // Reference Time
    referenceTime = this->readAttribute<std::string>("/","Reference UTC");
    boost::posix_time::ptime t0(boost::posix_time::time_from_string(referenceTime));

    // Get State Vector times
    attr = this->file->openGroup("/").openAttribute("State Vectors Times");
    space = attr.getSpace();
    ndims = space.getSimpleExtentDims(dims,NULL);
    t = new double[dims[0]];
    pos = new double[dims[0]*3];
    vel = new double[dims[0]*3];
    attr.read(H5::PredType::NATIVE_DOUBLE,t);

    // Read Position and velocity 
    // (this assumes that there are the same number of state vector times as positions and velocities)
    attr = this->file->openGroup("/").openAttribute("ECEF Satellite Position");
    attr.read(H5::PredType::NATIVE_DOUBLE,pos);
    attr = this->file->openGroup("/").openAttribute("ECEF Satellite Velocity");
    attr.read(H5::PredType::NATIVE_DOUBLE,vel);

    for (int i=0;i<dims[0];i++)
    {
        StateVector sv = StateVector();
	double thisPos[3];
	double thisVel[3];
    	for (int j=0;j<3;j++)
	{
	    thisPos[j] = pos[j+i*3];
	    thisVel[j] = vel[j+i*3];
	}
	sv.setTime(t0 + boost::posix_time::seconds(t[i]));
	sv.setPosition(thisPos);
	sv.setVelocity(thisVel);
	this->orbit->addStateVector(sv);
    }

    delete [] t;
    delete [] pos;
    delete [] vel;
}

void CSK::readLUT()
{
  int i, ndims;
  hsize_t dims[1];

  H5::Attribute attr = this->file->openGroup("/").openAttribute("Analog Signal Reconstruction Levels");
  H5::DataSpace space = attr.getSpace();
  ndims = space.getSimpleExtentDims(dims,NULL);

  this->lut = new double[dims[0]];
  attr.read(H5::PredType::NATIVE_DOUBLE,this->lut);

  for(i=0;i<(int)(dims[0]);i++)
    {
      if(std::isnan(this->lut[i]))
	{
	  this->lut[i] = 0.0;
        }
    }
}

int CSK::getProductType()
{
  int productTypeCode;

  std::string productType = this->readAttribute<std::string>("/","Product Type");

  if (productType.compare("RAW_B") == 0) {productTypeCode = RAW_B;}
  else if (productType.compare("SCS_B") == 0) {productTypeCode = SCS_B;}
  else {std::cerr << "Unknown product type: " << productType << std::endl; productTypeCode = UNK_B;}

  return productTypeCode;
}

int CSK::getAcquisitionMode()
{
   int acquisitionModeCode;

   std::string acquisitionMode = this->readAttribute<std::string>("/","Acquisition Mode");
   
   if (acquisitionMode.compare("HIMAGE") == 0) {acquisitionModeCode = HIMAGE;}
   else if (acquisitionMode.compare("ENHANCED SPOTLIGHT") == 0) {acquisitionModeCode = SPOTLIGHT;}
   else {std::cerr << "Unknown acquisition mode: " << acquisitionMode << std::endl; acquisitionModeCode = UNKNOWN;}
   this->scene->setBeam(acquisitionMode);

   return acquisitionModeCode;
}

SingleBandImage<std::complex<unsigned char> > *CSK::extractRawImage(const std::string outFile)
{
  unsigned char *IQ;
  int i,j,k,ndims;
  hsize_t dims[3],offset[3],count[3];
  hsize_t dimsm[1], offset_out[1],count_out[1];

  if (this->productType != RAW_B)
  {
      std::cerr << "Image does not appear to be of type RAW_B" << std::endl;
      return NULL;
  }

  readLUT();
  double max = calculateLUTRange();
  max = hypot(max,max);

  H5::DataSet dataset = this->file->openDataSet("/S01/B001");
  H5::DataType type = dataset.getDataType();

  H5::DataSpace dataspace = dataset.getSpace();
  ndims = dataspace.getSimpleExtentDims(dims,NULL);
  
  int height = dims[0];
  int width = dims[1];
  // Create image object
  SingleBandImage<std::complex<unsigned char> > *image = 
	  new SingleBandImage<std::complex<unsigned char> >(outFile.c_str(),"w",width,height);

  dimsm[0] = 2*dims[1];
  H5::DataSpace memspace = H5::DataSpace(1,dimsm,NULL);
  
  offset_out[0] = 0;
  count_out[0] = 2*dims[1];

  IQ = new unsigned char [2*width];
  std::complex<unsigned char> *cpx_buf = new std::complex<unsigned char>[width];
  ez::ezETAProgressBar pg((unsigned int)height/100);
  pg.start();
  for(k=0;k<(int)(height);k++)
    {
      if ((k%100) == 0)
	{
	  ++pg;
	}

      offset[0] = (hsize_t)k;
      offset[1] = 0;
      offset[2] = 0;
      count[0] = 1;
      count[1] = (hsize_t)dims[1];
      count[2] = 2;
      dataspace.selectHyperslab(H5S_SELECT_SET,count,offset,NULL,NULL);
      memspace.selectHyperslab(H5S_SELECT_SET,count_out,offset_out,NULL,NULL);
      dataset.read(IQ,type,memspace,dataspace,H5P_DEFAULT);
      
      for(i=0,j=0;i<(int)(2*dims[1]);i+=2,j++)
	{
	  double real = this->lut[(int)IQ[i]];
	  double imag = this->lut[(int)IQ[i+1]];
          unsigned char I = (unsigned char)(127.0*real/max+127.0);
          unsigned char Q = (unsigned char)(127.0*imag/max+127.0);
	  cpx_buf[j] = std::complex<unsigned char>(I,Q);
	}
      image->setRow(k,cpx_buf);
    }
  delete [] IQ;
  delete [] cpx_buf;

  // Extra line to flush the counter
  std::cout << std::endl;

  return image;
}

SingleBandImage<std::complex<short> > *CSK::extractSlc_short_Image(const std::string outFile)
{
  short *IQ;
  int i,j,k,ndims;
  hsize_t dims[3],offset[3],count[3];
  hsize_t dimsm[1], offset_out[1],count_out[1];

  if (this->productType != SCS_B)
  {
      std::cerr << "Image does not appear to be of type SCS_B" << std::endl;
      return NULL;
  }

  H5::DataSet dataset = this->file->openDataSet("/S01/SBI");
  H5::DataType type = dataset.getDataType();

  H5::DataSpace dataspace = dataset.getSpace();
  ndims = dataspace.getSimpleExtentDims(dims,NULL);

  int height = dims[0];
  int width = dims[1];
  
  // Create image object
  SingleBandImage<std::complex<short> > *image = 
	  new SingleBandImage<std::complex<short> >(outFile.c_str(),"w",width,height);

  dimsm[0] = 2*dims[1];
  H5::DataSpace memspace = H5::DataSpace(1,dimsm,NULL);
  
  offset_out[0] = 0;
  count_out[0] = 2*dims[1];

  IQ = new short [2*dims[1]];
  std::complex<short> *cpx_buf = new std::complex<short>[dims[1]];
  ez::ezETAProgressBar pg((unsigned int)dims[0]/100);
  pg.start();
  for(k=0;k<(int)(dims[0]);k++)
    {
      if ((k%100) == 0)
	{
		++pg;
	}

      offset[0] = (hsize_t)k;
      offset[1] = 0;
      offset[2] = 0;
      count[0] = 1;
      count[1] = (hsize_t)dims[1];
      count[2] = 2;
      dataspace.selectHyperslab(H5S_SELECT_SET,count,offset,NULL,NULL);
      memspace.selectHyperslab(H5S_SELECT_SET,count_out,offset_out,NULL,NULL);
      dataset.read(IQ,type,memspace,dataspace,H5P_DEFAULT);
 
      for(i=0,j=0;i<(int)(2*dims[1]);i+=2,j++)
	{
	  cpx_buf[j] = std::complex<short>((short)IQ[i],(short)IQ[i+1]);
	}
      image->setRow(k,cpx_buf);
    }
  delete [] IQ;
  delete [] cpx_buf;

  // Extra line to flush the counter
  std::cout << std::endl;

  return image;
}

SingleBandImage<std::complex<float> > *CSK::extractSlcImage(const std::string outFile)
{
  short *IQ;
  int i,j,k,ndims;
  hsize_t dims[3],offset[3],count[3];
  hsize_t dimsm[1], offset_out[1],count_out[1];

  if (this->productType != SCS_B)
  {
      std::cerr << "Image does not appear to be of type SCS_B" << std::endl;
      return NULL;
  }

  H5::DataSet dataset = this->file->openDataSet("/S01/SBI");
  H5::DataType type = dataset.getDataType();

  H5::DataSpace dataspace = dataset.getSpace();
  ndims = dataspace.getSimpleExtentDims(dims,NULL);

  int height = dims[0];
  int width = dims[1];
  
  // Create image object
  SingleBandImage<std::complex<float> > *image = 
	  new SingleBandImage<std::complex<float> >(outFile.c_str(),"w",width,height);

  dimsm[0] = 2*dims[1];
  H5::DataSpace memspace = H5::DataSpace(1,dimsm,NULL);
  
  offset_out[0] = 0;
  count_out[0] = 2*dims[1];

  IQ = new short [2*dims[1]];
  std::complex<float> *cpx_buf = new std::complex<float>[dims[1]];
  ez::ezETAProgressBar pg((unsigned int)dims[0]/100);
  pg.start();
  for(k=0;k<(int)(dims[0]);k++)
    {
      if ((k%100) == 0)
	{
		++pg;
	}

      offset[0] = (hsize_t)k;
      offset[1] = 0;
      offset[2] = 0;
      count[0] = 1;
      count[1] = (hsize_t)dims[1];
      count[2] = 2;
      dataspace.selectHyperslab(H5S_SELECT_SET,count,offset,NULL,NULL);
      memspace.selectHyperslab(H5S_SELECT_SET,count_out,offset_out,NULL,NULL);
      dataset.read(IQ,type,memspace,dataspace,H5P_DEFAULT);
 
      for(i=0,j=0;i<(int)(2*dims[1]);i+=2,j++)
	{
	  cpx_buf[j] = std::complex<float>((float)IQ[i],(float)IQ[i+1]);
	}
      image->setRow(k,cpx_buf);
    }
  delete [] IQ;
  delete [] cpx_buf;

  // Extra line to flush the counter
  std::cout << std::endl;

  return image;
}

double
CSK::calculateLUTRange() 
{
  int i;
  double min,max;
  
  max = DBL_MIN;
  min = DBL_MAX;
  for(i=0;i<256;i++)
    {
      if (this->lut[i] > max)
	{
	  max = this->lut[i];
	}
      if (this->lut[i] < min)
	{
	  min = this->lut[i];
	}
    }

  return max;
}

Instrument *
CSK::getInstrument()
{
	return this->instrument;
}

Scene *
CSK::getScene()
{
    return this->scene;
}

Orbit *
CSK::getOrbit()
{
	return this->orbit;
}

template <>
std::string
CSK::readAttribute<std::string>(const char *group, const char *attribute)
{
	std::string value;
        H5::Attribute attr;
        H5::DataType type;

	attr = this->file->openGroup(group).openAttribute(attribute);
	type = attr.getDataType();
	H5std_string stringbuffer("");
	attr.read(type,stringbuffer);
	value.assign(stringbuffer);

	return value;
}

template <>
std::string
CSK::readAttribute<std::string>(const char *group, const char *dataset, const char *attribute)
{
	std::string value;
        H5::Attribute attr;
        H5::DataType type;

	attr = this->file->openGroup(group).openDataSet(dataset).openAttribute(attribute);
	type = attr.getDataType();
	H5std_string stringbuffer("");
	attr.read(type,stringbuffer);
	value.assign(stringbuffer);

	return value;
}
