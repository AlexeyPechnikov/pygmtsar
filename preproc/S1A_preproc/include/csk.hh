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
#ifndef CSK_HH
#define CSK_HH 1

#include <complex>
#include "H5Cpp.h"
#include "metadata/Instrument.hh"
#include "metadata/Scene.hh"
#include "metadata/StateVector.hh"
#include "metadata/Orbit.hh"
#include "image/SingleBandImage.hh"



/**
 * A class for Cosmo-Skymed data
 */
class CSK
{
private:
  H5::H5File *file;
  double *lut;         //< Look up table for sensor values
  int productType;     
  int acquisitionMode; 
  Instrument *instrument;
  Scene *scene;
  Orbit *orbit;
  void populatePlatform();
  void populateInstrument();
  void populateScene();
  void populateOrbit();
  template <typename TType> TType readAttribute(const char *group,const char *attribute);
  template <typename TType> TType readAttribute(const char *group,const char *dataset,const char *attribute);
public:
  CSK(const std::string filename);
  ~CSK();
  /**
   * Get the data product type
   *
   * @return integer value representing the product type from the enum TProductType
   */
  int getProductType();
  /**
   * Get the data acquisition mode
   *
   * @return integer value representing the acquisition mode from the enum TAcquisitionMode
   */
  int getAcquisitionMode();
  /**
   * Read the image data lookup table.  Values are stored in *lut.
   */
  void readLUT();
  /**
   * Unpack the image data for a RAW_B product type
   *
   * @return SingleBandImage<std::complex<unsigned char> > object
   */
  SingleBandImage<std::complex<unsigned char> > *extractRawImage(const std::string outFile);
  /**
   * Unpack the image data for an SCS_B product type
   *
   * @return SingleBandImage<std::complex<short> > object
   */
  SingleBandImage<std::complex<short> > *extractSlc_short_Image(const std::string outFile);
  /**
   * Unpack the image data for an SCS_B product type
   *
   * @return SingleBandImage<std::complex<float> > object
   */
  SingleBandImage<std::complex<float> > *extractSlcImage(const std::string outFile);
  /**
   * Calculate the maximum value in the LUT
   *
   * @return the maximum value in the LUT
   */
  double calculateLUTRange();
  Instrument *getInstrument();
  Scene *getScene();
  Orbit *getOrbit();
  enum TProductType {RAW_B = 1, SCS_B = 2, UNK_B = -1};
  enum TAcquisitionMode {HIMAGE = 1, SPOTLIGHT = 2, UNKNOWN = -1};
};

template <typename TType>
TType
CSK::readAttribute(const char *group, const char *attribute)
{
	TType value;
        H5::Attribute attr;
        H5::DataType type;

	attr = this->file->openGroup(group).openAttribute(attribute);
	type = attr.getDataType();
	attr.read(type,&value);

	return value;
}

template <typename TType>
TType
CSK::readAttribute(const char *group, const char *dataset, const char *attribute)
{
	TType value;
        H5::Attribute attr;
        H5::DataType type;

	attr = this->file->openGroup(group).openDataSet(dataset).openAttribute(attribute);
	type = attr.getDataType();
	attr.read(type,&value);

	return value;
}

template <> std::string CSK::readAttribute(const char *group, const char *attibute);
template <> std::string CSK::readAttribute(const char *group, const char *dataset, const char *attibute);

#endif
