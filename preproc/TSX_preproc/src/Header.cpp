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
#include <fstream>
#include "byteswap.h"
#include "Header.hh"

Header::Header(bool isBigEndian)
{
  this->isBigEndian = isBigEndian;
}

void
Header::parse(std::istream &fin)
{
 fin.read((char *)(&this->bytesInBurst),sizeof(unsigned int));
 fin.read((char *)(&this->rangeSampleRelativeIndex),sizeof(int));
 fin.read((char *)(&this->rangeSamples),sizeof(int));
 fin.read((char *)(&this->azimuthSamples),sizeof(int));
 fin.read((char *)(&this->burstIndex),sizeof(int));
 fin.read((char *)(&this->rangelineTotalNumberOfBytes),sizeof(int));
 fin.read((char *)(&this->totalNumberOfLines),sizeof(int));
 fin.read((char *)(&this->format),4*sizeof(char));
 fin.read((char *)(&this->version),sizeof(int));
 fin.read((char *)(&this->oversamplingFactor),sizeof(int));
 fin.read((char *)(&this->inverseSPECANScalingRate),sizeof(double));

 if (!this->isBigEndian)
 {
  // Byte swap all of these
  bytesInBurst = bswap_32(bytesInBurst);
  rangeSampleRelativeIndex = bswap_32(rangeSampleRelativeIndex);
  rangeSamples = bswap_32(rangeSamples);
  azimuthSamples = bswap_32(azimuthSamples);
  burstIndex = bswap_32(burstIndex);
  rangelineTotalNumberOfBytes = bswap_32(rangelineTotalNumberOfBytes);
  totalNumberOfLines = bswap_32(totalNumberOfLines);
  version = bswap_32(version);
  oversamplingFactor = bswap_32(oversamplingFactor);
  inverseSPECANScalingRate = bswap_64(inverseSPECANScalingRate);
 }

 std::string formatCheck("CSAR");
 if (formatCheck.compare(format) != 0) {throw "Not a valid COSAR file";}
 // Skip to the end of the header line
 fin.seekg((rangelineTotalNumberOfBytes-48),std::ios_base::cur);
}

void
Header::print()
{
 std::cout << "Bytes In Burst " << bytesInBurst << std::endl;
 std::cout << "Range Sample Relative Index " << rangeSampleRelativeIndex << std::endl;
 std::cout << "Range Samples " << rangeSamples << std::endl;
 std::cout << "Azimuth Samples " << azimuthSamples << std::endl;
 std::cout << "Burst Index " << burstIndex << std::endl;
 std::cout << "Rangeline Total Number of Bytes " << rangelineTotalNumberOfBytes << std::endl;
 std::cout << "Total Number of Lines " << totalNumberOfLines << std::endl;
 std::cout << "Format " << format << std::endl;
 std::cout << "Version " << version << std::endl;
 std::cout << "Oversampling Factor " << oversamplingFactor << std::endl;
 // Then skip ahead 'Rangeline Total Number of Bytes - 48' to get to the end of the header
}

int
Header::getRangeSamples()
{
 return this->rangeSamples;
}

int
Header::getAzimuthSamples()
{
 return this->azimuthSamples;
}

int
Header::getRangelineTotalNumberOfBytes()
{
 return this->rangelineTotalNumberOfBytes;
}

int
Header::getTotalNumberOfLines()
{
 return this->totalNumberOfLines;
}

int
Header::getBytesInBurst()
{
 return this->bytesInBurst;
}
