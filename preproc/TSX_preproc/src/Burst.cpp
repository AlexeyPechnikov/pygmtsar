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
#include "Burst.hh"
#include "ezETAProgressBar.hpp"

Burst::Burst(int rangeSamples,int azimuthSamples,bool isBigEndian) 
{
 this->isBigEndian = isBigEndian;
 this->azimuthSamples = azimuthSamples;
 this->rangeSamples = rangeSamples;
 this->asri = new int[this->rangeSamples];
 this->asfv = new int[this->rangeSamples];
 this->aslv = new int[this->rangeSamples];
}

Burst::~Burst() 
{
 delete [] this->asri;
 delete [] this->asfv;
 delete [] this->aslv;
}

void
Burst::parse(std::istream &fin,std::ostream &fout)
{
 this->parseAzimuthHeader(fin);

 short *data = new short[2*this->rangeSamples];
 short *shortData = new short[2*this->rangeSamples];

 ez::ezETAProgressBar pg((unsigned int)this->azimuthSamples/1000);
 pg.start();
 for(int i=0;i<this->azimuthSamples;i++)
   {
    if ((i % 1000) == 0)
      {
	      ++pg;
      }
    this->parseRangeLine(fin,fout,i,data,shortData);
   }
  // Flush progress bar
  std::cout << std::endl;

  delete [] data;
  delete [] shortData;
}

void
Burst::parseAzimuthHeader(std::istream &fin)
{
 // For each of the three azimuth header lines, skip the first two 4-byte samples
 // Read 'Range Samples' number of 4 byte integers 
 fin.seekg(8, std::ios_base::cur);
 fin.read((char *)(this->asri),this->rangeSamples*sizeof(int));
 // again
 fin.seekg(8, std::ios_base::cur);
 fin.read((char *)(this->asfv),this->rangeSamples*sizeof(int));
 // and again
 fin.seekg(8, std::ios_base::cur);
 fin.read((char *)(this->aslv),this->rangeSamples*sizeof(int));

 if (!this->isBigEndian)
   {
     // Byte swap
     for(int i=0;i<this->rangeSamples;i++)
      {
       this->asri[i] = bswap_32(this->asri[i]);
       this->asfv[i] = bswap_32(this->asfv[i])-1;
       this->aslv[i] = bswap_32(this->aslv[i])-1;
      }
   }
}
void
Burst::parseRangeLine(std::istream &fin,std::ostream &fout,int lineNumber,short *data, short *shortData)
{
 //short *data;
 int rsfv,rslv;
 int asfv,aslv;
 //float *floatData;

 /*data = new short[2*this->rangeSamples];
 shortData = new short[2*this->rangeSamples];*/

 // Read line header
 fin.read((char*)(&rsfv),sizeof(int));
 fin.read((char*)(&rslv),sizeof(int));
 if (!this->isBigEndian)
   {
    // Byte swap
    rsfv = bswap_32(rsfv)-1;
    rslv = bswap_32(rslv)-1;
   }
 // Read data
 fin.read((char*)(data),2*this->rangeSamples*sizeof(short));
 // Byte swap data and mask out invalid points
 for(int rangeBin=0,j=0;rangeBin<this->rangeSamples;rangeBin++,j+=2)
   {
     asfv = this->asfv[rangeBin];
     aslv = this->aslv[rangeBin];
     // gdal_translate
     if ((lineNumber < asfv) || (lineNumber > aslv) || (rangeBin < rsfv) || (rangeBin >= rslv))
      {
       shortData[j] = 0;
       shortData[j+1] = 0;
      }
    else
      {
       if (!this->isBigEndian)
	 {
          data[j] = bswap_16(data[j]);
          data[j+1] = bswap_16(data[j+1]);
	 }
       shortData[j] = (short)data[j];
       shortData[j+1] = (short)data[j+1];
      }
   }
  fout.write((char*)shortData,2*this->rangeSamples*sizeof(short));
  //fout.write((char*)shortData,2*18800*sizeof(short));

  /*delete [] data;
  delete [] shortData;*/
}
