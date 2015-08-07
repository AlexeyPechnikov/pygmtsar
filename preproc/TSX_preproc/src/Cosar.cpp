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
#include "Cosar.hh"

// Cosar files are big-endian
// thus, we need to determine the endianness of the machine we are on
// and decide whether we need to swap bytes
Cosar::Cosar(std::string input, std::string output)
{
 // Check the endianness
 if (is_big_endian() == 1) 
   {
     std::cout << "Machine is Big Endian" << std::endl;
     this->isBigEndian = true;
   }
 else
   {
     std::cout << "Machine is Little Endian" << std::endl;
     this->isBigEndian = false;
   }
 this->fin.open(input.c_str(), std::ios::binary | std::ios::in);
 if (fin.fail()) { throw "Unable to open file " + input;}
 this->fout.open(output.c_str(), std::ios::binary | std::ios::out);
 if (fout.fail()) { throw "Unable to open file " + output;}
 try {
   this->header = new Header(this->isBigEndian);
 } catch(char *ex) {
     throw;
 }
}

Cosar::~Cosar()
{
  this->fin.close();
  this->fout.close();
}

void
Cosar::parse()
{
 this->header->parse(this->fin);
 int byteTotal = this->header->getRangelineTotalNumberOfBytes();
 int numLines = this->header->getTotalNumberOfLines();
 int burstSize = this->header->getBytesInBurst();
 int rangeSamples = this->header->getRangeSamples();
 int azimuthSamples = this->header->getAzimuthSamples();

 this->numberOfBursts = (int)(byteTotal*numLines)/burstSize;
 this->bursts = new Burst*[this->numberOfBursts];
 for(int i=0;i<this->numberOfBursts;i++)
   {
     this->bursts[i] = new Burst(rangeSamples,azimuthSamples,this->isBigEndian);
     this->bursts[i]->parse(this->fin,this->fout);
   }
}
