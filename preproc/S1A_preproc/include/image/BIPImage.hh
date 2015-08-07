/*

Copyright (C) 2007, 2008, 2009, 2011 Walter M. Szeliga

This file is part of roipac2grdfile.

roipac2grdfile is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

roipac2grdfile is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with roipac2grdfile.  If not, see <http://www.gnu.org/licenses/>.

*/

/**
 * A class to allow access to Band interleaved by pixel images using a memory map
 */
template <typename T>
class BIPImage : public MultiBandImage
{
public:
  T getValue(int x, int y, int band)
  {
    this->testBand(band);
    this->testCoordinates(x,y);
    return this->image[y*this->width*this->bands + this->bands*x + x];
  }
 
  void setValue(int x, int y, int band, T val)
  {
    this->testBand(band);
    this->testCoordinates(x,y);
    this->image[y*this->width*this->bands + this->bands*x + x] = val;
  }
};
