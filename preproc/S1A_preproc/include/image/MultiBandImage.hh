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

#include "Image.hh"

/**
 * A class to allow access to multi-band images using a memory map
 */
template <typename T>
class MultiBandImage : public Image<T>
{
protected:
  int bands;
  /**
   * Test whether the request band is valid
   */
  void testBand(int band)
  {
    if (band >= this->bands)
      {
	throw "band number out of bounds";
      }
  }
public:
  MultiBandImage(char *filename, const char *mode, int width, int height, int bands) : Image<T>(filename,mode,width,height)
  {
    this->bands = bands;
    this->size = (size_t)(this->width*this->height*this->bands*sizeof(T));
    this->createMap();
  }

  /**
   * @return the number of bands
   */
  int getBands()
  {
    return this->bands;
  }

  /**
   * Retrieve a value from the image.
   * @param x the pixel coordinate in the width direction
   * @param y the pixel coordinate in the height direction
   * @param band the band number
   * @return the value of the pixel at location (x,y)
   */ 
  virtual T getValue(int x, int y , int bands) =0;

  /**
   * Set the value of a pixel in an image
   * @param x the pixel coordinate in the width direction
   * @param y the pixel coordinate in the height direction 
   * @param band the band number
   * @param val the value to set the pixel at (x,y) to
   */
  virtual void setValue(int x, int y, int bands, T val) =0;
};
