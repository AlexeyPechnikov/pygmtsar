/*

Copyright (C) 2007, 2008, 2009, 2011, 2012 Walter M. Szeliga

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
#ifndef SINGLEBANDIMAGE_HH
#define SINGLEBANDIMAGE_HH 1

#include "Image.hh"

/**
 * A class to allow access to single banded images using a memory map
 */
template <typename T>
class SingleBandImage : public Image<T>
{
public:
  /**
   * The constructor
   * @param filename the name of the file containing the image
   * @param mode the mode with which to open the file, either 'r' for read or 'w' for write
   * @param width the width of the image in pixels
   * @param height the height of the image in pixels
   */
  SingleBandImage(const char *filename,const char *mode,int width,int height) : Image<T>(filename, mode, width, height)
  {
    this->size = (size_t)(this->width*this->height*sizeof(T));
    this->createMap();
  }
 
  /**
   * The destructor
   */ 
  ~SingleBandImage() {}

  /**
   * Retrieve a value from the image.
   * @param x the pixel coordinate in the width direction
   * @param y the pixel coordinate in the height direction
   * @return the value of the pixel at location (x,y)
   */ 
  T getValue(int x, int y)
  {
    this->testCoordinates(x,y);
    return this->image[y*this->width + x];
  }
 
  /**
   * Set the value of a pixel in an image
   * @param x the pixel coordinate in the width direction
   * @param y the pixel coordinate in the height direction 
   * @param val the value to set the pixel at (x,y) to
   */
  void setValue(int x, int y, T val)
  {
    this->testCoordinates(x,y);
    this->image[y*this->width + x] = val;
  }

  /**
   * Set the values for an entire image row
   * @param y the pixel row
   * @param row a pointer to width values
   */
 void setRow(int y, T *row)
 {
	 size_t n = (size_t)(this->width*sizeof(T));
	 memcpy(&(this->image[y*this->width]),row,n);
 }

};

#endif
