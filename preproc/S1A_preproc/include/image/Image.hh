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
#ifndef IMAGE_HH
#define IMAGE_HH 1

#include <sys/fcntl.h>
#include <sys/mman.h>
#include <string>

/**
 * A class to allow access to raster images using a memory map
 */
template <typename T>
class Image
{
protected:
  size_t size;
  int width;
  int height;
  const char *filename;
  int openFlags;
  int mapFlags;
  int fd;
  T *image;
  /**
   * Create the memory map of the underlying raster image
   */
  void createMap()
  {    
    // If we are creating this image for the first time, we need to "create" space
    // for it on the drive
    if ( this->openFlags == (O_RDWR | O_CREAT) )
      {
	this->fd = open(this->filename, this->openFlags, (mode_t)0600);
	int status = ftruncate(this->fd,this->size);
	if (status == -1) {throw "Unable to create file";}
      }
    else
      {
	this->fd = open(this->filename, this->openFlags);
      }
    this->image = (T *)mmap(0, this->size, this->mapFlags, MAP_SHARED, this->fd,0);
    if (this->image == MAP_FAILED)
      {
	throw "Memory mapping failed";
      }
  }
  /**
   * Test whether the requested coordinates are in bounds
   */
  void testCoordinates(int x, int y)
  {
    if (x >= this->width)
      {
	throw "X coordinate out of bounds";
      }
    if (y >= this->height)
      {
	throw "Y coordinate out of bounds";
      }
  }
public:
  /**
   * The constructor
   * @param filename the name of the file containing the image
   * @param mode the mode with which to open the file, either 'r' for read or 'w' for write
   * @param width the width of the image in pixels
   * @param height the height of the image in pixels
   */
  Image(const char *filename,const char *mode,int width,int height)
  {
    this->filename = filename;
    this->width = width;
    this->height = height;
    
    std::string read = "r";
    std::string write = "w";
    // Convert the mode to an oflag for open and a flag for mmap
    if (read.compare(mode) == 0)
      {
	this->openFlags = O_RDONLY;
	this->mapFlags = PROT_READ;
      }
    else if (write.compare(mode) == 0)
      {
	this->openFlags = (O_RDWR | O_CREAT);
	this->mapFlags = (PROT_READ | PROT_WRITE);
      }
  }
 
  /**
   * The destructor
   */ 
  ~Image()
  { 
    munmap(this->image,this->size);
    close(this->fd);
  }
 
  /**
   * @return the image width
   */ 
  int getWidth()
  {
    return this->width;
  }
 
  /**
   * @return the image height
   */ 
  int getHeight()
  {
    return this->height;
  }

  virtual T getValue(int x,int y) = 0;
  virtual void setValue(int x, int y, T val) = 0;
};

#endif
