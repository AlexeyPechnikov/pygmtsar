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
#ifndef DOPPLER_HH
#define DOPPLER_HH 1

#include <complex>
#include <vector>
#include "image/SingleBandImage.hh"

class Doppler
{
  private:
	  double constant;
	  double linear;
	  double quadratic;
	  std::vector<double> bin;
	  std::vector<double> frac;
          double *quadraticFit(std::vector<double> x, std::vector<double> y);
	  std::vector<double> matrixInverse(std::vector<double> A);
  public:
	  Doppler();
	  ~Doppler();
	  void calculateDoppler(SingleBandImage<std::complex<float> > image);
	  void calculateDoppler(SingleBandImage<std::complex<unsigned char> > image, double iBias, double qBias);
	  void print(const std::string output);
	  void fit();
	  double getConstant();
	  double getLinear();
	  double getQuadratic();
	  std::vector<double> testInverse(std::vector<double> A);
};

#endif
