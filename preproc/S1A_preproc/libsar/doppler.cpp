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
#include "Doppler.hh"
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include "ezETAProgressBar.hpp"

Doppler::Doppler() {
    constant = 0.0;
    linear = 0.0;
    quadratic = 0.0;
}
Doppler::~Doppler() {}

void
Doppler::calculateDoppler(SingleBandImage<std::complex<unsigned char> > image, double iBias, double qBias)
{
	std::complex<double> zero(0.0,0.0);
	std::vector<std::complex<double> > proj(image.getWidth(),zero);
	ez::ezETAProgressBar pg(image.getHeight());
	pg.start();
	// Iterate over the image grabbing two range lines
	// at a time.  Calculate the product of the conjugate 
	// of line one with the value of line two.  At that value
	// to a running sum for each range bin.
	for (int i=0;i<image.getHeight()-1;i++)
	{
		for (int j=0;j<image.getWidth();j++)
		{
			std::complex<unsigned char> a = image.getValue(j,i);
			std::complex<unsigned char> b = image.getValue(j,i+1);
			std::complex<double> ra((double)real(a)-iBias,(double)imag(a)-qBias);
			std::complex<double> rb((double)real(b)-iBias,(double)imag(b)-qBias);
			proj[j] = proj[j] + conj(ra)*rb;
		}
		++pg;
	}

	std::cout << std::endl;	
	// Then, iterate over this sum by range bin, and calculate
	// The complex phase angle, atan2(Im,Re)
	// and divide by (2*Pi) This is the doppler estimate in the
	// range direction
	this->frac.clear();
	this->bin.clear();
	this->frac.resize(image.getWidth());
	this->bin.resize(image.getWidth());
	for (int i=0;i<image.getWidth();i++)
	{
		this->bin[i] = (double)i;
		this->frac[i] = atan2(imag(proj[i]),real(proj[i]))/(2.0*M_PI);
	}
}

void
Doppler::calculateDoppler(SingleBandImage<std::complex<float> > image)
{
	std::complex<double> zero(0.0,0.0);
	std::vector<std::complex<double> > proj(image.getWidth(),zero);
	ez::ezETAProgressBar pg(image.getHeight());
	pg.start();
	// Iterate over the image grabbing two range lines
	// at a time.  Calculate the product of the conjugate 
	// of line one with the value of line two.  At that value
	// to a running sum for each range bin.
	for (int i=0;i<image.getHeight()-1;i++)
	{
		for (int j=0;j<image.getWidth();j++)
		{
			std::complex<float> a = image.getValue(j,i);
			std::complex<float> b = image.getValue(j,i+1);
			std::complex<double> ra((double)real(a),(double)imag(a));
			std::complex<double> rb((double)real(b),(double)imag(b));
			proj[j] = proj[j] + conj(ra)*rb;
		}
		++pg;
	}

	std::cout << std::endl;	
	// Then, iterate over this sum by range bin, and calculate
	// The complex phase angle, atan2(Im,Re)
	// and divide by (2*Pi) This is the doppler estimate in the
	// range direction
	this->frac.clear();
	this->bin.clear();
	this->frac.resize(image.getWidth());
	this->bin.resize(image.getWidth());
	for (int i=0;i<image.getWidth();i++)
	{
		this->bin[i] = (double)i;
		this->frac[i] = atan2(imag(proj[i]),real(proj[i]))/(2.0*M_PI);
	}
}

void
Doppler::print(const std::string output)
{
  std::ofstream out;
  out.open(output.c_str());
  for (int i=0;i<((signed)this->frac.size());i++)
  {
	  out << (signed)this->bin[i] << " " << this->frac[i] << std::endl;
  }
  out.close();
}

void
Doppler::fit()
{
	double *ans = quadraticFit(this->bin,this->frac);
        this->constant  = ans[0];
	this->linear    = ans[1];
        this->quadratic = ans[2];	
}

double *
Doppler::quadraticFit(std::vector<double> x,std::vector<double> y)
{
 std::vector<double> sumX(5,0.0);
 std::vector<double> sumYX(3,0.0);

 for (int i=0;i<(signed)x.size();i++)
 {
	 sumX[0] += 1.0;
	 sumX[1] += x[i];
	 sumX[2] += x[i]*x[i];
	 sumX[3] += x[i]*x[i]*x[i];
	 sumX[4] += x[i]*x[i]*x[i]*x[i];
	 sumYX[0] += y[i];
	 sumYX[1] += y[i]*x[i];
	 sumYX[2] += y[i]*x[i]*x[i];
 }
 std::vector<double> A(9,0.0);
 A[0] = sumX[0], A[1] = sumX[1], A[2] = sumX[2];
 A[3] = sumX[1], A[4] = sumX[2], A[5] = sumX[3];
 A[6] = sumX[2], A[7] = sumX[3], A[8] = sumX[4];

 std::vector<double> invA = matrixInverse(A);

 double *ans = new double[3];

 ans[0] = invA[0 + 0 * 3]*sumYX[0] + invA[1 + 0 * 3]*sumYX[1] + invA[2 + 0 * 3]*sumYX[2];
 ans[1] = invA[0 + 1 * 3]*sumYX[0] + invA[1 + 1 * 3]*sumYX[1] + invA[2 + 1 * 3]*sumYX[2];
 ans[2] = invA[0 + 2 * 3]*sumYX[0] + invA[1 + 2 * 3]*sumYX[1] + invA[2 + 2 * 3]*sumYX[2];


 return ans;
}

std::vector<double>
Doppler::testInverse(std::vector<double> A)
{
    return matrixInverse(A);
}

std::vector<double>
Doppler::matrixInverse(std::vector<double> A)
{
    std::vector<double> x;
    x = A; // Make of copy of A

    int N = (int)sqrt(x.size());

    for (int i=0;i<N;i++)
    {
	    x[i + i*N] = 1.0/x[i + i*N];
	    for (int j=0;j<N;j++)
	    {
		    if (j==i) {continue;}
		    x[i + j*N] = x[i + j*N]*x[i + i*N];
	    }

	    for (int j=0;j<N;j++)
	    {
	      for (int k=0;k<N;k++)
	      {
		      if (j==i) {continue;}
		      if (k==i) {continue;}
		      x[k + j*N] = x[k +j*N] - x[i + j*N]*x[k + i*N];
	      }
	    }

	    for (int k=0;k<N;k++)
	    {
		    if (k==i) {continue;}
		    x[k + i*N] = -1.0*x[i + i*N]*x[k + i*N];
	    }
    }

    return x;
}

double
Doppler::getConstant()
{
	return this->constant;
}

double
Doppler::getLinear()
{
	return this->linear;
}

double
Doppler::getQuadratic()
{
	return this->quadratic;
}
