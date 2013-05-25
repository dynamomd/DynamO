/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

    This program is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    version 3 as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#pragma once

#include <vector>
#include <cmath>
#include <limits>

namespace magnet {
  namespace image {

    /*! \brief Class which calculates the signed distances between
     * pixels inside a bitmap.
     *
     * This class uses the 3x3 dead reckoning algorithm, as presented
     * in "The ‘‘dead reckoning’’ signed distance transform," by
     * G. J. Grevera (doi:10.1016/j.cviu.2004.05.002).
     *
     * This function is a handy image processing algorithm to convert
     * bitmapped images into scalable forms. This is particularly
     * effective when trying to render vector graphics in an OpenGL
     * scene (e.g., see \ref CairoSurface).
     */
    class SignedDistanceTransform
    {
    public:
      SignedDistanceTransform(unsigned char* I, size_t width, size_t height):
	_width(width), _height(height)
      {
	//The border of the image must be 0
	for (size_t x(0); x < _width; ++x)
	  { I[iPos(x, 0)] = 0; I[iPos(x, _height - 1)] = 0; }
	for (size_t y(0); y < _height; ++y)
	  { I[iPos(0, y)] = 0; I[iPos(_width - 1, y)] = 0; }
	  
	double floatInf = std::numeric_limits<float>::max() / 2;

	p.resize(_width * _height);
	d.resize(_width * _height, floatInf);
	  
	//Iterate over the image's non-border pixels, finding color
	//edge pixels. 
	for (size_t y(1); y < _height - 1; ++y)
	  for (size_t x(1); x < _width - 1; ++x)
	    {
	      interpolateDistance(x, y, x - 1, y, I);
	      interpolateDistance(x, y, x + 1, y, I);
	      interpolateDistance(x, y, x, y - 1, I);
	      interpolateDistance(x, y, x, y + 1, I);
	    }

	//First "forward" pass
	for (int y(1); y < int(_height) - 1; ++y)
	  for (int x(1); x < int(_width) - 1; ++x)
	    {
	      int i2 = iPos(x,y);
	      check(x-1,y-1,std::sqrt(2), i2);
	      check(x,y-1,1, i2);
	      check(x+1,y-1,std::sqrt(2), i2);
	      check(x-1,y,1, i2);
	    }

	//Second "reverse" pass
	for (int y(_height - 2); y > 0; --y)
	  for (int x(_width - 2); x > 0; --x)
	    {
	      size_t i2 = iPos(x,y);
	      check(x+1,y,1, i2);
	      check(x-1,y+1,std::sqrt(2), i2);
	      check(x,y+1,1, i2);
	      check(x+1,y+1,std::sqrt(2), i2);
	    }

	//Transform the output to a scaled range
	for (size_t x(0); x < _width; ++x)
	  for (size_t y(0); y < _height; ++y)
	    {
	      size_t i = iPos(x,y);
	      double locd = (I[i] > 127) ?  d[i] : -d[i];
	      I[i] = std::min(255.0, std::max(0.0, 128 + locd));
	    }
      }
      
    protected:
      struct P: public std::array<int,2> 
      { P() { (*this)[0] = (*this)[1] = -1; } };

      size_t _width, _height;
      std::vector<P> p;
      std::vector<double> d;

      inline size_t iPos(size_t x, size_t y) { return y * _width + x; }
      
      inline void check(int x, int y, double delta, size_t i2)
      {
	int i1 = iPos(x,y);
	if (d[i1] + delta < d[i2]) 
	  {						      
	    p[i2] = p[i1];			              
	    const double deltaX = (p[i1][0] - x);
	    const double deltaY = (p[i1][1] - y);
	    d[i2] = std::sqrt(deltaX*deltaX + deltaY*deltaY);
	  }
      }
      
      void interpolateDistance(size_t x1, size_t y1, size_t x2, size_t y2, const unsigned char* I)
      {
	int i1 = iPos(x1, y1);
	if ((I[iPos(x2, y2)] > 127) != (I[i1] > 127))
	  {
	    //double interpDist = std::abs((128.0 - I[i1]) / (double(I[iPos(x2, y2)]) - I[i1]));
	    d[i1] = 1;//std::min(d[i1], interpDist);
	    p[i1][0] = x2;
	    p[i1][1] = y2;
	  }
      }
    };
  }
}
