/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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

#include <cmath>
#include <magnet/clamp.hpp>
#include <magnet/exception.hpp>

namespace magnet {
  namespace color {
    inline void HSVtoRGB(cl_uchar4& color, float h, float s = 1, float v = 1) 
    {
      float temp;
      h = std::modf(h, &temp);
      
      s = clamp(s, 0.0f, 1.0f);
      v = clamp(v, 0.0f, 1.0f);
      
      h = h * 6;
      
      unsigned int i = h;
      float f = h - i;
      float p = v * (1 - s);
      float q = v * (1 - s * f);
      float t = v * (1 - s * (1 - f));
      
      float r;
      float g;
      float b;
      
      switch(i) {
      case 0:	
	r = v;
	g = t;
	b = p;
	break;
	
      case 1:	
	r = q;
	g = v;
	b = p;
	break;
	
      case 2:
	r = p;
	g = v;
	b = t;
	break;
	
      case 3:
	r = p;
	g = q;
	b = v;
	break;
	
      case 4:
	r = t;
	g = p;
	b = v;
	break;
	
      case 5:
	r = v;
	g = p;
	b = q;
	break;
      }
      color.s[0] = r * 255;
      color.s[1] = g * 255;
      color.s[2] = b * 255;
      color.s[3] =     255;
    }

    inline std::string getOpenCLHSV()
    {
#define STRINGIFY(A) #A

      return STRINGIFY(
    void HSVtoRGB(__local uchar4* color, float h, float s, float v) 
{
  float temp;
  h = fract(h, &temp);

  s = clamp(s, 0.0, 1.0);
  v = clamp(v, 0.0, 1.0);

  h = h * 6;

  unsigned int i = h;
  float f = h - i;
  float p = v * (1 - s);
  float q = v * (1 - s * f);
  float t = v * (1 - s * (1 - f));

  float r;
  float g;
  float b;

  switch(i) {
  case 0:	
    r = v;
    g = t;
    b = p;
    break;

  case 1:	
    r = q;
    g = v;
    b = p;
    break;

  case 2:
    r = p;
    g = v;
    b = t;
    break;

  case 3:
    r = p;
    g = q;
    b = v;
    break;

  case 4:
    r = t;
    g = p;
    b = v;
    break;

  case 5:
    r = v;
    g = p;
    b = q;
    break;
  }
  *color = (uchar4)(r*255,g*255,b*255,255);
});
#undef STRINGIFY
    }
  }
}
