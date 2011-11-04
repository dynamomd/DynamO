/*  dynamo:- Event driven molecular dynamics simulator 
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
    template<class T>
    inline void RGBtoHSV(T color[4], T r, T g = 1, T b = 1, T alpha = 1) 
    {
      T min = std::min(r, std::min(g, b));
      T max = std::max(r, std::max(g, b));

      color[0] = 0; //h
      color[1] = 0; //s
      color[2] = max; //v
      color[3] = alpha; //a

      if (max == 0) return; //Early exit

      r /= max;
      g /= max;
      b /= max;
      min /= max;
      max = 1;

      T delta = max - min;
      color[1] = delta; //s

      if (delta == 0)
	return;

      r = (r - min) / delta;
      g = (g - min) / delta;
      b = (b - min) / delta;
      min = 0; max = 1;

      if ((r > g) && (r > b))
	color[0] = (g - b) / 6.0;
      else if ((g > r) && (g > b))
	color[0] = (1.0/3.0) + (b - r) / 6.0;
      else
	color[0] = (2.0/3.0) + (r - g) / 6.0;

      color[0] += (color[0] < 0);
    }


    template<class T>
    inline void HSVtoRGB(T color[4], T h, T s = 1, T v = 1, T alpha = 1) 
    {
      T temp;
      h = std::modf(clamp(h,0.0f, 1.0f), &temp) * 6;
      s = clamp(s, 0.0f, 1.0f);
      v = clamp(v, 0.0f, 1.0f);
      
      unsigned int i = h;
      T f = h - i;
      T p = v * (1 - s);
      T q = v * (1 - s * f);
      T t = v * (1 - s * (1 - f));
      
      T& r = color[0];
      T& g = color[1];
      T& b = color[2];
      r = g = b = 0;

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
      color[3] = alpha;
    }

    inline void HSVtoRGB(cl_uchar4& color, float h, float s = 1, float v = 1)
    {
      GLfloat floatcolor[4];
      HSVtoRGB(floatcolor, h, s, v);

      for (size_t i(0); i < 4; ++i)
	color.s[i] = 255 * floatcolor[i];
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
