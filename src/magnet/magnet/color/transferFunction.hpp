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

#include <magnet/math/spline.hpp>
#include <magnet/clamp.hpp>
#include <magnet/exception.hpp>
#include <magnet/color/HSV.hpp>
#include <cmath>
#include <tr1/array>
#include <stdint.h>

namespace magnet {
  namespace color {
    namespace detail {
      struct Knot { 
	Knot(double x, double h, double s, double v, double a):
	  _x(clamp(x, 0.0, 1.0)), 
	  _h(clamp(h, 0.0, 1.0)), 
	  _s(clamp(s, 0.0, 1.0)), 
	  _v(clamp(v, 0.0, 1.0)),
	  _a(clamp(a, 0.0, 1.0)) 
	{}
	
	const double& operator()(size_t i) const
	{ switch (i)
	    {
	    case 0: return _h;
	    case 1: return _s;
	    case 2: return _v;
	    case 3: return _a;
	    }
	  M_throw() << "Bad index";
	}

	double& operator()(size_t i)
	{ switch (i)
	    {
	    case 0: return _h;
	    case 1: return _s;
	    case 2: return _v;
	    case 3: return _a;
	    }
	  M_throw() << "Bad index";
	}
	
	bool operator<(const Knot& ok) const { return _x < ok._x; }
	
	double _x, _h, _s, _v, _a;
      };
    }

    class TransferFunction: protected std::vector<detail::Knot>
    {
    public:      
      typedef detail::Knot Knot;
      typedef std::vector<Knot> Base;
      typedef Base::const_iterator const_iterator;

      TransferFunction():_valid(false) 
      {
	for (size_t channel(0); channel < 4; ++channel)
	  {
	    spline[channel].setLowBC(math::Spline::FIXED_1ST_DERIV_BC, 0);
	    spline[channel].setHighBC(math::Spline::FIXED_1ST_DERIV_BC, 0);
	  }

	for (size_t channel(0); channel < 4; ++channel)
	  spline[channel].setType(math::Spline::LINEAR);

      }
      
      void addKnot(double x, double h, double s, double v, double a)
      { 
	push_back(Knot(x, h, s, v, a)); 
	_valid = false; 
      }

      const std::vector<float> getMap(size_t samples, float transmittanceScale)
      {
	if (!_valid) generate();

	std::vector<float> colorMap;
	colorMap.resize(4 * samples);

	for (size_t i(0); i < samples; ++i)
	  {
	    double transmittance = magnet::clamp(spline[3](double(i) / (samples - 1)), 0.0, 1.0);
	    transmittance *= transmittance;
	    transmittance *= transmittance;
	    transmittance *= transmittance;
	    transmittance *= transmittanceScale;

	    HSVtoRGB<float>(&(colorMap[4 * i]),
			    magnet::clamp(spline[0](double(i) / (samples - 1)), 0.0, 1.0),
			    magnet::clamp(spline[1](double(i) / (samples - 1)), 0.0, 1.0), 
			    magnet::clamp(spline[2](double(i) / (samples - 1)), 0.0, 1.0),
			    transmittance);
	  }
	
	return colorMap;
      }
      
      /*! \brief Calculate the pre-integrated color map.
       */
      std::vector<float> getPreIntegratedMap(size_t samples, float transmittanceScale)
      {
	if (!_valid) generate();

	std::vector<float> colorMap = getMap(samples, transmittanceScale);
	std::vector<float> integral;
	integral.resize(4 * samples);

	//Initial value
	for (size_t c(0); c < 4; ++c)
	  integral[4 * 0 + c] = colorMap[4 * 0 + c];

	for (size_t i(1); i < samples; ++i)
	  {
	    for (size_t c(0); c < 3; ++c)
	      integral[4 * i + c] = integral[4 * (i - 1) + c] 
		+ colorMap[4 * i + c] * colorMap[4 * i + 3];

	    integral[4 * i + 3] = integral[4 * (i - 1) + 3] + colorMap[4 * i + 3];
	  }

	//Normalization
	for (size_t i(0); i < samples; ++i)
	  for (size_t c(0); c < 4; ++c)
	    integral[4 * i + c] /= (samples-1);

	return integral;
      }

      void addInterpolatedKnot(float x)
      {
	addKnot(x, spline[0](x), spline[1](x), spline[2](x), spline[3](x));
      }

      std::tr1::array<float,4> getValue(float x)
      {
	std::tr1::array<float,4> retval = {{spline[0](x), spline[1](x), spline[2](x), spline[3](x)}};
	return retval;
      }

      void setKnot(const const_iterator& it, const Knot& val)
      {
	*(Base::begin() + (it - begin())) = val;
	_valid = false; 
      }

      void eraseKnot(const const_iterator& it)
      {
	erase(Base::begin() + (it - begin()));
	_valid = false; 
      }

      void setKnot(size_t index, const Knot& val)
      {
	*(Base::begin() + index) = val;
      }

      const_iterator begin() const { return Base::begin(); }
      const_iterator end() const { return Base::end(); }
      void clear() { _valid = false; Base::clear(); }
      size_t size() const { return Base::size(); }
      size_t max_size() const { return Base::max_size(); }
      size_t capacity() const { return Base::capacity(); }
      bool empty() const { return Base::empty(); }


    protected:
      void generate()
      {
	for (size_t channel(0); channel < 4; ++channel)
	  {
	    spline[channel].clear();
	    
	    for (const_iterator iPtr = begin(); iPtr != end(); ++iPtr)
	      spline[channel].addPoint(iPtr->_x, (*iPtr)(channel));
	  }
	
	_valid = true;
      }

      bool _valid;
      math::Spline spline[4];
    };
  }
}
