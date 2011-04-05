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

#include <magnet/math/spline.hpp>
#include <magnet/clamp.hpp>
#include <magnet/exception.hpp>
#include <cmath>
#include <stdint.h>

namespace magnet {
  namespace color {
    class TransferFunction
    {
    public:
      TransferFunction():_valid(false) {}
      
      void addKnot(double x, double r, double g, double b, double a)
      { _knots.push_back(Knot(x,r,g,b,a)); _valid = false; }

      const std::vector<uint8_t>& getColorMap()
      {
	if (!_valid) generate();
	return _colorMap;
      }

    protected:
      void generate()
      {
	_colorMap.resize(4 * 256);

	math::Spline spline;
	spline.setLowBC(math::Spline::FIXED_1ST_DERIV_BC, 0);
	spline.setHighBC(math::Spline::FIXED_1ST_DERIV_BC, 0);

	for (size_t channel(0); channel < 4; ++channel)
	  {
	    spline.clear();

	    for (std::vector<Knot>::const_iterator iPtr = _knots.begin();
		 iPtr != _knots.end(); ++iPtr)
	      spline.addPoint(iPtr->_x, (*iPtr)(channel));

	    for (size_t i(0); i < 256; ++i)
	      _colorMap[4 * i + channel] = clamp(255 * spline(double(i) / 255.0), 0.0, 255.0);
	  }
	
	_valid = true;
      }
      

      struct Knot { 
	Knot(double x, double r, double g, double b, double a):
	  _x(clamp(x, 0.0, 1.0)), 
	  _r(clamp(r, 0.0, 1.0)), 
	  _g(clamp(g, 0.0, 1.0)), 
	  _b(clamp(b, 0.0, 1.0)),
	  _a(clamp(a, 0.0, 1.0)) 
	{}
	
	double operator()(size_t i) const
	{ switch (i)
	    {
	    case 0: return _r;
	    case 1: return _g;
	    case 2: return _b;
	    case 3: return _a;
	    }
	  M_throw() << "Bad index";
	}

	double _x, _r, _g, _b, _a;
      };

      std::vector<Knot> _knots;
      std::vector<uint8_t> _colorMap;
      bool _valid;
    };
  }
}
