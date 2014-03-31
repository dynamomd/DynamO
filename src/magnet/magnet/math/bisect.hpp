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

namespace magnet {
  namespace math {
#ifdef MAGNET_DEBUG
    namespace {
      inline bool comparesign(double a, double b)
      { return (a < 0) == (b < 0); }
    }
#endif

    template<class Functor>
    inline double bisect(const Functor& func, double t1, double t2, double rootthreshold, const size_t nIt = 5000)
    {
#ifdef MAGNET_DEBUG
      if (comparesign(func(t1), func(t2)))
	M_throw() << "No sign change in the interval!";
#endif
      bool negative_min = func(t1) < 0;

      for (size_t i(0); i < nIt; ++i)
	{
	  double tm = 0.5 * (t1 + t2);
	  double f = func(tm);
	    
	  if (std::abs(f) < rootthreshold) return tm;

	  if ((f < 0.0) == negative_min)
	    t1 = tm;
	  else 
	    t2 = tm;
	}

      return t1;
    }
  }
}
