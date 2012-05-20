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
    template<class Functor, class T = double>
    struct Bisect: public Functor {
      inline T bisectRoot(double t1, double t2, double rootthreshold, const size_t nIt = 500)
      {
#ifdef MAGNET_DEBUG
	if ((Functor::operator()(t1) < 0) == (Functor::operator()(t2) < 0))
	  M_throw() << "No sign change in the interval!";
	
	if ((Functor::operator()(t1) < 0))
	  M_throw() << "bisecting from negative to positive!";
#endif

	for(size_t i = 0; i < nIt; ++i)
	  {
	    T tm = 0.5 * (t1 + t2);
	    T f = Functor::operator()(tm);
	    if ((std::abs(f) < rootthreshold) && f > 0.0)
	      {
	        t1 = tm;
	        break;
	      }
	    
	    if (f < 0.0)
	      t2 = tm;
	    else 
	      t1 = tm;
	  }

	return t1;
      }
    };
  }
}
