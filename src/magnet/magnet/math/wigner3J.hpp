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

#include <boost/math/special_functions.hpp>

namespace magnet {
  namespace math {
    double wignerThreej(const int & la, const int & lb, 
			const int & lc, const int & ma, 
			const int & mb, const int & mc)
    {
      //Check the sum of the m's is zero
      if (ma + mb + mc != 0) return 0.0;
      
      int numin = 0, numax = la-ma;
      
      {
	int tmp = -ma+lb-lc;
	if (tmp > numin) numin = tmp;
	
	tmp = la+mb-lc;
	if (tmp > numin) numin = tmp;
	
	tmp = lb+mb;
	if (tmp < numax) numax = tmp;
	
	tmp = la+lb-lc;
	if (tmp < numax) numax = tmp;
      }
      
      int sign(1);
      if ((numin % 2)) sign = -1;
      
      double retval = 0;
      
      using namespace boost::math;
      
      for (int nu(numin); nu <= numax; nu++)
	{
	  retval += sign 
	    / (factorial<double>(la-ma-nu)*factorial<double>(lc-lb+ma+nu)
	       *factorial<double>(lb+mb-nu)*factorial<double>(lc-la-mb+nu)
	       *factorial<double>(nu)*factorial<double>(la+lb-lc-nu));
	  sign = -sign;
	}
      
      retval *= std::sqrt(factorial<double>(la+lb-lc)
			  *factorial<double>(la+lc-lb)
			  *factorial<double>(lb+lc-la)
			  / factorial<double>(la+lb+lc+1));
      
      retval *= std::sqrt(factorial<double>(la+ma)*factorial<double>(lb+mb)
			  *factorial<double>(lc+mc)*factorial<double>(la-ma)
			  *factorial<double>(lb-mb)*factorial<double>(lc-mc));
      
      if ((abs(la-lb-mc) % 2) != 0) retval = -retval;
      
      return retval;
    }
    
  }
}
