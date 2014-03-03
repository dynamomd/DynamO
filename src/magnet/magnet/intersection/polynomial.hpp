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
#include <algorithm>

namespace magnet {
  namespace intersection {
    namespace detail {
      inline double firstOrder(const double f0, const double f1) {
	if (f1 >= 0) return HUGE_VAL;
	return std::max(0.0, - f0 / f1);
      }
      
      inline double secondOrderPositive(const double f0, const double f1, const double f2)
      {
	double arg = f1 * f1 - 2 * f2 * f0;
	if ((f1 >= 0) || arg <= 0) return HUGE_VAL;
	return std::max(0.0, 2 * f0 / (-f1 + std::sqrt(arg)));
      }

      inline double secondOrderNegative(const double f0, const double f1, const double f2)
      {
	double arg = f1 * f1 - 2 * f2 * f0;
	if (arg <= 0) return std::max(0.0, -f1 / f2);
	if (f1 > 0)
	  return std::max(0.0, (-f1 - std::sqrt(arg)) / f2);
	else
	  return std::max(0.0, 2 * f0 / (-f1 + std::sqrt(arg)));
      }

      inline double secondOrder(const double f0, const double f1, const double f2)
      {
	if (f2 == 0) return firstOrder(f0, f1);
	else if (f2 < 0) return secondOrderNegative(f0,f1,f2);
	else /*(f2 > 0)*/ return secondOrderPositive(f0,f1,f2);
      }

      /* class F { template<size_t deriv> double eval(double deltat); template<size_t deriv, bool decreasing> nextRoot(double deltat); } */

      //template<class F>
      //inline double stableAlgorithm(const F& f)
      //{
      //	if (f.eval<0>(0) > 0) return f.template nextRoot<0, true>(0);
      //  if (f.eval<1>(0) < 0) return 0;
      //	const double tderiv = f.template nextRoot<1, true>(0);
      //	if (tderiv == HUGE_VAL) return HUGE_VAL;
      //	if (f.eval<0>(tderiv) > 0) return f.template nextRoot<0, true>(tderiv);
      //}
    }
  }
}
