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
#include <magnet/math/cubic.hpp>

namespace magnet {
  namespace intersection {
    namespace detail {
      inline double firstOrder(const double f0, const double f1) {
	if (f1 >= 0) return HUGE_VAL;
	return std::max(0.0, - f0 / f1);
      }
      
      inline double secondOrder(const double f0, const double f1, const double f2)
      {
	if (f2 == 0) return firstOrder(f0, f1);
	const double arg = f1 * f1 - 2 * f2 * f0;
	if (f2 < 0) {
	  if (arg <= 0) return std::max(0.0, -f1 / f2);
	  if (f1 > 0)
	    return std::max(0.0, (-f1 - std::sqrt(arg)) / f2);
	  else
	    return std::max(0.0, 2 * f0 / (-f1 + std::sqrt(arg)));
	} else /*(f2 > 0)*/ {
	  const double arg = f1 * f1 - 2 * f2 * f0;
	  if ((f1 >= 0) || arg <= 0) return HUGE_VAL;
	  return std::max(0.0, 2 * f0 / (-f1 + std::sqrt(arg)));
	}
      }


      inline double thirdOrder(const double f0, const double f1, const double f2, const double f3)
      {
	if (f3 == 0) return secondOrder(f0, f1, f2);

	//Calculate and sort the roots of the overlap function
	std::array<double, 3> roots;
	const size_t nroots = magnet::math::cubicSolve(6 * f2 / (2 * f3), 6 * f1 / f3, 6 * f0 / f3, roots[0], roots[1], roots[2]);
	std::sort(roots.begin(), roots.begin() + nroots);

	//Calculate and sort the roots of the overlap function's
	//derivative
	std::array<double, 2> derivroots;
	const size_t nderivroots
	  = magnet::math::quadSolve(f3, f2, f1, derivroots[0], derivroots[1]);
	if (derivroots[1] < derivroots[0])
	  std::swap(derivroots[0], derivroots[1]);
	
	if (f3 > 0) {
	  if (nderivroots == 0) return HUGE_VAL;
	  if (derivroots[1] < 0) return HUGE_VAL;
	  if ((nroots==1) && (roots[0] < derivroots[0])) return HUGE_VAL;
	  return std::max(0.0, (nroots==3) ? roots[1]: derivroots[0]);
	} else /*(f3 < 0)*/ {
	  if (nderivroots == 0) return std::max(0.0, roots[0]);
	  if ((derivroots[0] > 0) && (roots[0] < derivroots[0]))
	    return std::max(0.0, roots[0]);
	  return std::max(0.0, std::max(derivroots[1], roots[nroots-1]));
	}
      }

      inline double fourthOrder(const double f0, const double f1, const double f2, const double f3, const double f4, double f0char, double precision=1e-16)
      {
	if (f4 == 0) return thirdOrder(f0, f1, f2, f3);
	
	//Determine and sort the roots of the derivative
	std::array<double, 3> roots;
	const size_t rootCount = magnet::math::cubicSolve(3 * f3 / f4, 6 * f2 / f4, 6 * f1 / f4, roots[0], roots[1], roots[2]);
	std::sort(roots.begin(), roots.begin() + rootCount);

	//The overlap function
	auto f = [&] (double t) { return (((f4 * t / 4 + f3) * t / 3 + f2) * t /2 + f1) * t + f0; };
	
	const double rootthreshold = f0char * precision;
	
	if (f4 > 0) {
	  if ((roots[0] > 0) && (f(roots[0]) < 0))
	    {
	      if (f(0) <= 0)
		return 0;
	      else
		return magnet::math::bisect(f, 0, roots[0], rootthreshold);
	    }
	  
	  if ((rootCount == 3) && (roots[2] > 0) && (f(roots[2]) < 0))
	    {
	      double tmin = std::max(0.0, roots[1]);
	      if (f(tmin) <= 0)
		return tmin;
	      else
		return magnet::math::bisect(f, tmin, roots[2], rootthreshold);
	    }
	  return HUGE_VAL;
	} else /*(f4 < 0)*/ {
	  if ((rootCount == 3) && (roots[1] > 0) && (f(roots[1]) < 0))
	    {
	      double tmin = std::max(0.0, roots[0]);
	      if (f(tmin) <= 0)
		return tmin;
	      else
		return magnet::math::bisect(f, tmin, roots[1], rootthreshold);
	    }
      
	  double tlast = roots[rootCount - 1];
      
	  if (f(tlast) <= 0) return std::max(0.0, tlast);

	  if ((tlast < 0) && (f(0) <= 0)) return 0;
      
	  double t0 = std::max(0.0, tlast);
      
	  double deltate = std::pow(- 24 * f0char / f4, 0.25);
      
	  while (f(t0 + deltate) > 0)
	    {
	      t0 += deltate;
	      deltate *= 2;
	    }
	  return magnet::math::bisect(f, t0, t0+deltate, rootthreshold);
	}
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
