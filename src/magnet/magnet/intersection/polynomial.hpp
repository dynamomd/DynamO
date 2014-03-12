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
      template<size_t Order>
      struct PolynomialFunction: public std::array<double, Order + 1>
      {
	template <typename... T> 
	PolynomialFunction(T... ts) : std::array<double, Order + 1>{{ts...}} {} 

	void timeShift(double dt) {
	  std::array<double, Order + 1> newCoeffs;
	  for (size_t i(0); i < Order + 1; ++i)
	    newCoeffs[i] = eval(dt, i);
	}

	double eval(double dt, size_t derivative = 0) const {
	  double accum(0);
	  for (size_t i(Order); i > derivative; --i) {
	    accum = (*this)[i] + (dt * accum) / (i + 1 - derivative);
	  }
	  return accum * dt + (*this)[derivative];
	}
	
	void flipSign() {
	  for (size_t i(0); i < Order + 1; ++i)
	    (*this)[i] = -(*this)[i];
	}
      };
      
      template<size_t Order>
      PolynomialFunction<Order-1> lowerOrder(const PolynomialFunction<Order>& f) {
	PolynomialFunction<Order-1> retval;
	for (size_t i(0); i < Order; ++i)
	  retval[i] = f[i];
	return retval;
      }

      inline double nextEvent(const PolynomialFunction<1>& f) {
	if (f[1] >= 0) return HUGE_VAL;
	return std::max(0.0, - f[0] / f[1]);
      }

      inline double nextEvent(const PolynomialFunction<2>& f) {
	if (f[2] == 0) return nextEvent(lowerOrder(f));
	const double arg = f[1] * f[1] - 2 * f[2] * f[0];
	if (f[2] < 0) {
	  if (arg <= 0) return std::max(0.0, -f[1] / f[2]);
	  if (f[1] > 0)
	    return std::max(0.0, (-f[1] - std::sqrt(arg)) / f[2]);
	  else
	    return std::max(0.0, 2 * f[0] / (-f[1] + std::sqrt(arg)));
	} else /*(f[2] > 0)*/ {
	  const double arg = f[1] * f[1] - 2 * f[2] * f[0];
	  if ((f[1] >= 0) || arg <= 0) return HUGE_VAL;
	  return std::max(0.0, 2 * f[0] / (-f[1] + std::sqrt(arg)));
	}
      }

      inline double nextEvent(const PolynomialFunction<3>& f)
      {
	if (f[3] == 0) return nextEvent(lowerOrder(f));

	//Calculate and sort the roots of the overlap function
	std::array<double, 3> roots;
	const size_t nroots = magnet::math::cubicSolve(6 * f[2] / (2 * f[3]), 6 * f[1] / f[3], 6 * f[0] / f[3], roots[0], roots[1], roots[2]);
	std::sort(roots.begin(), roots.begin() + nroots);

	//Calculate and sort the roots of the overlap function's
	//derivative
	std::array<double, 2> derivroots;
	const size_t nderivroots = magnet::math::quadSolve(f[3], f[2], f[1], derivroots[0], derivroots[1]);
	if (derivroots[1] < derivroots[0])
	  std::swap(derivroots[0], derivroots[1]);
	
	if (f[3] > 0) {
	  if (nderivroots == 0) return HUGE_VAL;
	  if (derivroots[1] < 0) return HUGE_VAL;
	  if ((nroots==1) && (roots[0] < derivroots[0])) return HUGE_VAL;
	  return std::max(0.0, (nroots==3) ? roots[1] : derivroots[0]);
	} else /*(f[3] < 0)*/ {
	  if (nderivroots == 0) return std::max(0.0, roots[0]);
	  if ((derivroots[0] > 0) && (roots[0] < derivroots[0]))
	    return std::max(0.0, roots[0]);
	  return std::max(0.0, std::max(derivroots[1], roots[nroots-1]));
	}
      }

      inline double nextEvent(const PolynomialFunction<4>& f, double f0char, double precision=1e-16)
      {
	if (f[4] == 0) return nextEvent(lowerOrder(f));
	
	//Determine and sort the roots of the derivative
	std::array<double, 3> roots;
	const size_t rootCount = magnet::math::cubicSolve(3 * f[3] / f[4], 6 * f[2] / f[4], 6 * f[1] / f[4], roots[0], roots[1], roots[2]);
	std::sort(roots.begin(), roots.begin() + rootCount);

	const double rootthreshold = f0char * precision;
	
	auto bisectFunc = [&] (double t) { return f.eval(t); };

	if (f[4] > 0) {
	  if ((roots[0] > 0) && (f.eval(roots[0]) < 0))
	    {
	      if (f.eval(0) <= 0)
		return 0;
	      else
		return magnet::math::bisect(bisectFunc, 0, roots[0], rootthreshold);
	    }
	  
	  if ((rootCount == 3) && (roots[2] > 0) && (f.eval(roots[2]) < 0))
	    {
	      double tmin = std::max(0.0, roots[1]);
	      if (f.eval(tmin) <= 0)
		return tmin;
	      else
		return magnet::math::bisect(bisectFunc, tmin, roots[2], rootthreshold);
	    }
	  return HUGE_VAL;
	} else /*(f4 < 0)*/ {
	  if ((rootCount == 3) && (roots[1] > 0) && (f.eval(roots[1]) < 0))
	    {
	      double tmin = std::max(0.0, roots[0]);
	      if (f.eval(tmin) <= 0)
		return tmin;
	      else
		return magnet::math::bisect(bisectFunc, tmin, roots[1], rootthreshold);
	    }
      
	  double tlast = roots[rootCount - 1];
      
	  if (f.eval(tlast) <= 0) return std::max(0.0, tlast);

	  if ((tlast < 0) && (f.eval(0) <= 0)) return 0;
      
	  double t0 = std::max(0.0, tlast);
      
	  double deltate = std::pow(- 24 * f0char / f[4], 0.25);
      
	  while (f.eval(t0 + deltate) >= 0)
	    {
	      t0 += deltate;
	      deltate *= 2;
	    }
	  return magnet::math::bisect(bisectFunc, t0, t0+deltate, rootthreshold);
	}
      }
    }
  }
}
