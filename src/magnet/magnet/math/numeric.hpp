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
#include <magnet/math/operators.hpp>
#include <magnet/math/precision.hpp>
#include <magnet/exception.hpp>
#include <boost/math/tools/roots.hpp>

namespace magnet {
  namespace math {
    namespace detail {
      /*! \brief Update and checking of safeguards, called after an
	iterative step is taken towards a root.

	Returns 0 on error, 1 on successful, and 2 on converged.
      */
      template<class F, size_t Derivatives, class Real>
      inline int process_iterative_step(const F& f, std::array<Real, Derivatives>& curr_state,
					Real& x, Real new_x, Real& low_bound, Real& high_bound, Real x_precision) {

	//std::cout << "new_x = " << new_x << std::endl;
	//std::cout << "[" << low_bound << "," << high_bound << "]" << std::endl;

	if ((new_x >= high_bound) || (new_x <= low_bound)) {
	  //std::cout << "Out of bounds" << std::endl;
	  //Out of bounds step, failed!
	  return 0;
	}
	
	//Re evalulate the derivatives
	auto new_state = f(new_x);

	//Check for convergence
	Real delta = new_x - x;

	if ((std::abs(delta) < std::abs(x_precision * new_x)) || (new_state[0] == Real(0))) {
	  //std::cout << "Converged" << std::endl;
	  //We've converged
	  x = new_x;
	  curr_state = new_state;
	  return 2;
	}
	
	//Check if the function has increased
	if (std::abs(new_state[0]) > std::abs(curr_state[0])) {
	  //std::cout << "Function increased!" << std::endl;
	  //The function increased! method has failed
	  return 0;
	}

	//Not converged or failed, update bounds and continue
	if (delta >= 0)
	  low_bound = x;

	if (delta <= 0)
	  high_bound = x;

	curr_state = new_state;
	x = new_x;

	//std::cout << "new bounds [" << low_bound << "," << high_bound << "]" << std::endl;

	return 1;
      }
    }
      
    /*! \brief A single step of the Newton-Raphson method for finding roots.
     */
    template<class F, size_t Derivatives, class Real>
    inline int newton_raphson_step(const F& f, std::array<Real, Derivatives>& curr_state,  
				    Real& x, Real& low_bound, Real& high_bound, Real x_precision) {
      static_assert(Derivatives >= 1, "Can only perform newton raphson if at least 1 derivative is available");

      //std::cout << "NR Method" << std::endl;
      if (curr_state[1] == 0) {
	//std::cout << "Zero derivative!" << std::endl;
	//Zero derivatives cause x to diverge, so abort
	return 0;
      }

      return detail::process_iterative_step(f, curr_state, x, x - curr_state[0] / curr_state[1], low_bound, high_bound, x_precision);
    }
    
    /*! \brief A single step of Halley's method for finding roots.
     */
    template<class F, size_t Derivatives, class Real>
    inline int halley_step(const F& f, std::array<Real, Derivatives>& curr_state,  
			    Real& x, Real& low_bound, Real& high_bound, Real x_precision) {
      static_assert(Derivatives >= 2, "Can only perform Halley iteration if at least 1 derivative is available");

      //std::cout << "Halley Method" << std::endl;
      Real numerator = 2 * curr_state[0] * curr_state[1];
      Real denominator = 2 * curr_state[1] * curr_state[1] - curr_state[0] * curr_state[2];
      
      if ((denominator == 0) || !std::isfinite(denominator)) {
	//std::cout << "Halley Overflow!" << std::endl;
	//Cannot proceed with a zero or infinite denominator, or if the div
	//has overflowed (+inf)
	return 0;
      }
      
      Real delta = - numerator / denominator;
      Real deltaNR = - curr_state[0] / curr_state[1];
      //std::cout << "delta = " << delta << std::endl;
      //std::cout << "deltaNR = " << deltaNR << std::endl;
      if (std::signbit(delta) != std::signbit(deltaNR)) {
	//std::cout << "Halley != NR" << std::endl;
	//The Halley and Newton Raphson iterations would proceed in
	//opposite directions. This happens near multiple roots where
	//the second derivative causes overcompensation. Fail so NR is
	//used instead.
	return 0;
      }

      return detail::process_iterative_step(f, curr_state, x, x + delta, low_bound, high_bound, x_precision);
    }

    /*! \brief A single step of Schroeder's method for finding roots.
     */
    template<class F, size_t Derivatives, class Real>
    inline int schroeder_step(const F& f, std::array<Real, Derivatives>& curr_state,  Real& x, Real& low_bound, Real& high_bound, Real x_precision) {
      static_assert(Derivatives >= 2, "Can only perform Halley iteration if at least 1 derivative is available");

      if (curr_state[1] == 0)
	//Cannot proceed with a zero first derivative
	return 0;

      return detail::process_iterative_step(f, curr_state, x, x - 2 * curr_state[0] * curr_state[1] - curr_state[2] * curr_state[0] * curr_state[0] / (2 * curr_state[1] * curr_state[1] * curr_state[1]), low_bound, high_bound, x_precision);
    }

    /*! \brief A single bisection step.
     */
    template<class F, size_t Derivatives, class Real>
    inline int bisection_step(const F& f, std::array<Real, Derivatives>& curr_state, Real& x, Real& low_bound, Real& high_bound, Real x_precision) {
      if ((low_bound >= high_bound) || std::isinf(low_bound) || std::isinf(high_bound))
	//This is not a valid interval
	return 0;

      //std::cout << "Bisection" << std::endl;
      auto f_low = f(low_bound);
      auto f_high = f(high_bound);
      
      if (std::signbit(f_low[0]) == std::signbit(f_high[0])) {
	//std::cout << "No sign change!" << std::endl;
	//No sign change in the interval
	return 0;
      }

      Real x_mid = (low_bound + high_bound) / 2;
      auto new_state = f(x_mid);

      //std::cout << "new_x = " << x_mid << std::endl;

      Real delta = x_mid - x;
      
      if ((std::abs(delta) < std::abs(x_precision * x_mid)) || (new_state[0] == Real(0))) {
	//std::cout << "Converged!" << std::endl;
	//We've converged
	x = x_mid;
	curr_state = new_state;
	return 2;
      }

      //Update bounds and continue
      if (std::signbit(new_state[0]) == std::signbit(f_high[0])) {
	high_bound = x_mid;
      } else {
	low_bound = x_mid;
      }
      
      //std::cout << "new bounds [" << low_bound << "," << high_bound << "]" << std::endl;
      x = x_mid;
      curr_state = new_state;
      
      return 1;
    }
    
    /*! \brief Safeguarded newton_raphson method for detecting a root.

      This returns false if the method is not converging or if the
      number of iterations was exceeded.
     */
    template<class F, class Real>
    bool newton_raphson(const F& f, Real& x, size_t iterations = 20, 
			Real low_bound = -HUGE_VAL, Real high_bound = +HUGE_VAL, 
			int digits = std::numeric_limits<Real>::digits / 2)
    {
      auto last_state = f(x);
      static_assert(last_state.size() > 1, "Require one derivative of the objective function for Halley's method");

      Real x_precision = static_cast<Real>(ldexp(1.0, 1 - digits));

      if (last_state[0] == 0)
	return true;
    
      int status = 1;
      while ((--iterations) && (status != 2)) {
	status = newton_raphson_step(f, last_state, x, low_bound, high_bound, x_precision);
	if (!status) {
	  status = bisection_step(f, last_state, x, low_bound, high_bound, x_precision);
	  if (!status)
	    return false;
	}
      }

      return (status == 2);
    }

    /*! \brief Safeguarded Halley's method for detecting a root.

      This returns false if the method is not converging or if the
      number of iterations was exceeded.
     */
    template<class F, class Real>
    bool halleys_method(const F& f, Real& x, size_t iterations = 20, Real low_bound = -HUGE_VAL, 
			Real high_bound = +HUGE_VAL, int digits = std::numeric_limits<Real>::digits / 2)
    {
      auto last_state = f(x);

      static_assert(last_state.size() > 2, "Require two derivatives of the objective function for Halley's method");

      Real x_precision = static_cast<Real>(ldexp(1.0, 1 - digits));

      if (last_state[0] == 0)
	return true;

      //std::cout << "x0=" << x << std::endl;
      int status = 1;
      while ((--iterations) && (status != 2)) {
	status = halley_step(f, last_state, x, low_bound, high_bound, x_precision);
	if (!status) {
	  status = newton_raphson_step(f, last_state, x, low_bound, high_bound, x_precision);
	  if (!status) {
	    status = bisection_step(f, last_state, x, low_bound, high_bound, x_precision);
	    if (!status)
	      return false;
	  }
	}
      }
      return (status == 2);
    }
  }
}
