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
      */
      template<class F, size_t Derivatives, class Real>
      inline bool process_iterative_step(const F& f, std::array<Real, Derivatives>& curr_state,
					 Real& x, Real new_x, Real& low_bound, Real& high_bound) {
	if ((new_x >= high_bound) || (new_x <= low_bound))
	  //Out of bounds step, abort
	  return false;
	
	//Re evalulate the derivatives
	auto new_state = f(new_x);
	
	if (std::abs(new_state[0]) >= std::abs(curr_state[0]))
	  //The function did not decrease, abort 
	  return false;
      
	//Accept this step
	low_bound = std::min(low_bound, x);
	high_bound = std::max(high_bound, x);
	curr_state = new_state;
	x = new_x;

	return true;
      }
    }
      
    /*! \brief A single step of the Newton-Raphson method for finding roots.
     */
    template<class F, size_t Derivatives, class Real>
    inline bool newton_raphson_step(const F& f, std::array<Real, Derivatives>& curr_state,  
				    Real& x, Real& low_bound, Real& high_bound) {
      static_assert(Derivatives >= 1, "Can only perform newton raphson if at least 1 derivative is available");
      if (curr_state[1] == 0)
	//Zero derivatives cause x to diverge, so abort
	return false;

      return detail::process_iterative_step(f, curr_state, x, x - curr_state[0] / curr_state[1], low_bound, high_bound);
    }
    
    /*! \brief A single step of Halley's method for finding roots.
     */
    template<class F, size_t Derivatives, class Real>
    inline bool halley_step(const F& f, std::array<Real, Derivatives>& curr_state,  
			    Real& x, Real& low_bound, Real& high_bound) {
      static_assert(Derivatives >= 2, "Can only perform Halley iteration if at least 1 derivative is available");

      Real denominator = 2 * curr_state[1] * curr_state[1] - curr_state[0] * curr_state[2];
      
      if ((denominator == 0) || !std::isfinite(denominator))
	//Cannot proceed with a zero denominator, or if it has overflowed (+inf)
	return false;
      
      //Take a step
      return detail::process_iterative_step(f, curr_state, x, x - 2 * curr_state[0] * curr_state[1] / denominator, low_bound, high_bound);
    }

    /*! \brief A single step of Schroeder's method for finding roots.
     */
    template<class F, size_t Derivatives, class Real>
    inline bool schroeder_step(const F& f, std::array<Real, Derivatives>& curr_state,  Real& x, Real& low_bound, Real& high_bound) {
      static_assert(Derivatives >= 2, "Can only perform Halley iteration if at least 1 derivative is available");

      if (curr_state[1] == 0)
	//Cannot proceed with a zero first derivative
	return false;

      //Take a step
      return detail::process_iterative_step(f, curr_state, x, x - 2 * curr_state[0] * curr_state[1] - curr_state[2] * curr_state[0] * curr_state[0] / (2 * curr_state[1] * curr_state[1] * curr_state[1]), low_bound, high_bound);
    }

    /*! \brief A single bisection step.
     */
    template<class F, size_t Derivatives, class Real>
    inline bool bisection_step(const F& f, std::array<Real, Derivatives>& curr_state, Real& low_bound, Real& high_bound) {
      if ((low_bound >= high_bound) || std::isinf(low_bound) || std::isinf(high_bound))
	//This is not a valid interval
	return false;

      auto f_low = f(low_bound);
      auto f_high = f(high_bound);
      
      if (std::signbit(f_low[0]) == std::signbit(f_high[0]))
	//No sign change in the interval
	return false;
      
      Real x_mid = (low_bound + high_bound) / 2;
      auto last_state = f(x_mid);
      
      if (std::signbit(last_state[0]) == std::signbit(f_high[0])) {
	high_bound = x_mid;
      } else {
	low_bound = x_mid;
      }
      
      return true;
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

      Real factor = static_cast<Real>(ldexp(1.0, 1 - digits));
   
      //Bound where the function can go
      if (last_state[1] > 0)
	high_bound = std::min(high_bound, x);
      if (last_state[1] < 0)
	low_bound = std::max(low_bound, x);
      
      Real old_x = std::numeric_limits<Real>::max();

      //std::cout << "x0= " << x << std::endl;
      while (--iterations) {
	//std::cout << "it=" << iterations << std::endl;
	if (!newton_raphson_step(f, last_state, x, low_bound, high_bound)) {
	  //std::cout << "NR fail!" << std::endl;
	  //Newton Raphson failed, try a bisection step
	  if (bisection_step(f, last_state, low_bound, high_bound))
	    x = (low_bound + high_bound) / 2;
	  else
	    return false;
	}

	Real delta = x - old_x;
	old_x = x;

	//std::cout << "x= " << x << std::endl;
	//std::cout << "precision=" << std::abs(x * factor) << std::endl;
	//std::cout << "delta=" << std::abs(delta) << std::endl;	

	if ((last_state[0] == Real(0)) || (std::abs(x * factor) > std::abs(delta)))
	  return true;
      }
      
      //std::cout << "Too many iterations" << std::endl;
      //Convergence failure
      return false;
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

      Real factor = static_cast<Real>(ldexp(1.0, 1 - digits));

      if (last_state[0] == 0)
	return true;

      if (last_state[1] != 0) {
	Real direction = - last_state[0] / last_state[1];
	//Bound where the function can go
	if (direction < 0)
	  high_bound = std::min(high_bound, x);
	else
	  low_bound = std::max(low_bound, x);
      }
      
      Real old_x = std::numeric_limits<Real>::max();

      //std::cout << "x0= " << x << std::endl;
      while (--iterations) {	
	//std::cout << "it=" << iterations << std::endl;
	if (!halley_step(f, last_state, x, low_bound, high_bound)) {
	  //std::cout << "H fail!" << std::endl;
	  //Halley step failed, try a Newton-Raphson step
	  if (!newton_raphson_step(f, last_state, x, low_bound, high_bound)) {
	    //std::cout << "NR fail!" << std::endl;
	    //Newton-Raphson failed! Try a bisection
	    if (bisection_step(f, last_state, low_bound, high_bound))
	      x = (low_bound + high_bound) / 2;
	    else
	      return false;
	  }
	}
	
	Real delta = x - old_x;
	old_x = x;

	//std::cout << "x= " << x << std::endl;
	//std::cout << "precision=" << std::abs(x * factor) << std::endl;
	//std::cout << "delta=" << std::abs(delta) << std::endl;	

	if ((last_state[0] == Real(0)) || (std::abs(x * factor) > std::abs(delta)))
	  return true;
      }
      
      //std::cout << "Too many iterations" << std::endl;
      return false;
    }
  }
}
