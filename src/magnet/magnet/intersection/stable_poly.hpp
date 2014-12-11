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
#include <magnet/math/polynomial.hpp>
#include <magnet/math/bisect.hpp>

namespace magnet {
  namespace intersection {

    /*!  \brief Implementation of the generic stable event-detection
      algorithm.

      For this generic implementation to work, the function f must
      have implementations of the functions \ref
      magnet::math::derivative, \ref magnet::math::shift_function, and
      \ref magnet::math::next_root defined.
    */
    template<class Function>
    double nextEvent(const Function& f)
    {
      //Determine the derivative
      const auto df = ::magnet::math::derivative(f);
      
      const double f_start = eval(f, 0.0);
      const double df_start = eval(df, 0.0);

      double f_shift = 0.0;
      //Check if starting overlapped
      if (f_start <= 0) {

	//If we're approaching then the current time is the time of
	//the next event
	if (df_start < 0)
	  return 0;

	//We need to find when the derivative next turns
	//negative. Here we recurse, as this will correctly check that
	//the function indeed turns negative.
	const double next_df_root = nextEvent(df);

	//If the overlap function never becomes negative, then there
	//is never an event.
	if (next_df_root == HUGE_VAL)
	  return HUGE_VAL;

	//If it turns around while still overlapped/in contact, then
	//the turning point is the next event.
	if (eval(f, next_df_root) <= 0)
	  return next_df_root;

	//The turning point of the derivative is in the positive
	//region, so use this as the starting point of a normal event
	//search.
	f = ::magnet::math::shift_function(f, next_df_root);
	f_shift = next_df_root;
      }
      
      return f_shift + ::magnet::math::next_root(f);
    }

    //Use a template overload
    template<size_t Order>
    double nextEvent(const ::magnet::math::Polynomial<Order>& f) {
      //If the polynomial is lower order, drop to that solution (it
      //may be solvable by radicals)
      if (f[Order] == 0)
	return nextEvent(::magnet::math::change_order<Order-1>(f));

      //Call the generic template implementation
      nextEvent<::magnet::math::Polynomial<Order> >(f);
    }
    
    /*! \brief Calculate the interval until the 1st order Polynomial
        is negative and has a negative derivative. 
	
	\param f The Polynomial under consideration.
    */
    inline double nextEvent(const ::magnet::math::Polynomial<1>& f) {
      //If the gradient is not negative now, it never will be
      if (f[1] >= 0) return HUGE_VAL;
      //Return the time of the root, or now if we're past the root
      return std::max(0.0, - f[0] / f[1]);
    }
    
    /*! \brief Calculate the interval until the 2nd order Polynomial
        is negative and has a negative derivative.
	
	\param f The Polynomial under consideration.
    */
    inline double nextEvent(const ::magnet::math::Polynomial<2>& f) {
      //If the polynomial is linear, drop to that solution
      if (f[2] == 0) return nextEvent(::magnet::math::change_order<1>(f));
      
      const double arg = f[1] * f[1] - 4 * f[2] * f[0];

      if (f[2] < 0) {
	//Polynomial limits towards overlap at t -> +inf
	
	//If there are no roots, it never escapes overlap, return the
	//time of the turning point or now if it is in the past
	if (arg <= 0) return std::max(0.0, -f[1] / (2 * f[2]));

	//There are roots. select a stable form of the quadratic to
	//compute the largest root.
	if (f[1] > 0)
	  return std::max(0.0, (-f[1] - std::sqrt(arg)) / (2 * f[2]));
	else
	  return std::max(0.0, 2 * f[0] / (-f[1] + std::sqrt(arg)));
      } 
      /* else (f[2] > 0)*/

      //Interactions only happen if there are roots and we're in the
      //region between the first root and the turning point
      if ((f[1] >= 0) || arg <= 0) return HUGE_VAL;

      //Return the time of the root using a stable quadratic formula
      return std::max(0.0, 2 * f[0] / (-f[1] + std::sqrt(arg)));
    }

    /*! \brief Calculate the interval until the 3rd order Polynomial
        is negative and has a negative derivative.
	
	\param f The Polynomial under consideration.
    */
    inline double nextEvent(const ::magnet::math::Polynomial<3>& f)
    {
      //If the polynomial is quadratic, drop to that solution
      if (f[3] == 0) return nextEvent(::magnet::math::change_order<2>(f));
      
      //Calculate and sort the roots of the overlap function
      auto roots = solve_roots(f);
      std::sort(roots.begin(), roots.end());
      
      //Calculate and sort the roots of the overlap function's
      //derivative
      auto derivroots = solve_roots(derivative(f));
      std::sort(derivroots.begin(), derivroots.end());
      
      if (f[3] > 0) {
	//Cubic limits away from overlap at t -> +inf.
	
	//The only region where a cubic overlap function with
	//positive f[3] is decreasing is between the two turning
	//points (if they exist).
	
	//Check there are two turning points, and that we are before
	//the second one.
	if (derivroots.size() != 2) return HUGE_VAL;
	if (derivroots[1] < 0) return HUGE_VAL;
	
	//If there is only one root, it must come after the first
	//turning point (otherwise the two turning points are
	//outside the overlapped (negative f) zone.
	if ((roots.size()==1) && (roots[0] < derivroots[0])) return HUGE_VAL;
	
	//There is definitely an event. 
	
	if ((roots.size() == 1) || (roots.size() == 2))
	  //If there's one root, the event happens after the first
	  //turning point. If there's two roots, the first turning
	  //point must be a point of inflection, and it still
	  //happens after the first turning point.
	  return std::max(0.0, derivroots[0]);
	
	//If there are three roots it happens after the first event.
	return std::max(0.0, roots[1]);
      } 
      
      /* else (f[3] < 0)*/ 
      //Cubic limits towards overlap at t -> +inf
      
      if (derivroots.size() == 0) {
	//If there are no turning points, then the function is always
	//decreasing, so it must have one root and the collision
	//happens after this point.
	return std::max(0.0, roots[0]);
      }
      
      //There must be one or more turning points.
      if ((derivroots[0] > 0) && (roots[0] < derivroots[0]))
	return std::max(0.0, roots[0]);
      return std::max(0.0, std::max(derivroots[1], roots[roots.size()-1]));
    }
  }
}
