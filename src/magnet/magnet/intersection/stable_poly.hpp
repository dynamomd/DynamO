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
    /*! \brief Implementation of \ref next_root which utilises \ref
      solve_roots.
      
      This is mainly useful for simple functions with fast root
      detection algorithms. For example, the linear, quadratic, and
      cubic functions are solved via radicals, so this function is
      used to provide an implementation of \ref next_root in these
      cases.
    */
    template<class Real, size_t Order, char Letter>
    typename std::enable_if<(Order < 5), Real>::type 
    next_negative(const ::magnet::math::Polynomial<Order, Real, Letter>& f) {

      double first_root = HUGE_VAL;
      double second_root = HUGE_VAL;
      auto roots = magnet::math::solve_roots(f);
      auto it = roots.begin();

      //skip negative roots
      while ((it != roots.end()) && (*it < 0)) ++it;
      
      //Grab the next two positive roots (if they exist)
      if (it != roots.end()) {
	first_root = *it;
	if ((++it) != roots.end())
	  second_root = *it;
      }

      //There are two numerical cases to be avoided. 
      //
      //First, although f(0)>0 on entry to this function, there may be
      //a sign change numerically at x=0 due to finite precision. The
      //only way to test for this is to actually test if an unexpected
      //sign change has happened sometime before the next root (or
      //infinity if no root has happened). If there is no root change,
      //then x=1 is a logical point to sample the sign of the
      //root. Otherwise we must sample in-between now and the next
      //root.
      const double sample_point = (first_root == HUGE_VAL) ? 1.0 : first_root / 2;
      if (magnet::math::eval(f, magnet::math::Variable<Letter>() == sample_point) < 0)
	//They're already approaching, return an instant collision
	return 0.0;

      //Second, the root detected in the future may have an odd or
      //even multiplicity. Odd roots are actual sign transitions
      //whereas even roots are not. Therefore it is not enough to
      //simply detect the root, we must test that it actually changes
      //sign. Numerically determining the multiplicity of a root of a
      //floating point polynomial (or other function) is futile. We
      //only care if a sign change occurs.
      
      if (first_root != HUGE_VAL)
	do {
	  const double sample_point = (second_root == HUGE_VAL) ? first_root + 1.0 : (first_root + second_root) / 2;
	  if (magnet::math::eval(f, magnet::math::Variable<Letter>() == sample_point) < 0)
	    return first_root;
	  first_root = second_root;
	  second_root = ((++it) != roots.end()) ? *it : HUGE_VAL;
	} while (first_root != HUGE_VAL);
      
      return HUGE_VAL;
    }
      
    //template<class Real, size_t Order, char Letter>
    //typename std::enable_if<(Order > 3), Real>::type 
    //next_negative(const Polynomial<Order, Real, Letter>& f) {
    //  //Drop down to a lower order solver if available
    //  if (f[Order] == 0)
    //	return next_root(change_order<Order-1>(f));
    //  
    //  //Test if we're currently turning negative
    //  if (f[0] < 0)
    //	return 0.0;
    //
    //  auto pos_root_bounds = VAS_real_root_bounds(f);
    //  if (pos_root_bounds.size() != 0) {
    //	std::sort(pos_root_bounds.begin(), pos_root_bounds.end());
    //	const Real& a = pos_root_bounds[0].first;
    //	const Real& b = pos_root_bounds[0].second;
    //	boost::uintmax_t iter = 100;
    //	auto root = boost::math::tools::toms748_solve([&](Real x) { return eval(f, x); }, a, b, boost::math::tools::eps_tolerance<Real>(100), iter);
    //	return (root.first + root.second) / 2;
    //  }
    //
    //  return std::numeric_limits<Real>::infinity();
    //}

    /*!  \brief Implementation of the generic stable event-detection
      algorithm.

      For this generic implementation to work, the function f must
      have implementations of the functions \ref
      magnet::math::derivative, \ref magnet::math::shift_function, and
      \ref magnet::math::next_root defined.
    */
    template<class Function>
    double nextEvent(Function f)
    {
      ::magnet::math::Variable<'t'> t;

      //Determine the derivative
      const auto df = ::magnet::math::derivative(f, t);
      
      const double f_start = eval(f, t == 0.0);
      const double df_start = eval(df, t == 0.0);

      //Check if starting was not overlapped.
      if (f_start >= 0) {
	return next_negative(f);
      }
      
      //If we're approaching then the current time is the time of
      //the next event
      if (df_start < 0)
	return 0;
      
      //We need to find when the derivative next turns negative.
      const double df_next_root = next_negative(df);
      
      //If the overlap function never becomes negative, then there is
      //never an event.
      if (df_next_root == HUGE_VAL)
	return HUGE_VAL;
      
      //If it turns around while still overlapped/in contact, then
      //the turning point is the next event. 
      if (eval(f, t==df_next_root) <= 0)
	return df_next_root;
      
      //Ok, try searching after the maxima for events
      auto shift_f = shift_function(f, df_next_root);
      return df_next_root + next_negative(shift_f);
    }

    /*! \brief Calculate the interval until the 1st order Polynomial
        is negative and has a negative derivative.

	This is a specialisation for linear polynomials to close the
	recursive definition of the general nextEvent implementation.
	
	\param f The Polynomial under consideration.
    */
    template<char Var>
    inline double nextEvent(const ::magnet::math::Polynomial<1, double, Var>& f) {
      //If the gradient is not negative now, it never will be
      if (f[1] >= 0) return HUGE_VAL;
      //Return the time of the root, or now if we're past the root
      return std::max(0.0, - f[0] / f[1]);
    }
    
    /*! \brief Calculate the interval until the 2nd order Polynomial
        is negative and has a negative derivative.
	
	This is an optimised case, as many interactions use quadratic
	overlap functions.
	
	\param f The Polynomial under consideration.
    */
    template<char Var>
    inline double nextEvent(const ::magnet::math::Polynomial<2, double, Var>& f) {
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
  }
}
