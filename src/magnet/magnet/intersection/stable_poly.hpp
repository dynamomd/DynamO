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
      if (f[2] == 0) return nextEvent(::magnet::math::Polynomial<1>(f));
      
      const double arg = f[1] * f[1] - 2 * f[2] * f[0];

      if (f[2] < 0) {
	//Polynomial limits towards overlap at t -> +inf
	
	//If there are no roots, it never escapes overlap, return the
	//time of the turning point or now if it is in the past
	if (arg <= 0) return std::max(0.0, -f[1] / f[2]);

	//There are roots. select a stable form of the quadratic to
	//compute the largest root.
	if (f[1] > 0)
	  return std::max(0.0, (-f[1] - std::sqrt(arg)) / f[2]);
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
      if (f[3] == 0) return nextEvent(::magnet::math::Polynomial<2>(f));
      
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
	if (roots.size() != 1)
	  M_throw() << "Unexpected solution in cubic function";
	return std::max(0.0, roots[0]);
      }
      
      //There must be one or more turning points.
      if ((derivroots[0] > 0) && (roots[0] < derivroots[0]))
	return std::max(0.0, roots[0]);
      return std::max(0.0, std::max(derivroots[1], roots[roots.size()-1]));
    }
    
    inline double nextEvent(const ::magnet::math::Polynomial<4>& f, double f0char, double precision=1e-16)
    {
      //If the polynomial is cubic, drop to that solution
      if (f[4] == 0) return nextEvent(::magnet::math::Polynomial<3>(f));

      const double rootthreshold = f0char * precision;
      
      //Determine and sort the roots of the derivative
      const auto df = derivative(f);
      auto droots = solve_roots(df);
      std::sort(droots.begin(), droots.end());
      
      if (droots.size() == 0)
	M_throw() << "Unexpected, cubic with no roots?";

      double last_point = 0;
      double last_f = eval(f, last_point);
      double last_df = eval(df, last_point);

      //Check if the start time satisfies the conditions
      if ((last_f <= 0) && (last_df < 0)) {
	return 0;
      }

      auto bisectFunc = [&] (double t) { return eval(f, t); };
	
      //Look through the roots of the derivative to see if we can
      //bracket a root, or identify a turning point
      for (double droot : droots) 
	if (droot > last_point) {
	  const double new_f = eval(f, droot);
	  const double new_df = eval(df, droot);
	  
	  if (last_f >= 0) {
	    //Looking for a root
	    if (new_f < 0)
	      //Root bisected between [last_point, droot]
	      return magnet::math::bisect(bisectFunc, last_point, droot, rootthreshold);
	    //There is no root within the bracket, update the location
	    //and carry on.
	  } else {
	    //We're looking to see if the curve escapes before turning
	    //back on itself
	    if (new_f < 0)
	      //Its curved back while overlapped. This is the start of
	      //an unstable condition
	      return new_f;
	  }
	  
	  //Update our current location
	  last_point = droot;
	  last_f = new_f;
	  last_df = new_df;
	}

      
      //We've checked all derivative roots to try to bracket any
      //roots. Only thing left to check
      if (f[4] < 0) {	
	double delta_t = std::pow(- 24 * f0char / f[4], 0.25);
	last_point += delta_t;
	//Move forward exponentially, looking to bracket the root
	while (eval(f, last_point) >= 0)
	  last_point += (delta_t *= 2);
	return magnet::math::bisect(bisectFunc, last_point - delta_t, last_point, rootthreshold);
      }
      return HUGE_VAL;
    }
  }
}
