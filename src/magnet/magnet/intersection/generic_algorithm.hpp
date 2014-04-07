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
#include <magnet/math/frenkelroot.hpp>
#include <boost/math/tools/roots.hpp>

namespace magnet {
  namespace intersection {
    namespace detail {
      /*! \brief A class which takes the derivative of an overlap
	function.
       */
      template <class Base, int derivative = 1> class FDerivative: public Base
      {
      public:
	FDerivative(const Base& b):Base(b) {}
	template<size_t d> double eval(double dt) const { return Base::template eval<d+derivative>(dt); }
	template<size_t d> double max() const { return Base::template max<d+derivative>(); }
      };

      template<class F> std::pair<bool, double>
      halleySearch(const F& f, double t_guess, const double t_min, const double t_max, const int digits = 12, size_t iterations = 50)
      {
	const double digitfactor = std::ldexp(1.0, 1 - digits);
	do {
	  const double f0 = f.template eval<0>(t_guess);
	  //Check if we're at a root already
	  if (f0 == 0) return std::pair<bool,double>(true, t_guess); 
	  //Begin method
	  const double f1 = f.template eval<1>(t_guess);
	  const double f2 = f.template eval<2>(t_guess);
	    
	  //If we have zero derivatives, just abort as the algorithm
	  //will try again from somewhere else
	  if ((f1 == 0) && (f2==0)) break;

	  double delta = 0;
	  if (f2 != 0) 
	    {
	      const double denom = 2 * f0;
	      const double num = 2 * f1 - f0 * (f2 / f1);
	      if ((std::abs(num) < 1) && (std::abs(denom) >= std::abs(num) * std::numeric_limits<double>::max()))
		//Possible overflow, switch to newton step
		delta = f0 / f1;
	      else
		delta = denom / num;
	      
	      if (delta * f1 / f0 < 0)
		//Possible cancellation error, switch to newton step
		delta = f0 / f1;
	    }
	  else
	    delta = f0 / f1;
	  
	  t_guess -= delta;
	  
	  //Check we've not gone out of range
	  if ((t_guess < t_min) || (t_guess > t_max)) break;

	  //Check if we've converged
	  if (std::abs(t_guess * digitfactor) >= std::abs(delta))
	    return std::pair<bool,double>(true, t_guess);
	} while (--iterations);
	
	//Failed! We're not at a root, and the best lower bound we can do is t_min.
	return std::pair<bool,double>(false, HUGE_VAL);
      }
      

      template<class F> std::pair<bool, double> nextDecreasingRoot(const F& f, double t_min, double t_max, const int digits=15)
      {
	enum {
	  LOW = 0,
	  HIGH = 1,
	};

	int active_boundary = LOW;
	const double f2max = f.template max<2>();
	
	//Loop while we still have a valid search window
	while (t_min < t_max)
	  {
	    double& t_current = (active_boundary == HIGH) ? t_max : t_min;
	    const double f0 = f.template eval<0>(t_current);
	    const double f1 = f.template eval<1>(t_current);
	    
	    //Improve the boundary first
	    auto boundary_roots = math::quadraticEquation(- 0.5 * std::copysign(f2max, f0), f1, f0);
	    if (active_boundary == LOW)
	      t_min += std::max(boundary_roots.first, boundary_roots.second);
	    else /*(active_boundary == HIGH)*/
	      t_max += std::min(boundary_roots.first, boundary_roots.second);
	    
	    //Switch the working boundary for the next iteration
	    if (t_max != HUGE_VAL)
	      active_boundary = !active_boundary;

	    //Now search for a root:
	    auto search_result = halleySearch(f, t_current, t_min, t_max);
	    
	    //If searching failed, restart from the other bound
	    if (!search_result.first) continue;
	    
	    {
	      //Searching succeeded check for earlier roots in this
	      //interval using recursion.
	      const double current_root = search_result.second;
	      const double f1 = f.template eval<1>(current_root);
	      const double inner_t_max = current_root - 2.0 * std::abs(f1 / f2max);
	      auto check_result = nextDecreasingRoot(f, t_min, inner_t_max, digits);
	      if (check_result.second != HUGE_VAL)
		//We have found an earlier root. Return this one instead.
		return check_result;
		
	      //There are no earlier roots, current_root is the
	      //next root in the region. Return it if its
	      //approaching
	      if (f.template eval<1>(current_root) < 0)
		return std::pair<bool,double>(true, current_root);
	      
	      //Use it as a new lower bound if its a receeding
	      //root and carry on
	      t_min = current_root + 2.0 * std::abs(f1 / f2max);
	    }
	  }
	return std::pair<bool,double>(false, HUGE_VAL);
      }
    }
      
    /*! \brief A generic implementation of the stable EDMD algorithm
      which uses Frenkel's root finder on overlap functions.
      
      \tparam T The type of the overlap function which is being solved.
      \param f The overlap function.
      \param tol The maximum error on the frenkel root finder
    */
    template<class T> std::pair<bool, double> nextEvent(const T& f, const double t_min, const double t_max, const double tol)
    {
      const double f0 = f.template eval<0>(t_min);
      const double f1 = f.template eval<1>(t_min);

      //If particles are not in contact, just search for the next contact
      if (f0 > 0) return detail::nextDecreasingRoot(f, t_min, t_max, tol);
      
      //Particles are either in contact or overlapped. Check if they're approaching
      if (f1 < 0) return std::pair<bool, double>(true, 0.0);

      //Overlapping but they're moving away from each other. Determine
      //when they reach their next maximum separation (it may be the
      //current time if f1==0).
      detail::FDerivative<T> fprime(f);
      std::pair<bool, double> derivroot = detail::nextDecreasingRoot(fprime, t_min, t_max, tol);      
      //Check if they just keep retreating from each other, which means that they never interact
      if (derivroot.second == HUGE_VAL) return std::pair<bool, double>(false, HUGE_VAL);
      
      //If they are still overlapping at this time, it doesn't matter
      //if derivroot is a virtual (recalculate) event or an actual
      //turning point. We can just return it and either have a
      //collision then or recalculate then.
      if (f.template eval<0>(derivroot.second) < 0) return derivroot;
      
      //Real or virtual, the derivroot contains a time before the
      //next interaction which is outside the invalid state, we just
      //use this as our lower bound and carry on the search.
      return detail::nextDecreasingRoot(f, derivroot.second, t_max, tol);
    }
  }
}
