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
	template<size_t d> double eval() const { return Base::template eval<d+derivative>(); }
	template<size_t d> double max() const { return Base::template max<d+derivative>(); }
	bool test_root() const { return true; }
      };

      template<class F> std::pair<bool, double>
      halleysMethod(const F& f, double t_guess, const double t_min, const double t_max, size_t iterations = 100)
      {
	double delta = 0;
	do
	  {
	    const double f0 = f.template eval<0>(t_guess);
	    //Check if we're at a root already
	    if (f0==0) return std::pair<bool,double>(true, t_guess); 
	    //Begin method
	    const double f1 = f.template eval<1>(t_guess);
	    const double f2 = f.template eval<2>(t_guess);
	    
	    //If we have zero derivatives, just abort as the algorithm
	    //will try again from somewhere else
	    if ((f1 == 0) && (f2==0)) break;

	    if (f2 == 0) //Use a newton step
	      delta = f0 / f1;
	    else
	      {
		const double denom = 2 * f0;
		const double num = 2 * f1 - f0 * (f2 / f1);
		if ((std::abs(num) < 1) && (std::abs(denom) >= std::abs(num) * std::
	      }
	    else
	    delta = ;
	    
	  }
	while (--iterations);
	
	//Failed! We're not at a root, and the best lower bound we can do is t_min.
	return std::pair<bool,double>(false, t_min);
      }

      template<class F> std::pair<bool, double> nextDecreasingRoot(F f, double t_min, double t_max, const double tol)
      {
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
