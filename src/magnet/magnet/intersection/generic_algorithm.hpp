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

namespace magnet {
  namespace intersection {
    namespace detail {
      /*! \brief A class which takes the derivative of an overlap
	function.
       */
      template <class Base, int derivative = 1> class OFDerivative: public Base
      {
      public:
	OFDerivative(const Base& b):Base(b) {}
	template<size_t d> double eval() const { return Base::template eval<d+derivative>(); }
	template<size_t d> double max() const { return Base::template max<d+derivative>(); }
	bool test_root() const { return true; }
      };
    }

    /*! \brief A generic implementation of the stable EDMD algorithm
      which uses Frenkel's root finder on overlap functions.
      
      \tparam T The type of the overlap function which is being solved.
      \param f The overlap function.
      \param err The maximum error on the frenkel root finder
    */
    template<class T> std::pair<bool, double> generic_algorithm(T f, double t_max, double err)
    {
      double f0 = f.template eval<0>();
      double f1 = f.template eval<1>();

      //First treat overlapping or in contact particles which are approaching
      if ((f0 <= 0) && (f1 < 0)) return std::pair<bool, double>(true, 0.0);
    
      //Now treat overlapping particles which are not approaching
      if (f0 < 0)
	{
	  //Overlapping but they're moving away from each
	  //other. Determine when they reach their next maximum
	  //separation.
	  detail::OFDerivative<T> fprime(f);

	  std::pair<bool, double> derivroot = math::frenkelRootSearch(fprime, 0, t_max, err);

	  //Check if they just keep retreating from each other, which means that they never interact
	  if (derivroot.second == HUGE_VAL) return std::pair<bool, double>(false, HUGE_VAL);

	  //Check if the time returned is not overlapping
	  T froot(f);
	  froot.stream(derivroot.second);

	  //If they are still overlapping at this time, it doesn't
	  //matter if derivroot is a virtual (recalculate) event or
	  //an actual turning point. We can just return it
	  if (froot.template eval<0>() < 0) return derivroot;

	  //Real or virtual, the derivroot contains a time before the
	  //next interaction which is outside the invalid state, we just
	  //use this as our lower bound
	  return math::frenkelRootSearch(f, derivroot.second, t_max, err);
	}

      //If the particles are in contact, but not approaching, we need to
      //skip this initial root
      double t_min = (f0 == 0) ? 2.0 * std::abs(f.template eval<1>()) / f.template max<2>() : 0;
      return math::frenkelRootSearch(f, t_min, t_max, err);
    }
  }
}
