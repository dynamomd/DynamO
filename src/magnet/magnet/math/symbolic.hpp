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
#include <complex>

namespace magnet {
  namespace math {
    /*! \brief Returns the empty product of a type.
      
      The empty product is a term whose multiplicative action is null
      (can be ignored).
    */
    template<class T>
    constexpr T empty_product(const T&) { 
      static_assert(std::is_arithmetic<T>::value, "Nullary products must be defined for this type");
      return T(1); 
    }

    /*! \brief Returns the empty product of a std::complex type.
      
      The empty product is a term whose multiplicative action is null
      (can be ignored).
    */
    template<class T>
    constexpr std::complex<T> empty_product(const std::complex<T>&) { 
      return std::complex<T>(T(1), T()); 
    }

    /*! \brief Returns the empty sum of a type.
      
      The empty sum is a term whose additive (and typically its
      subtractive) action is null (can be ignored).
    */
    template<class T> 
    constexpr T empty_sum(const T&) { 
      static_assert(std::is_arithmetic<T>::value, "Nullary products must be defined for this type");
      return T(); 
    }

    /*! \brief Returns the empty sum of std::complex types.
      
      The empty sum is a term whose additive (and typically its
      subtractive) action is null (can be ignored).
    */
    template<class T>
    constexpr std::complex<T> empty_sum(const std::complex<T>&) { 
      return std::complex<T>(T(), T()); 
    }

  }
}
