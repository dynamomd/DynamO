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
#include <magnet/math/vector.hpp>
#include <complex>

namespace magnet {
  namespace math {
    /*!\brief Compile-time symbolic representation of zero.
     */
    struct NullSymbol {
      operator int () const { return 0; }
    };

    /*!\brief Compile-time symbolic representation of one.
     */
    struct UnitySymbol {
      operator int () const { return 1; }
    };

    /*!\brief Symbolic representation of a variable.
     */
    struct Var {
    };

    namespace detail {
      /*!\brief Type trait to determine if a certain type is a constant.

	This is used to enable the derivative operation to convert
	these types to NullSymbol types. It is also to apply a
	specialised eval and minmax function to these types.
       */
      template<class T>
      struct IsConstant {
	static const bool value = std::is_arithmetic<T>::value;
      };

      template<class T, size_t N>
      struct IsConstant<NVector<T, N> > {
	static const bool value = true;
      };

      template<class T>
      struct IsConstant<std::complex<T> > {
	static const bool value = true;
      };

      template<>
      struct IsConstant<NullSymbol> {
	static const bool value = true;
      };

      template<>
      struct IsConstant<UnitySymbol> {
	static const bool value = true;
      };
    }

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

    /*! \brief Returns the empty product of a NVector type.
      
      The empty product is a term whose multiplicative action is null
      (can be ignored).
    */
    template<class T, size_t N>
    constexpr NVector<T, N> empty_product(const NVector<T, N>&) { 
      return NVector<T, N>(T(1));
    }

    /*! \brief Returns the empty sum of a type.
      
      The empty sum is a term whose additive (and typically its
      subtractive) action is null (can be ignored). This definition
      only applies for selected types and assumes that their default
      constructors create the empty sum.
    */
    template<class T>
    constexpr auto empty_sum(const T&) -> typename std::enable_if<detail::IsConstant<T>::value, T>::type { 
      return T();
    }

    /*! \brief Provides expansion (and simplification) of symbolic
      functions.

      The purpose of this function is to reduce the complexity of
      symbolic expressions to accelerate any successive
      evaluations. This should not change the calculated values, but
      should optimise for use under repeated evaluations.

      The default operation is to do nothing.
    */
    template<class T> 
    const T& expand(const T& f) { return f; }

    /*! \brief Evaluates a symbolic expression at a given point.

      This generic implementation is used for constant terms.
    */
    template<class T, class Real> 
    auto eval(const T& f, const Real& x) -> typename std::enable_if<detail::IsConstant<T>::value, const T&>::type
    { return f; }

    /*! \brief Evaluates a symbolic Var variable expression at a given point.
    */
    template<class Real> Real eval(const Var& f, const Real& x) { return x; }
    
    /*! \brief Output operator for NullSymbol types. */
    inline std::ostream& operator<<(std::ostream& os, const NullSymbol&) {
      os << "Null";
      return os;
    }

    /*! \brief Output operator for UnitySymbol types. */
    inline std::ostream& operator<<(std::ostream& os, const UnitySymbol&) {
      os << "Unity";
      return os;
    }

    /*! \brief Output operator for Var types. */
    inline std::ostream& operator<<(std::ostream& os, const Var&) {
      os << "v";
      return os;
    }
    
    /*! \} */
  
    /*! \brief Determine the derivative of a symbolic expression.
      
      This default implementation only applies for arithmetic
      types. As most symbols are constants or arithmetic types (int,
      float, double) by default their derivatives are zero.
    */
    template<class T, typename = typename std::enable_if<detail::IsConstant<T>::value> >
    NullSymbol derivative(const T&) { return NullSymbol();}

    /*! \brief Determine the derivative of a symbolic expression.
     */
    UnitySymbol derivative(const Var&) { return UnitySymbol(); }
    
    /*! \brief Determine the minimum and maximum of a symbolic
      expression within some bounds of x.
      
      For arithmetic types (int, float, double) their bounds are
      their actual values.
    */
    template<class T, class Real, typename = typename std::enable_if<detail::IsConstant<T>::value>::type>
    std::pair<T, T> minmax(const T& f, const Real& x_min, const Real& x_max)
    { return std::pair<T, T>(f, f); }

    /*! \brief Determine the minimum and maximum of a Var expression
      within some bounds of x.
    */
    template<class Real>
    std::pair<Real, Real> minmax(const Var& f, const Real& x_min, const Real& x_max)
    { return std::pair<Real, Real>(x_min, x_max); }
  }
}
