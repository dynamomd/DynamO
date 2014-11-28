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
#include <magnet/math/operators.hpp>

namespace magnet {
  namespace math {
    namespace detail {
      typedef enum {
	SIN,
	COS
      } Function_t;
    }
    /*! \brief Symbolic representation of non-polynomial functions.
    */
    template<class Arg, detail::Function_t Func>
    class Function {
    public:
      Arg _arg;

      Function(const Arg& a): _arg(a) {}
      
      template<class Real>
      auto operator()(Real x) const ->decltype(_arg(x)) {
	switch (Func) {
	case detail::SIN: return std::sin(_arg(x));
	case detail::COS: return std::cos(_arg(x));
	}
      }

      auto operator-() const -> decltype(multiply(-1, *this)) {
	return multiply(-1, *this);
      }

      template<class RHS>
      auto operator+(const RHS& r) const -> decltype(add(*this, r)) {
	return add(*this, r);
      }

      template<class RHS>
      auto operator*(const RHS& r) const -> decltype(multiply(*this, r)) {
	return multiply(*this, r);
      }
    };

    /*! \relates Function
      \brief Helper function for creating sine Function types
    */
    template<class A>
    Function<A, detail::SIN> sin(const A& a) { return Function<A, detail::SIN>(a); }

    /*! \relates Function
      \brief Helper function for creating cosine Function types.
    */
    template<class A>
    Function<A, detail::COS> cos(const A& a) { return Function<A, detail::COS>(a); }
    
    /*! \relates Function 
      
      \brief Writes a human-readable representation of the Function to
      the output stream.
    */
    template<class A, detail::Function_t Func>
    inline std::ostream& operator<<(std::ostream& os, const Function<A, Func>& s) {
      switch (Func) {
      case detail::SIN: os << " sin("; break;
      case detail::COS: os << " cos("; break;
      }
      os << s._arg << ")";
      return os;
    }

    /*! \relates Function
      
      \brief Here we explicitly avoid distributing Function types over
      the coefficients of Polynomial types.

      Functions are usually expensive to evaluate, so we keep them factored out
    */
    template<class Real, size_t N, class A, detail::Function_t Func>
    auto operator*(const Polynomial<N, Real>& poly, const Function<A, Func>& f) -> decltype(multiply(poly, f))
    { return multiply(poly, f); }

    /*! \relates Function
      \name Function calculus
      \{
    */
    /*! \brief Derivative of sine functions.*/
    template<class A>
    auto derivative(const Function<A, detail::SIN>& f) -> decltype(derivative(f._arg) * cos(f._arg))
    { return derivative(f._arg) * cos(f._arg); }

    /*! \brief Derivative of cosine functions.*/
    template<class A>
    auto derivative(const Function<A, detail::COS>& f) -> decltype(-derivative(f._arg) * sin(f._arg))
    { return -derivative(f._arg) * sin(f._arg); }
    /*! \} */
  }
}
