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
#include <magnet/math/symbolic.hpp>
#include <magnet/math/polynomial.hpp>
#include <magnet/math/operators.hpp>

namespace magnet {
  namespace math {
    namespace detail {
      /*! \relates Function
	\brief Function type template parameter.
       */
      typedef enum {
	SIN,
	COS
      } Function_t;
    }
    /*! \brief Function types are symbolic representation of
        non-polynomial functions.
    */
    template<class Arg, detail::Function_t Func>
    class Function {
    public:
      /*! \brief The symbolic expression which is the argument of the
          Function. 
      */
      Arg _arg;
      
      /*! \brief Construct a Function with its argument. */
      Function(const Arg& a): _arg(a) {}
      
      /*! \brief Evaluate a Function at a given value of \f$x\f$. */
      template<class Real>
      auto operator()(Real x) const ->decltype(_arg(x)) {
	switch (Func) {
	case detail::SIN: return std::sin(_arg(x));
	case detail::COS: return std::cos(_arg(x));
	}
      }
    };

    /*! \relates Function
      \name Function creation helper functions
      \{
    */
    /*! \brief Helper function for creating sine Function types. */
    template<class A>
    Function<A, detail::SIN> sin(const A& a) { return Function<A, detail::SIN>(a); }

    /*! \brief Helper function for creating cosine Function types. */
    template<class A>
    Function<A, detail::COS> cos(const A& a) { return Function<A, detail::COS>(a); }
    /*! \} */

    /*! \relates Function
      \name Function input/output operators
      \{
    */    
    /*! \brief Writes a human-readable representation of the Function
      to the output stream.
    */
    template<class A, detail::Function_t Func>
    inline std::ostream& operator<<(std::ostream& os, const Function<A, Func>& s) {
      switch (Func) {
      case detail::SIN: os << "sin("; break;
      case detail::COS: os << "cos("; break;
      }
      os << s._arg << ")";
      return os;
    }
    /*! \} */

    /*! \relates Function
      \name Function algebraic operators
     */    
    /*! \brief Unary negation operator for Function types.
     */
    template<class A, detail::Function_t Func>
    auto operator-(const Function<A, Func>& f) -> decltype(multiply(-1, f))
    { return multiply(-1, f); }

    /*! \brief Generic left-handed addition operator for Function types.
     */
    template<class A, detail::Function_t Func, class RHS>
    auto operator+(const Function<A, Func>& f, const RHS& r) -> decltype(add(f, r)) 
    { return add(f, r); }

    /*! \brief Generic right-handed addition operator for Function types.
     */
    template<class Real, class A, detail::Function_t Func>
    auto operator+(const Real& r, const Function<A, Func>& f) -> decltype(add(r, f))
    { return add(r, f); }
    
    /*! \brief Generic left-handed multiplication operator for Function types.
     */
    template<class A, detail::Function_t Func, class RHS>
    auto operator*(const Function<A, Func>& f, const RHS& r) -> decltype(multiply(f, r))
    { return multiply(f, r); }
    
    /*! \brief Generic right-handed multiplication operator for Function types.
     */
    template<class Real, class A, detail::Function_t Func>
    auto operator*(const Real& r, const Function<A, Func>& f) -> decltype(multiply(r, f))
    { return multiply(r, f); }

    /*! \brief Specific multiplication operator for two Function types
        to prevent ambiguous overloads.
     */
    template<class A, detail::Function_t Func,class A2, detail::Function_t Func2>
    auto operator*(const Function<A, Func>& f1, const Function<A2, Func2>& f2) -> decltype(multiply(f1, f2))
    { return multiply(f1, f2); }
    
    /*! \brief Specialisation for multiplication of Function and
      Polynomial types to prevent polynomials distributing Functions
      over all coefficients.
    */
    template<class A, detail::Function_t Func, class Real, size_t N>
    auto operator*(const Function<A, Func>& f1, const Polynomial<N, Real>& poly) -> decltype(multiply(f1, poly))
    { return multiply(f1, poly); }
    
    /*! \brief Specialisation for multiplication of Polynomial and
      Function types to prevent polynomials distributing Functions
      over all coefficients.

      This operator is to explicitly avoid distributing Function types
      over the coefficients of Polynomial types. Functions are usually
      expensive to evaluate, so we keep them factored out
    */
    template<class Real, size_t N, class A, detail::Function_t Func>
    auto operator*(const Polynomial<N, Real>& poly, const Function<A, Func>& f) -> decltype(multiply(poly, f))
    { return multiply(poly, f); }
    /*! \} */
    
    /*! \relates Function
      \name Function calculus
      \{
    */
    /*! \brief Generic derivative of sine functions.*/
    template<class A>
    auto derivative(const Function<A, detail::SIN>& f) -> decltype(derivative(f._arg) * cos(f._arg))
    { return derivative(f._arg) * cos(f._arg); }

    /*! \brief Generic derivative of cosine functions.*/
    template<class A>
    auto derivative(const Function<A, detail::COS>& f) -> decltype(-derivative(f._arg) * sin(f._arg))
    { return -derivative(f._arg) * sin(f._arg); }
    /*! \} */

    /*! \relates Function
      \name Function bounds
      \{
    */
    /*! \brief The maximum absolute value of a Function in a defined range.

      For sine and cosine, we return the trivial value 1 without
      attempting to specialise for sections of the oscillation.
     */
    template<class Arg, detail::Op_t Op, class Real>
    inline Real max_abs_val(const Function<Arg, Op>& f, const Real tmin, const Real tmax)
    {
      switch (Op) {
      case detail::SIN: return 1;
      case detail::COS: return 1;
      };
    }
    /* \} */
  }
}
