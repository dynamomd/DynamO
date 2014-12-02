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
    };

    /*! \brief Evaluates a sine Function expression at a given point.
    */
    template<class Arg, class Real>
    auto eval(const Function<Arg, detail::SIN>& f, const Real& x) -> decltype(std::sin(x))
    { return std::sin(x); }

    /*! \brief Evaluates a cosine Function expression at a given point.
    */
    template<class Arg, class Real>
    auto eval(const Function<Arg, detail::COS>& f, const Real& x) -> decltype(std::cos(x))
    { return std::cos(x); }

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
    auto operator-(const Function<A, Func>& f) -> decltype( -1 * f)
    { return -1 * f; }

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

    /*! \brief Generic left-handed subtraction operator for Function types.
     */
    template<class A, detail::Function_t Func, class RHS>
    auto operator-(const Function<A, Func>& f, const RHS& r) -> decltype(subtract(f, r)) 
    { return subtract(f, r); }

    /*! \brief Generic right-handed subtraction operator for Function types.
     */
    template<class Real, class A, detail::Function_t Func>
    auto operator-(const Real& r, const Function<A, Func>& f) -> decltype(subtract(r, f))
    { return subtract(r, f); }
    
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
    /*! \} */
    
    /*! \relates Function
      \name Function calculus
      \{
    */
    /*! \brief Generic derivative of sine functions.*/
    template<class A>
    auto derivative(const Function<A, detail::SIN>& f) -> decltype(derivative(f._arg) * cos(f._arg))
    { return derivative(f._arg) * cos(f._arg); }

    /*! \brief Derivative of a sine function with a constant argument.*/
    template<class T>
    Polynomial<0,T> derivative(const Function<Polynomial<0,T>, detail::SIN>& f)
    { return Polynomial<0, T>{}; }

    /*! \brief Generic derivative of cosine functions.*/
    template<class A>
    auto derivative(const Function<A, detail::COS>& f) -> decltype(-derivative(f._arg) * sin(f._arg))
    { return -derivative(f._arg) * sin(f._arg); }

    /*! \brief Derivative of a cosine function with a constant argument.*/
    template<class T>
    Polynomial<0,T> derivative(const Function<Polynomial<0,T>, detail::COS>& f)
    { return Polynomial<0, T>{}; }

    /*! \} */

    /*! \relates Function
      \name Function bounds
      \{
    */
    /*! \brief Estimates of the maximum and minimum value of a
        Function in a defined range.

      For sine and cosine, we return the trivial value [-1, +1]
      without attempting to specialise for sections of the
      oscillation.
     */
    template<class Arg, detail::Function_t Op, class Real>
    inline std::pair<double, double> minmax(const Function<Arg, Op>& f, const Real x_min, const Real x_max)
    {
      switch (Op) {
      case detail::SIN:
      case detail::COS: 
	return std::pair<double, double>(-1, 1);
      };
    }
    /* \} */
  }
}
