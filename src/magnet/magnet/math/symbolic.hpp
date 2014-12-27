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
#include <magnet/containers/stack_vector.hpp>
#include <complex>

namespace magnet {
  namespace math {
    /*!\brief Compile-time symbolic representation of zero.
     */
    struct NullSymbol {
      operator int () const { return 0; }

      /*! \brief Unary negation on NullSymbol has no action.*/
      NullSymbol operator-() const {
	return NullSymbol();
      }

      /*! \brief Unary positive on NullSymbol has no action.*/
      NullSymbol operator+() const {
	return NullSymbol();
      }
    };

    /*!\brief Compile-time symbolic representation of one.
     */
    struct UnitySymbol {
      operator int () const { return 1; }
    };

    /*!\brief Compile-time symbolic representation of a variable
      substitution.
    */
    template<char Letter, class Arg> struct VariableSubstitution {
      VariableSubstitution(const Arg& val):_val(val) {}
      Arg _val;
    };
    
    /*!\brief Symbolic representation of a variable.

      This class is used to denote a variable. The template argument
      is a single ASCII character which represents this variable and
      is used to identify it during symbolic actions and output.
    */
    template<char Letter> struct Variable {
      template<class Arg>
      VariableSubstitution<Letter, Arg> operator==(const Arg& a) const {
	return VariableSubstitution<Letter, Arg>(a);
      }
    };

    namespace detail {
      /*!\brief Type trait to determine if a certain type is a constant.

	This is used to enable the derivative operation to convert
	these types to NullSymbol types. It is also to apply a
	specialised functions to these types.
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

      The operation on constant terms is to do nothing.
    */
    template<class T, typename = typename std::enable_if<detail::IsConstant<T>::value> ::type>
    T expand(const T& f) { return f; }

    /*! \brief Provides expansion (and simplification) of symbolic
      functions.

      The operation on Variable terms is to do nothing.
    */
    template<char Letter>
    Variable<Letter> expand(const Variable<Letter>& f) { return f; }

    /*! \brief Evaluates a symbolic expression by substituting a
      variable for another expression.
	
      If a arithmetic type is substituted, this will likely cause a
      numerical evaluation of the expression. This "helper"
      implementation converts to a substitution for the variable
      'x'.
    */
    template<class T, class VarArg>
    auto eval(const T& f, const VarArg& xval) -> decltype(substitution(f, Variable<'x'>() == xval))
    { return substitution(f, Variable<'x'>() == xval); }

    /*! \brief Evaluates a symbolic expression using a substitution.
      
      This is just a synonym for substitution.
    */
    template<class T, char Letter, class Arg>
    auto eval(const T& f, const VariableSubstitution<Letter, Arg>& x)  -> decltype(substitution(f, x))
    { return substitution(f,x); }
    
    /*! \brief Default implementation of substitution of a symbolic
      expression at a given point.
      
      This implementation only applies if the term is a constant term.
    */
    template<class T, char Letter, class Arg,
	     typename = typename std::enable_if<detail::IsConstant<T>::value>::type >
    T substitution(const T& f, const VariableSubstitution<Letter, Arg>&)
    { return f; }

    /*! \brief Evaluates a symbolic Variable at a given point.
      
      This is only used if the Variable is the correct letter for the
      substitution.
    */
    template<char Letter, class Arg>
    Arg substitution(const Variable<Letter>& f, const VariableSubstitution<Letter, Arg>& x)
    { return x._val; }

    /*! \brief Evaluates a symbolic Variable at a given point.
      
      This is only used if the Variable is not the correct letter for the
      substitution.
    */
    template<char Letter1, class Arg, char Letter2,
	     typename = typename std::enable_if<Letter1 != Letter2>::type>
    Variable<Letter1> substitution(const Variable<Letter1>& f, const VariableSubstitution<Letter2, Arg>& x)
    { return f; }
    
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

    /*! \brief Output operator for Variable types. */
    template<char Letter>
    inline std::ostream& operator<<(std::ostream& os, const Variable<Letter>&) {
      os << Letter;
      return os;
    }

    /*! \brief Output operator for VariableSubstitution types. */
    template<char Letter, class Arg>
    inline std::ostream& operator<<(std::ostream& os, const VariableSubstitution<Letter, Arg>& sub) {
      os << Letter << " <- " << sub._val;
      return os;
    }
    
    /*! \} */
  
    /*! \brief Determine the derivative of a symbolic expression.
      
      This default implementation gives all consants
      derivative of zero.
    */
    template<class T, char Letter,
	     typename = typename std::enable_if<detail::IsConstant<T>::value>::type>
    NullSymbol derivative(const T&, Variable<Letter>) { return NullSymbol(); }

    /*! \brief Determine the derivative of a variable.

      If the variable is the variable in which a derivative is being
      taken, then this overload should be selected to return
      UnitySymbol.
    */
    template<char Letter1, char Letter2,
	     typename = typename std::enable_if<Letter1 == Letter2>::type>
    UnitySymbol derivative(Variable<Letter1>, Variable<Letter2>) { return UnitySymbol(); }

    inline containers::StackVector<double, 0> solve_real_roots(NullSymbol f) {
      return containers::StackVector<double, 0>();
    }

    /*! \brief Determine the derivative of a variable by another variable.

      If the variable is NOT the variable in which a derivative is
      being taken, then this overload should be selected to return
      NullSymbol.
    */
    template<char Letter1, char Letter2,
	     typename = typename std::enable_if<Letter1 != Letter2>::type>
    NullSymbol derivative(Variable<Letter1>, Variable<Letter2>) { return NullSymbol(); }
 
    /*! \brief Shift a function forward. It returns \f$g(x)=f(x+a)\f$

      For constant terms, these remain the same so this generic
      implementation does nothing.
    */
    template<class F, class Real,
	     typename = typename std::enable_if<detail::IsConstant<F>::value>::type>
    inline F shift_function(const F& f, const Real t) {
      return f;
    }
    
    /*! \brief Calculate the next real root of a symbolic function.

      For constant terms, +inf is returned to indicate no root was
      found.
    */
    template<class F, 
	     typename = typename std::enable_if<detail::IsConstant<F>::value>::type>
    inline double next_root(const F& f) {
      return HUGE_VAL;
    }

    /*! \brief Estimate the error in evaluating a function at a given time.
     */
    template<class F, class Real,
	     typename = typename std::enable_if<detail::IsConstant<F>::value>::type>
    inline double precision(const F& f, const Real) {
      return 0.0;
    }

    template<size_t Order, class Real = double, char Letter = 'x'> class Polynomial;

    /*! \brief Symbolic Factorial function.
     
      This template implementation returns UnitySymbol for 0! and 1!,
      allowing simplification of symbolic expressions.
     */
    template<size_t i> struct Factorial {
      static size_t eval() { return i * Factorial<i-1>::eval(); }
    };
    
    template<> struct Factorial<1> {
      static UnitySymbol eval() { return UnitySymbol(); }
    };

    template<> struct Factorial<0> {
      static UnitySymbol eval() { return UnitySymbol(); }
    };

    /*! \brief Symbolic Inverse factorial function.
     
      This template implementation returns UnitySymbol for 1/0! and 1/1!,
      allowing simplification of symbolic expressions.
     */
    template<size_t i> struct InvFactorial {
      static double eval() { return 1.0 / Factorial<i>::eval(); }
    };
    
    template<> struct InvFactorial<1> {
      static UnitySymbol eval() { return UnitySymbol(); }
    };

    template<> struct InvFactorial<0> {
      static UnitySymbol eval() { return UnitySymbol(); }
    };
    
    namespace detail {
      template<size_t State, size_t max_Order, char Letter>
      struct TaylorSeriesWorker {
	template<class Real>
	static NullSymbol eval(const NullSymbol& f, const Real& a)
	{ return NullSymbol(); }

	template<class F, class Real>
	static auto eval(const F& f, const Real& a) -> decltype(InvFactorial<State>::eval() * substitution(f, Variable<Letter>() == a) + (Variable<Letter>() - a) * TaylorSeriesWorker<State+1, max_Order, Letter>::eval(derivative(f, Variable<Letter>()), a))
	{
	  return InvFactorial<State>::eval() * substitution(f, Variable<Letter>() == a) + (Variable<Letter>()-a) * TaylorSeriesWorker<State+1, max_Order, Letter>::eval(derivative(f, Variable<Letter>()), a);
	}
      };

      template<size_t max_Order, char Letter>
      struct TaylorSeriesWorker<max_Order,max_Order,Letter> {
	template<class F, class Real>
	static auto eval(const F& f, const Real& a) -> decltype(InvFactorial<max_Order>::eval() * substitution(f, Variable<Letter>() == a))
	{ return InvFactorial<max_Order>::eval() * substitution(f, Variable<Letter>() == a); }

	template<class Real>
	static NullSymbol eval(const NullSymbol& f, const Real& a)
	{ return NullSymbol(); }
      };
    }

    /*! \brief Generate a Taylor series representation of a Symbolic
        expression.
     */
    template<size_t Order, char Letter, class F, class Real>
    auto taylor_series(const F& f, Real a) -> decltype(detail::TaylorSeriesWorker<0, Order, Letter>::eval(f, a))
    { return detail::TaylorSeriesWorker<0, Order, Letter>::eval(f, a); }
  }
}
