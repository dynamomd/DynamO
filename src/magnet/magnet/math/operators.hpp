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

namespace magnet {
  namespace math {
    namespace detail {
      /*! \brief Symbolic representation of a binary symbolic operation. 
       */
      template<class LHStype, class RHStype>
      struct BinaryOp {
	LHStype _l;
	RHStype _r;
	BinaryOp(const LHStype& l, const RHStype& r): _l(l), _r(r) {}
      };
    }

    /*! \brief Type trait which enables symbolic operators for
        algebraic operations (*+-). */
    template<class T> struct SymbolicOperators {
      static const bool value = false;
    };

    template<> struct SymbolicOperators<UnitySymbol> {
      static const bool value = true;
    };

     
    template<char Letter> struct SymbolicOperators<Variable<Letter> > {
      static const bool value = true;
    };

    template<> struct SymbolicOperators<NullSymbol> {
      static const bool value = true;
    };
     
    /*! \brief Type trait which denotes if operations should be
        reordered to bring these types together. 

	This is true for all arithmetic types, as operations on these
	generally can be collapsed into a single term.
    */
    template<class T1, class T2> struct Reorder {
      static const bool value = std::is_arithmetic<T1>::value && std::is_arithmetic<T2>::value;
    };
    
#define CREATE_BINARY_OP(HELPERNAME, CLASSNAME, OP)			\
    template<class LHStype, class RHStype>				\
    struct CLASSNAME : public detail::BinaryOp<LHStype, RHStype> {	\
      typedef detail::BinaryOp<LHStype, RHStype> Base;			\
      CLASSNAME(const LHStype& l, const RHStype& r): Base(l, r) {}	\
    };									\
									\
    template<class LHS, class RHS> struct SymbolicOperators<CLASSNAME<LHS,RHS> > { \
      static const bool value = true;					\
    };									\
									\
    template<class LHS, class RHS, class Real>				\
    auto eval(const CLASSNAME<LHS, RHS>& f, const Real& x)		\
      -> decltype(eval(f._l, x) OP eval(f._r, x))			\
    { return eval(f._l, x) OP eval(f._r, x); }				\
    template<class LHS, class RHS>					\
    /*! \brief Helper function for creating this BinaryOp. */		\
    CLASSNAME<LHS, RHS> HELPERNAME(const LHS& l, const RHS& r)		\
    { return CLASSNAME<LHS, RHS>(l, r); }				\
									\
    /*! \brief Pass the expand operator to the arguments of the operation */ \
    template<class LHS, class RHS>					\
    auto expand(const CLASSNAME<LHS, RHS>& f) -> decltype(expand(f._l) OP expand(f._r)) { \
      return expand(f._l) OP expand(f._r);				\
    }									\
									\
    /*! \brief Helper function which reorders (A*B)*C to (A*C)*B operations. */	\
    template<class T1, class T2, class T3,				\
	     typename = typename std::enable_if<Reorder<T1, T3>::value>::type>	\
    auto HELPERNAME(const CLASSNAME<T1, T2>& l, const T3& r)		\
      -> CLASSNAME<decltype(l._l OP r), T2>				\
    { return HELPERNAME(l._l OP r, l._r); }				\
    									\
    /*! \brief Helper function which reorders (A*B)*C to (B*C)*A operations. */	\
    template<class T1, class T2, class T3,				\
	     typename = typename std::enable_if<Reorder<T2, T3>::value>::type>	\
    auto HELPERNAME(const CLASSNAME<T1, T2>& l, const T3& r)		\
      -> CLASSNAME<decltype(l._r OP r), T1>				\
    { return HELPERNAME(l._r OP r, l._l); }				\
									\
    /*! \brief Helper function which reorders A*(B*C) to (A*B)*C operations. */	\
    template<class T1, class T2, class T3,				\
	     typename = typename std::enable_if<Reorder<T1, T2>::value>::type>	\
    auto HELPERNAME(const T1& l, const CLASSNAME<T2, T3>& r)		\
      -> CLASSNAME<decltype(l OP r._l), T3>				\
    { return HELPERNAME(l OP r._l, r._r); }				\
									\
    /*! \brief Helper function which reorders A*(B*C) to (A*C)*B operations. */	\
    template<class T1, class T2, class T3,				\
	     typename = typename std::enable_if<Reorder<T1, T3>::value>::type>	\
    auto HELPERNAME(const T1& l, const CLASSNAME<T2, T3>& r)		\
      -> CLASSNAME<decltype(l OP r._r), T2>				\
    { return HELPERNAME(l OP r._r, r._l); }				
    
    CREATE_BINARY_OP(add, AddOp, +)
    CREATE_BINARY_OP(subtract, SubtractOp, -)
    CREATE_BINARY_OP(multiply, MultiplyOp, *)
    CREATE_BINARY_OP(divide, DivideOp, /)

    /*! \relates BinaryOp
      \name BinaryOp optimisations
    */
    
    /*! \brief Optimised multiply if the LHS term is a NullSymbol. */
    template<class RHS>
    NullSymbol multiply(const NullSymbol& l, const RHS& r)
    { return NullSymbol(); }
    
    /*! \brief Optimised multiply if the RHS term is a NullSymbol. */
    template<class LHS>
    NullSymbol multiply(const LHS& l, const NullSymbol& r)
    { return NullSymbol(); }
    
    /*! \brief Optimised multiply if both terms are NullSymbol types. */
    NullSymbol multiply(const NullSymbol& l, const NullSymbol& r)
    { return NullSymbol(); }

    /*! \brief Optimised multiply if the LHS term is a UnitySymbol. */
    template<class RHS>
    RHS multiply(const UnitySymbol& l, const RHS& r)
    { return r; }
    
    /*! \brief Optimised multiply if the RHS term is a UnitySymbol. */
    template<class LHS>
    LHS multiply(const LHS& l, const UnitySymbol& r)
    { return l; }
    
    /*! \brief Optimised multiply if both terms are UnitySymbol types. */
    UnitySymbol multiply(const UnitySymbol&, const UnitySymbol&)
    { return UnitySymbol(); }

    /*! \brief Optimised addition if the LHS term is a NullSymbol. */
    template<class RHS>
    RHS add(const NullSymbol& l, const RHS& r)
    { return r; }
    
    /*! \brief Optimised addition if the RHS term is a NullSymbol. */
    template<class LHS>
    LHS add(const LHS& l, const NullSymbol& r)
    { return l; }
    
    /*! \brief Optimised addition if both terms are NullSymbol types. */
    NullSymbol add(const NullSymbol& l, const NullSymbol& r)
    { return NullSymbol(); }

    /*! \brief Optimised addition if the LHS term is a NullSymbol. */
    template<class RHS>
    auto subtract(const NullSymbol& l, const RHS& r) -> decltype(-r)
    { return -r; }
    
    /*! \brief Optimised addition if the RHS term is a NullSymbol. */
    template<class LHS>
    LHS subtract(const LHS& l, const NullSymbol& r)
    { return l; }
    
    /*! \brief Optimised addition if both terms are NullSymbol types. */
    NullSymbol subtract(const NullSymbol& l, const NullSymbol& r)
    { return NullSymbol(); }


    /*! \brief Expand multiplication and addition operations. */
    template<class LHS1, class RHS1, class RHS>
    auto expand(const MultiplyOp<AddOp<LHS1, RHS1>, RHS>& f)
      -> decltype(expand(f._l._l * f._r + f._l._r * f._r)) {
      return expand(f._l._l * f._r + f._l._r * f._r);
    }

    /*! \brief Expand multiplication and addition operations. */
    template<class LHS1, class RHS1, class RHS>
    auto expand(const MultiplyOp<RHS, AddOp<LHS1, RHS1> >& f)
      -> decltype(expand(f._l * f._r._l + f._l * f._r._r)) {
      return expand(f._l * f._r._l + f._l * f._r._r);
    }

    /*! \brief Expand multiplication and addition operations. */
    template<class LHS1, class RHS1, class LHS2, class RHS2>
    auto expand(const MultiplyOp<AddOp<LHS1, RHS1>,  AddOp<LHS2, RHS2> >& f)
      -> decltype(expand(f._l._l * f._r._l + f._l._l * f._r._r + f._l._r * f._r._l + f._l._r * f._r._r)) {
      return expand(f._l._l * f._r._l + f._l._l * f._r._r + f._l._r * f._r._l + f._l._r * f._r._r);
    }

    /*! \} */

    /*! \name Symbolic algebra
      \{
    */

    /*! \brief Simple combination rule to enable Symbolic operations,
        but avoid redundantly specifying where two SymbolicOperator
        classes are operated on. */
    template<class LHS, class RHS> 
    struct ApplySymbolicOps {
      static const bool value = SymbolicOperators<LHS>::value || (!SymbolicOperators<LHS>::value && SymbolicOperators<RHS>::value);
    };

    /*! \brief Symbolic addition operator. */
    template<class LHS, class RHS,
	     typename = typename std::enable_if<ApplySymbolicOps<LHS, RHS>::value>::type>
    auto operator+(const LHS& l, const RHS& r) -> decltype(add(l, r))
    { return add(l,r); }

    /*! \brief Symbolic multiplication operator. */
    template<class LHS, class RHS,
	     typename = typename std::enable_if<ApplySymbolicOps<LHS, RHS>::value>::type>
    auto operator*(const LHS& l, const RHS& r) -> decltype(multiply(l, r))
    { return multiply(l,r); }

    /*! \brief Symbolic subtraction operator. */
    template<class LHS, class RHS,
	     typename = typename std::enable_if<ApplySymbolicOps<LHS, RHS>::value>::type>
    auto operator-(const LHS& l, const RHS& r) -> decltype(subtract(l, r))
    { return subtract(l,r); }

    /*! \brief Symbolic divide operator. */
    template<class LHS, class RHS,
	     typename = typename std::enable_if<ApplySymbolicOps<LHS, RHS>::value>::type>
    auto operator/(const LHS& l, const RHS& r) -> decltype(divide(l, r))
    { return divide(l,r); }

    /*! \brief Derivatives of AddOp operations.
     */
    template<char dVariable, class LHS, class RHS>
    auto derivative(const AddOp<LHS, RHS>& f, Variable<dVariable>) -> decltype(derivative(f._l, Variable<dVariable>()) + derivative(f._r, Variable<dVariable>()))
    { return derivative(f._l, Variable<dVariable>()) + derivative(f._r, Variable<dVariable>()); }

    /*! \brief Derivatives of SubtractOp operations.
     */
    template<char dVariable, class LHS, class RHS>
    auto derivative(const SubtractOp<LHS, RHS>& f, Variable<dVariable>) -> decltype(derivative(f._l, Variable<dVariable>()) - derivative(f._r, Variable<dVariable>()))
    { return derivative(f._l, Variable<dVariable>()) - derivative(f._r, Variable<dVariable>()); }

    /*! \brief Derivatives of MultiplyOp operations.
     */
    template<char dVariable, class LHS, class RHS>
    auto derivative(const MultiplyOp<LHS, RHS>& f, Variable<dVariable>) -> decltype(derivative(f._l, Variable<dVariable>()) * f._r + f._l * derivative(f._r, Variable<dVariable>()))
    { return derivative(f._l, Variable<dVariable>()) * f._r + f._l * derivative(f._r, Variable<dVariable>()); }

    /*! \} */

    /*! \brief Determines the max and min over a certain range. */
    template<class LHS, class RHS, class Real>
    auto minmax(const AddOp<LHS, RHS>& f, const Real x_min, const Real x_max) -> std::pair<decltype(minmax(f._l, x_min, x_max).first + minmax(f._r, x_min, x_max).first), 
											   decltype(minmax(f._l, x_min, x_max).second + minmax(f._r, x_min, x_max).second)>
    {
      typedef std::pair<decltype(minmax(f._l, x_min, x_max).first + minmax(f._r, x_min, x_max).first), 
			decltype(minmax(f._l, x_min, x_max).second + minmax(f._r, x_min, x_max).second)> 
	RetType;
      auto lv = minmax(f._l, x_min, x_max);
      auto rv = minmax(f._r, x_min, x_max);
      return RetType(lv.first + rv.first, lv.second + rv.second);
    }

    /*! \brief Determines the max and min over a certain range. */
    template<class LHS, class RHS, class Real>
    auto minmax(const MultiplyOp<LHS, RHS>& f, const Real x_min, const Real x_max) -> std::pair<decltype(minmax(f._l, x_min, x_max).first * minmax(f._r, x_min, x_max).first), 
												decltype(minmax(f._l, x_min, x_max).second * minmax(f._r, x_min, x_max).second)>
    {
      typedef std::pair<decltype(minmax(f._l, x_min, x_max).first * minmax(f._r, x_min, x_max).first), 
			decltype(minmax(f._l, x_min, x_max).second * minmax(f._r, x_min, x_max).second)> RetType;
      auto lv = minmax(f._l, x_min, x_max);
      auto rv = minmax(f._r, x_min, x_max);
      return RetType(lv.first * rv.first, lv.second * rv.second);
    }
    /*! \} */

    /*! \relates BinaryOp
      \name BinaryOp input/output operators
      \{
    */
    /*! \brief Writes a human-readable representation of the AddOp to the output stream. */
    template<class LHS, class RHS>
    inline std::ostream& operator<<(std::ostream& os, const AddOp<LHS, RHS>& op) {
      os << "(" << op._l << ") + (" << op._r << ")";
      return os;
    }

    /*! \brief Writes a human-readable representation of the SubtractOp to the output stream. */
    template<class LHS, class RHS>
    inline std::ostream& operator<<(std::ostream& os, const SubtractOp<LHS, RHS>& op) {
      os << "(" << op._l << ") - (" << op._r << ")";
      return os;
    }

    /*! \brief Writes a human-readable representation of the MultiplyOp to the output stream. */
    template<class LHS, class RHS>
    inline std::ostream& operator<<(std::ostream& os, const MultiplyOp<LHS, RHS>& op) {
      os << "(" << op._l << ") * (" << op._r << ")";
      return os;
    }

    /*! \brief Writes a human-readable representation of the DivideOp to the output stream. */
    template<class LHS, class RHS>
    inline std::ostream& operator<<(std::ostream& os, const DivideOp<LHS, RHS>& op) {
      os << "(" << op._l << ") / (" << op._r << ")";
      return os;
    }
    /*! \} */

    namespace {
      /*! \brief Generic implementation of the eval routine for PowerOp.
	
	As the types of non-arithmetic arguments to PowerOp might
	change with each round of multiplication, we must be careful
	to accommodate this using templated looping. This class
	achieves this.
      */
      template<size_t Power>
      struct PowerOpEval {
	template<class Arg_t>
	static auto eval(Arg_t x) -> decltype(PowerOpEval<Power-1>::eval(x) * x) {
	  return PowerOpEval<Power-1>::eval(x) * x;
	}
	static UnitySymbol eval(UnitySymbol) { return UnitySymbol(); }
      };

      template<>
      struct PowerOpEval<1> {
	template<class Arg_t>
	static Arg_t eval(Arg_t x) {
	  return x;
	}
	static UnitySymbol eval(UnitySymbol) { return UnitySymbol(); }
      };

      template<>
      struct PowerOpEval<0> {
	template<class Arg_t>
	static UnitySymbol eval(Arg_t x) {
	  return UnitySymbol();
	}
      };
    }

    /*! \brief Symbolic representation of a (positive) power operator.
     */
    template<class Arg, size_t Power>
    struct PowerOp {
      Arg _arg;
      PowerOp(Arg a): _arg(a) {}
    };
    
    /*! \brief Evaluate PowerOp symbol at a value of x.
      
      This operator is only used if the result of evaluating the
      argument is an arithmetic type. If this is the case, the
      evaluation is passed to std::pow.
    */
    template<class Arg, size_t Power, class Real>
    auto eval(const PowerOp<Arg, Power>& f, const Real& x) -> typename std::enable_if<std::is_arithmetic<decltype(eval(f._arg, x))>::value, decltype(std::pow(eval(f._arg, x), Power))>::type
    { return std::pow(eval(f._arg, x), Power); }

    /*! \brief Evaluate PowerOp symbol at a value of x.
      
      This is the general implementation for PowerOp.
    */
    template<class Arg, size_t Power, class Real>
    auto eval(const PowerOp<Arg, Power>& f, const Real& x) -> typename std::enable_if<!std::is_arithmetic<decltype(f._arg(x))>::value, decltype(PowerOpEval<Power>::eval(eval(f._arg, x)))>
    { return PowerOpEval<Power>::eval(eval(f._arg, x)); }


    /*! \relates PowerOp
      \name PowerOp helper functions.
    */
    /*! \brief Helper function for creating PowerOp types. */
    template<size_t N, class Arg>
    PowerOp<Arg, N> pow(const Arg& f)
    { return PowerOp<Arg, N>(f); }

    /*! \} */

    /*! \relates PowerOp
      \name PowerOp input/output operators.
    */
    /*! \brief Writes a human-readable representation of the BinaryOp to the output stream. */
    template<class Arg, size_t Power>
    inline std::ostream& operator<<(std::ostream& os, const PowerOp<Arg, Power>& p) {
      os << "(" << p._arg << ")^" << Power;
      return os;
    }
    /*! \} */

    
    /*! \relates PowerOp
      \name PowerOp operations
      \{
      \brief Enablement of the symbolic operators for the PowerOp type.
    */
    template<class Arg, size_t Power>
    struct SymbolicOperators<PowerOp<Arg, Power> > {
      static const bool value = true;
    };

    /*! \brief Expansion operator for PowerOp types. */
    template<class Arg, size_t Power> 
    auto expand(const PowerOp<Arg, Power>& f) -> decltype(PowerOpEval<Power>::eval(f._arg))
    { return PowerOpEval<Power>::eval(f._arg); }


    /*! \brief Derivatives of PowerOp operations.
     */
    template<char dVariable, class Arg, size_t Power>
    auto derivative(const PowerOp<Arg, Power>& f, Variable<dVariable>) -> decltype(derivative(f._arg, Variable<dVariable>()) * PowerOp<Arg, Power-1>(f._arg))
    { return Power * derivative(f._arg, Variable<dVariable>()) * PowerOp<Arg, Power-1>(f._arg); }

    template<char dVariable, class Arg>
    auto derivative(const PowerOp<Arg, 1>& f, Variable<dVariable>) -> decltype(derivative(f._arg, Variable<dVariable>()))
    { return derivative(f._arg, Variable<dVariable>()); }

    template<char dVariable, class Arg>
    auto derivative(const PowerOp<Arg, 2>& f, Variable<dVariable>) -> decltype(derivative(f._arg, Variable<dVariable>()) * f._arg)
    { return 2 * derivative(f._arg, Variable<dVariable>()) * f._arg; }

    template<char dVariable, class Arg>
    NullSymbol derivative(const PowerOp<Arg, 0>& f, Variable<dVariable>)
    { return NullSymbol(); }

    /*! \brief The maximum and minimum values of the PowerOp over a specifed range.
      \param f The PowerOp operation to evaluate.
      \param x_min The minimum x bound.
      \param x_max The maximum x bound.
    */
    template<class Arg, size_t Power, class Real>
    inline auto minmax(const PowerOp<Arg, Power>& f, const Real x_min, const Real x_max) -> std::pair<decltype(PowerOpEval<Power>::eval(minmax(f._arg, x_min, x_max).first)),
												      decltype(PowerOpEval<Power>::eval(minmax(f._arg, x_min, x_max).second))>
    {
      typedef std::pair<decltype(PowerOpEval<Power>::eval(minmax(f._arg, x_min, x_max).first)), decltype(PowerOpEval<Power>::eval(minmax(f._arg, x_min, x_max).second))> RetType;
      auto val = minmax(f._arg, x_min, x_max);
      
      auto min_pow = PowerOpEval<Power>::eval(val.first);
      auto max_pow = PowerOpEval<Power>::eval(val.second);
      
      //For odd powers, sign is preserved, so the arguments min^Power
      //is always less than the arguments max^Power
      if (Power % 2)
	return RetType(min_pow, max_pow);
      else {
	auto min = std::min(min_pow, max_pow);
	auto max = std::max(min_pow, max_pow);
	//If min-max range spans zero, we must include it.
	if ((val.first < 0) && (val.second > 0))
	  min = std::min(min, 0.0);

	return RetType(min, max);
      }
    }
    /*! \} */
  }
}
    
