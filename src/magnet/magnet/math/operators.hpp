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

    /*! \brief Type trait which denotes BinaryOp types. */
    template<class T> struct IsOp {
      static const bool value = false;
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
    template<class LHS, class RHS> struct IsOp<CLASSNAME<LHS,RHS> > {	\
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

    /*! \} */

    /*! \relates BinaryOp
      \name BinaryOp algebra
      \{
    */

    /*! \brief Left-handed addition operator for BinaryOp types. */
    template<class Op, class RHS>
    auto operator+(const Op& l, const RHS& r) -> typename std::enable_if<IsOp<Op>::value, decltype(add(l, r))>::type
    { return add(l,r); }

    /*! \brief Right-handed addition operator for BinaryOp types. */
    template<class Op, class LHS>
    auto operator+(const LHS& l, const Op& r) -> typename std::enable_if<IsOp<Op>::value && ! IsOp<LHS>::value, decltype(add(l, r))>::type
    { return add(l,r); }

    /*! \brief Left-handed multiplication operator for BinaryOp types. */
    template<class Op, class RHS>
    auto operator*(const Op& l, const RHS& r) -> typename std::enable_if<IsOp<Op>::value, decltype(multiply(l, r))>::type
    { return multiply(l,r); }

    /*! \brief Right-handed multiplication operator for BinaryOp types. */
    template<class Op, class LHS>
    auto operator*(const LHS& l, const Op& r) -> typename std::enable_if<IsOp<Op>::value && ! IsOp<LHS>::value, decltype(multiply(l, r))>::type
    { return multiply(l,r); }

    /*! \brief Left-handed subtraction operator for BinaryOp types. */
    template<class Op, class RHS>
    auto operator-(const Op& l, const RHS& r) -> typename std::enable_if<IsOp<Op>::value, decltype(subtract(l, r))>::type
    { return subtract(l,r); }

    /*! \brief Right-handed subtraction operator for BinaryOp types. */
    template<class Op, class LHS>
    auto operator-(const LHS& l, const Op& r) -> typename std::enable_if<IsOp<Op>::value && ! IsOp<LHS>::value, decltype(subtract(l, r))>::type
    { return subtract(l,r); }

    /*! \brief Left-handed subtraction operator for BinaryOp types. */
    template<class Op, class RHS>
    auto operator/(const Op& l, const RHS& r) -> typename std::enable_if<IsOp<Op>::value, decltype(divide(l, r))>::type
    { return divide(l,r); }

    /*! \brief Right-handed subtraction operator for BinaryOp types. */
    template<class Op, class LHS>
    auto operator/(const LHS& l, const Op& r) -> typename std::enable_if<IsOp<Op>::value && ! IsOp<LHS>::value, decltype(divide(l, r))>::type
    { return divide(l,r); }

    /*! \brief Expand addition BinaryOp types.

      If the classes have specialised operators for addition, then the
      decltype lookup will succeed and the addition is shunted to
      those classes. If not, this lookup will fail to expand the
      addition and it is instead carried out by the BinaryOp class.
    */
    template<class LHS, class RHS>
    auto expand(const AddOp<LHS, RHS>& f) -> decltype(expand(f._l) + expand(f._r)) {
      return expand(f._l) + expand(f._r);
    }

    /*! \brief Expand multiplication BinaryOp types.

      If the classes have specialised operators for multiplication, then the
      decltype lookup will succeed and the addition is shunted to
      those classes. If not, this lookup will fail to expand the
      addition and it is instead carried out by the BinaryOp class.
    */
    template<class LHS, class RHS>
    auto expand(const MultiplyOp<LHS, RHS>& f) -> decltype(expand(f._l) * expand(f._r)) {
      return expand(f._l) * expand(f._r);
    }

    /*! \brief Derivatives of AddOp operations.
     */
    template<class LHS, class RHS>
    auto derivative(const AddOp<LHS, RHS>& f) -> decltype(derivative(f._l) + derivative(f._r))
    { return derivative(f._l) + derivative(f._r); }

    /*! \brief Derivatives of SubtractOp operations.
     */
    template<class LHS, class RHS>
    auto derivative(const SubtractOp<LHS, RHS>& f) -> decltype(derivative(f._l) - derivative(f._r))
    { return derivative(f._l) - derivative(f._r); }

    /*! \brief Derivatives of MultiplyOp operations.
     */
    template<class LHS, class RHS>
    auto derivative(const MultiplyOp<LHS, RHS>& f) -> decltype(derivative(f._l) * f._r + derivative(f._r) * f._l)
    { return derivative(f._l) * f._r + derivative(f._r) * f._l; }


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
      os << op._l << " + " << op._r ;
      return os;
    }

    /*! \brief Writes a human-readable representation of the SubtractOp to the output stream. */
    template<class LHS, class RHS>
    inline std::ostream& operator<<(std::ostream& os, const SubtractOp<LHS, RHS>& op) {
      os << op._l << " - " << op._r ;
      return os;
    }

    /*! \brief Writes a human-readable representation of the MultiplyOp to the output stream. */
    template<class LHS, class RHS>
    inline std::ostream& operator<<(std::ostream& os, const MultiplyOp<LHS, RHS>& op) {
      os << "{" << op._l << " * " << op._r << "}";
      return os;
    }

    /*! \brief Writes a human-readable representation of the DivideOp to the output stream. */
    template<class LHS, class RHS>
    inline std::ostream& operator<<(std::ostream& os, const DivideOp<LHS, RHS>& op) {
      os << op._l << " / {" << op._r << "}" ;
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
      };

      template<>
      struct PowerOpEval<1> {
	template<class Arg_t>
	static Arg_t eval(Arg_t x) {
	  return x;
	}
      };

      template<>
      struct PowerOpEval<0> {
	template<class Arg_t>
	static double eval(Arg_t x) {
	  return 1;
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
    
//    /*! \brief Evaluate PowerOp symbol at a value of x.
//      
//      This operator is only used if the result of evaluating the
//      argument is an arithmetic type. If this is the case, the
//      evaluation is passed to std::pow.
//    */
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
      \name PowerOp algebra operations
    */

    /*! \brief Left-handed multiplication operator for PowerOp types. */
    template<class Arg, size_t Power> 
    auto expand(const PowerOp<Arg, Power>& f) -> decltype(PowerOpEval<Power>::eval(f._arg))
    { return PowerOpEval<Power>::eval(f._arg); }

    /*! \brief Left-handed multiplication operator for PowerOp types. */
    template<class Arg, size_t Power, class RHS>
    auto operator*(const PowerOp<Arg, Power>& l, const RHS& r) -> decltype(multiply(l, r))
    { return multiply(l, r); }
    
    /*! \brief Right-handed multiplication operator for PowerOp types. */
    template<class Arg, size_t Power, class RHS>
    auto operator*(const RHS& l, const PowerOp<Arg, Power>& r) -> decltype(multiply(l, r))
    { return multiply(l, r); }
    
    /*! \brief Multiplication operator for two PowerOp types. */
    template<class Arg1, size_t Power1, class Arg2, size_t Power2>
    auto operator*(const PowerOp<Arg1, Power1>& l, const PowerOp<Arg2, Power2>& r) -> decltype(multiply(l, r))
    { return multiply(l, r); }

    /*! \brief Left-handed addition operator for PowerOp types. */
    template<class Arg, size_t Power, class RHS>
    auto operator+(const PowerOp<Arg, Power>& l, const RHS& r) -> decltype(add(l, r))
    { return add(l, r); }
    
    /*! \brief Right-handed addition operator for PowerOp types. */
    template<class Arg, size_t Power, class RHS>
    auto operator+(const RHS& l, const PowerOp<Arg, Power>& r) -> decltype(add(l, r))
    { return add(l, r); }
    
    /*! \brief Addition operator for two PowerOp types. */
    template<class Arg1, size_t Power1, class Arg2, size_t Power2>
    auto operator+(const PowerOp<Arg1, Power1>& l, const PowerOp<Arg2, Power2>& r) -> decltype(add(l, r))
    { return add(l, r); }

    /*! \brief Left-handed subtraction operator for PowerOp types. */
    template<class Arg, size_t Power, class RHS>
    auto operator-(const PowerOp<Arg, Power>& l, const RHS& r) -> decltype(subtract(l, r))
    { return subtract(l, r); }
    
    /*! \brief Right-handed subtraction operator for PowerOp types. */
    template<class Arg, size_t Power, class RHS>
    auto operator-(const RHS& l, const PowerOp<Arg, Power>& r) -> decltype(subtract(l, r))
    { return subtract(l, r); }
    
    /*! \brief Subtraction operator for two PowerOp types. */
    template<class Arg1, size_t Power1, class Arg2, size_t Power2>
    auto operator-(const PowerOp<Arg1, Power1>& l, const PowerOp<Arg2, Power2>& r) -> decltype(subtract(l, r))
    { return subtract(l, r); }

    /*! \brief Derivatives of PowerOp operations.
     */
    template<class Arg, size_t Power>
    auto derivative(const PowerOp<Arg, Power>& f) -> decltype(derivative(f._arg) * PowerOp<Arg, Power-1>(f._arg))
    { return Power * derivative(f._arg) * PowerOp<Arg, Power-1>(f._arg); }

    template<class Arg>
    auto derivative(const PowerOp<Arg, 1>& f) -> decltype(derivative(f._arg))
    { return derivative(f._arg); }

    template<class Arg>
    auto derivative(const PowerOp<Arg, 2>& f) -> decltype(derivative(f._arg) * f._arg)
    { return 2 * derivative(f._arg) * f._arg; }

    template<class Arg>
    double derivative(const PowerOp<Arg, 0>& f)
    { return 0; }

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
    
