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

namespace magnet {
  namespace math {
    namespace detail {
      typedef enum {
	ADD,
	MULTIPLY,
	DIVIDE,
      } Op_t;
    }
    /*! \brief Symbolic representation of a binary operator. 
      
      When dealing with multiple symbols (Polynomial or Sin terms), it
      is convenient to have a representation of operators between
      them. This class represents these operations.
     */
    template<class LHStype, class RHStype, detail::Op_t Op>
    struct BinaryOp {
      LHStype _l;
      RHStype _r;
      
      BinaryOp(LHStype l, RHStype r): _l(l), _r(r) {}
      
      template<class R>
      auto operator()(const R& x) const -> decltype(_l(x) + _r(x)) {
	switch (Op){
	case detail::ADD: return _l(x) + _r(x);
	case detail::MULTIPLY: return _l(x) * _r(x);
	case detail::DIVIDE: return _l(x) / _r(x);
	}
      }
    };

    /*! \brief Provides expansion (and simplification) of symbolic
        functions.

      The purpose of this function is to reduce the complexity of
      symbolic expressions to accelerate any successive
      evaluations. This should not change the calculated values, but
      should optimise for use under repeated evaluations.

      The default operation is to do nothing.
     */
    template<class T> const T& expand(const T& f) { return f; }

    /*! \relates BinaryOp
      \brief Helper function for creation of addition BinaryOp types.
    */
    template<class LHS, class RHS>
    BinaryOp<LHS, RHS, detail::ADD> add(const LHS& l, const RHS& r)
    { return BinaryOp<LHS, RHS, detail::ADD>(l, r); }

    /*! \relates BinaryOp
      \brief Helper function for creation of subtraction BinaryOp types.
    */
    template<class LHS, class RHS>
    BinaryOp<LHS, RHS, detail::ADD> subtract(const LHS& l, const RHS& r)
    { return BinaryOp<LHS, RHS, detail::ADD>(l, -r); }

    /*! \relates BinaryOp
      \brief Helper function for creation of multiply BinaryOp types.
    */
    template<class LHS, class RHS>
    BinaryOp<LHS, RHS, detail::MULTIPLY> multiply(const LHS& l, const RHS& r)
    { return BinaryOp<LHS, RHS, detail::MULTIPLY>(l, r); }

    /*! \relates BinaryOp
      \brief Helper function for creation of multiply BinaryOp types.
    */
    template<class LHS, class RHS>
    BinaryOp<LHS, RHS, detail::ADD> divide(const LHS& l, const RHS& r)
    { return BinaryOp<LHS, RHS, detail::DIVIDE>(l, r); }

    /*! \relates BinaryOp
      \name BinaryOp algebra
      \{
    */

    /*! \brief Left-handed multiplication operator for BinaryOp types. */
    template<class LHS, class RHS, detail::Op_t Op, class RRHS>
    auto operator*(const BinaryOp<LHS, RHS, Op>& l, const RRHS& r) -> decltype(multiply(l,r))
    { return multiply(l,r); }

    /*! \brief Right-handed multiplication operator for BinaryOp types. */
    template<class LHS, class RHS, detail::Op_t Op, class LLHS>
    auto operator*(const LLHS& l, const BinaryOp<LHS, RHS, Op>& r) -> decltype(multiply(l,r))
    { return multiply(l,r); }

    /*! \brief Multiplication operator for two BinaryOp types. */
    template<class LHS1, class RHS1, detail::Op_t Op1, class LHS2, class RHS2, detail::Op_t Op2>
    auto operator*(const BinaryOp<LHS1, RHS1, Op1>& l, const BinaryOp<LHS2, RHS2, Op2>& r) -> decltype(multiply(l,r))
    { return multiply(l,r); }

    /*! \brief Left-handed addition operator for BinaryOp types. */
    template<class LHS, class RHS, detail::Op_t Op, class RRHS>
    auto operator+(const BinaryOp<LHS, RHS, Op>& l, const RRHS& r) -> decltype(add(l,r))
    { return add(l,r); }

    /*! \brief Right-handed addition operator for BinaryOp types. */
    template<class LHS, class RHS, detail::Op_t Op, class LLHS>
    auto operator+(const LLHS& l, const BinaryOp<LHS, RHS, Op>& r) -> decltype(add(l,r))
    { return add(l,r); }

    /*! \brief Addition operator for two BinaryOp types. */
    template<class LHS1, class RHS1, detail::Op_t Op1, class LHS2, class RHS2, detail::Op_t Op2>
    auto operator+(const BinaryOp<LHS1, RHS1, Op1>& l, const BinaryOp<LHS2, RHS2, Op2>& r) -> decltype(add(l,r))
    { return add(l,r); }

    /*! \brief Left-handed subtraction operator for BinaryOp types. */
    template<class LHS, class RHS, detail::Op_t Op, class RRHS>
    auto operator-(const BinaryOp<LHS, RHS, Op>& l, const RRHS& r) -> decltype(subtract(l,r))
    { return subtract(l,r); }

    /*! \brief Right-handed subtraction operator for BinaryOp types. */
    template<class LHS, class RHS, detail::Op_t Op, class LLHS>
    auto operator-(const LLHS& l, const BinaryOp<LHS, RHS, Op>& r) -> decltype(subtract(l,r))
    { return subtract(l,r); }

    /*! \brief Subtraction operator for two BinaryOp types. */
    template<class LHS1, class RHS1, detail::Op_t Op1, class LHS2, class RHS2, detail::Op_t Op2>
    auto operator-(const BinaryOp<LHS1, RHS1, Op1>& l, const BinaryOp<LHS2, RHS2, Op2>& r) -> decltype(subtract(l,r))
    { return subtract(l,r); }

    /*! \brief Expand addition BinaryOp types.

      If the classes have specialised operators for addition, then the
      decltype lookup will succeed and the addition is shunted to
      those classes. If not, this lookup will fail to expand the
      addition and it is instead carried out by the BinaryOp class.
    */
    template<class LHS, class RHS>
    auto expand(const BinaryOp<LHS, RHS, detail::ADD>& f) -> decltype(expand(f._l) + expand(f._r)) {
      return expand(f._l) + expand(f._r);
    }

    /*! \brief Expand multiplication BinaryOp types.

      If the classes have specialised operators for multiplication, then the
      decltype lookup will succeed and the addition is shunted to
      those classes. If not, this lookup will fail to expand the
      addition and it is instead carried out by the BinaryOp class.
    */
    template<class LHS, class RHS>
    auto expand(const BinaryOp<LHS, RHS, detail::MULTIPLY>& f) -> decltype(expand(f._l) * expand(f._r)) {
      return expand(f._l) * expand(f._r);
    }

    /*! \brief Derivatives of Addition operations.
    */
    template<class LHS, class RHS>
    auto derivative(const BinaryOp<LHS, RHS, detail::ADD>& f) -> decltype(derivative(f._l) + derivative(f._r))
    { return derivative(f._l) + derivative(f._r); }

    /*! \brief Derivatives of Multiplication operations.
    */
    template<class LHS, class RHS>
    auto derivative(const BinaryOp<LHS, RHS, detail::MULTIPLY>& f) -> decltype(add(multiply(derivative(f._l),f._r),multiply(derivative(f._r), f._l)))
    { return add(multiply(derivative(f._l),f._r),multiply(derivative(f._r), f._l)); }
    /*! \} */

    /*! \relates BinaryOp
      \name BinaryOp input/output operators
      \{
    */
    /*! \brief Writes a human-readable representation of the BinaryOp to the output stream. */
    template<class LHS, class RHS, detail::Op_t Op>
    inline std::ostream& operator<<(std::ostream& os, const BinaryOp<LHS, RHS, Op>& op) {
      os << "{" << op._l;
      switch (Op){
      case detail::ADD:      os << " + "; break;
      case detail::MULTIPLY: os << " * "; break;
      case detail::DIVIDE:   os << " / "; break;
      }
      os << op._r << "}";
      return os;
    }
    /*! \} */

    /*! \brief Symbolic representation of a (positive) power operator.
     */
    template<class Arg, size_t Power>
    struct PowerOp {
      Arg _arg;
      
      PowerOp(Arg a): _arg(a) {}
      
      /*! \brief Evaluate symbol at a value of x.
	
	This operator is only used if the result of evaluating the
	argument is an arithmetic type. If this is the case, the
	evaluation is passed to std::pow.
       */
      template<class R>
      auto operator()(const R& x) const -> typename std::enable_if<std::is_arithmetic<decltype(_arg(x))>::value, decltype(std::pow(_arg(x), Power))>::type {
	return std::pow(_arg(x), Power);
      }

      template<class R>
      auto operator()(const R& x) const -> typename std::enable_if<!std::is_arithmetic<decltype(_arg(x))>::value, decltype(std::pow(_arg(x), Power))>::type {
	static_assert(R(), "Not implemented yet");
      }
    };

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

    /*! \} */
    
  }
}
