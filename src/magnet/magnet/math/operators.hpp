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

#include <magnet/math/polynomial.hpp>

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

    /*! \brief Provides simplification of symbolic functions.

      The purpose of this function is to reduce the complexity of
      symbolic expressions to accelerate any successive
      evaluations. This should not change the 
      
      This should remove operations and simplify types wherever
      possible.
     */
    template<class T> T simplify(const T& f) { return f; }

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

    /*! \brief Simplify addition BinaryOp types.

      If the classes have specialised operators for addition, then the
      decltype lookup will succeed and the addition is shunted to
      those classes. If not, this lookup will fail to simplify the
      addition and it is instead carried out by the BinaryOp class.
    */
    template<class LHS, class RHS>
    auto simplify(const BinaryOp<LHS, RHS, detail::ADD>& f) -> decltype(simplify(f._l) + simplify(f._r)) {
      return simplify(f._l) + simplify(f._r);
    }

    /*! \brief Simplify multiplication BinaryOp types.

      If the classes have specialised operators for multiplication, then the
      decltype lookup will succeed and the addition is shunted to
      those classes. If not, this lookup will fail to simplify the
      addition and it is instead carried out by the BinaryOp class.
    */
    template<class LHS, class RHS>
    auto simplify(const BinaryOp<LHS, RHS, detail::MULTIPLY>& f) -> decltype(simplify(f._l) * simplify(f._r)) {
      return simplify(f._l) * simplify(f._r);
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
    /*! \brief Writes a human-readable representation of the Polynomial to the output stream. */
    template<class LHS, class RHS, detail::Op_t Op>
    inline std::ostream& operator<<(std::ostream& os, const BinaryOp<LHS, RHS, Op>& op) {
      os << " (" << op._l << ") ";
      switch (Op){
      case detail::ADD:      os << "+"; break;
      case detail::MULTIPLY: os << "*"; break;
      case detail::DIVIDE:   os << "/"; break;
      }
      os << " (" << op._r << ")";
      return os;
    }

    /*! \} */

  }
}
