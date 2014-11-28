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
    /*! \brief Symbolic representation of \f$f(x) *
      sin\left(g(x)\right) function.
     */
    template<class Prefactor, class Arg>
    class SinSymbol {
    public:
      Prefactor _prefactor;
      Arg _arg;

      SinSymbol(const Prefactor& p, const Arg& a):
	_prefactor(p),
	_arg(a)
      {}
      
      template<class Real>
      auto operator()(Real x) const ->decltype(_prefactor(x) * std::sin(_arg(x))) {
	return _prefactor(x) * std::sin(_arg(x));
      }

      template<class P>
      auto operator*(const P& p) -> SinSymbol<decltype(_prefactor * p), Arg> {
	return SinSymbol<decltype(_prefactor * p), Arg>(p * _prefactor, _arg);
      }
    };
    
    /*! \relates SinSymbol
      \brief Helper function for creating SinSymbols.
    */
    template<class A>
    SinSymbol<Polynomial<0>, A> Sin(const A& a) {
      return SinSymbol<Polynomial<0>, A>(Polynomial<0>{1.0}, a);
    }
    
    /*! \relates SinSymbol
      \brief Optimised Polynomial and SinSymbol multiplication.
      
      When multiplying a polynomial with a SinSymbol, it is preferred
      to prevent the distribution of the SinSymbol over the
      coefficients of the polynomial, but instead merge this with the
      Prefactor term.
    */
    template<class A, class P, class PReal, size_t POrder>
    auto operator*(const Polynomial<POrder, PReal>& p, const SinSymbol<P, A>& s) 
      -> SinSymbol<decltype(s._prefactor * p), A> {
      return SinSymbol<decltype(s._prefactor * p), A>(p * s._prefactor, s._arg);
    }

    /*! \relates SinSymbol
      \name SinSymbol input/output operations
      \{
    */
    /*! \brief Writes a human-readable representation of the SinSymbol
        to the output stream. 
    */
    template<class P, class A>
    inline std::ostream& operator<<(std::ostream& os, const SinSymbol<P, A>& s) {
      os << "(" << s._prefactor << ") * sin(" << s._arg << ")";
      return os;
    }
    /*! \} */
  }
}
