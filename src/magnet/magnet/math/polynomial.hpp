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

#include <stdexcept>
#include <ostream>
#include <array>

namespace magnet {
  namespace math {
    template<size_t Order, class Real = double> class Polynomial;
    namespace detail {
      constexpr size_t max_order(size_t N, size_t M) {
	return N > M ? N : M;
      }
    }

    template<size_t Order, class Real>
    class Polynomial : public std::array<Real, Order+1>
    {
      typedef std::array<Real, Order+1> Base;
    public:
      using Base::operator[];
      Polynomial() { Base::fill(Real()); };

      Polynomial(const Polynomial<Order+1, Real>& poly) {
#ifdef MAGNET_DEBUG
	if (poly[Order] != 0) 
	  M_throw() << "Trying to reduce the order of a polynomial with non-zero highest order coefficients!";
#endif
	std::copy(poly.begin(), poly.end()-1, Base::begin());
      };
      
      Polynomial(Real val) { Base::fill(Real()); Base::operator[](0) = val; }
      Polynomial(std::initializer_list<Real> _list) {
	if (_list.size() > Order+1)
	  throw std::length_error("initializer list too long");
      
	size_t i = 0; 
	auto it = _list.begin();
	for (; it != _list.end(); ++i, ++it)
	  Base::operator[](i) = *it;

	for (; i < Order+1; ++i)
	  Base::operator[](i) = 0.0;
      }

      template<size_t N, class Real2>
	Polynomial(const Polynomial<N, Real2>& poly) {
	static_assert(N <= Order, "Can only promote to higher order polynomials");
	size_t i(0);
	for (; i <= N; ++i)
	  Base::operator[](i) = poly[i];
	for (; i <= Order + 1; ++i)
	  Base::operator[](i) = Real();
      }

      Polynomial<Order> operator-() const {
	Polynomial<Order> retval;
	for (size_t i(0); i <= Order; ++i)
	  retval[i] = -Base::operator[](i);
	return retval;
      }

      /*! \brief Evaluate the polynomial at x. */
      Real eval(Real x) const {
	Real sum = Base::operator[](Order);
	for(int i = Order - 1; i >= 0; --i)
	  {
	    sum *= x;
	    sum += Base::operator[](i);
	  }
	return sum;
      }

      Real operator()(Real x) const {
	return eval(x);
      }
    };

    /* For all operations below we do not assume that we have a
       closure. For example, a vector multiplied by a vector is a
       scalar. If the polynomial coefficients are vectors, we must
       accomodate this.*/
    /*          Addition       */
    /* We assume that addition is commutative.*/
    template<class Real1, class Real2, size_t N>
    auto operator+(const Real1& r, const Polynomial<N,Real2>& poly)->Polynomial<N, decltype(poly[0] + r)>
    { return poly + r; }

    template<class Real1, class Real2, size_t N>
    auto operator+(const Polynomial<N,Real1>& poly, const Real2& r)->Polynomial<N, decltype(poly[0] + r)>
    {
      Polynomial<N, decltype(poly[0] + r)> retval(poly);
      retval[0] += r;
      return retval;
    }
    template<size_t M, size_t N, class Real1, class Real2>
    auto operator+(const Polynomial<M, Real1>& poly1, const Polynomial<N, Real2>& poly2)->Polynomial<detail::max_order(M, N), decltype(poly1[0] + poly2[0])>
    {
      Polynomial<detail::max_order(M, N), decltype(poly1[0] + poly2[0])> retval(poly1);
      for (size_t i(0); i <= N; ++i)
	retval[i] += poly2[i];
      return retval;
    }

    /*          Subtraction       */
    /* We assume that we have a closure under subtraction and that it is
       not commutative.*/
    template<class Real1, class Real2, size_t N>
    auto operator-(const Real1& r, const Polynomial<N, Real2>& poly)->Polynomial<N,decltype((-poly)[0]+r)>
    {
      Polynomial<N, decltype((-poly)[0]+r)> retval = -poly;
      retval[0] += r;
      return retval;  
    }
    template<class Real1, class Real2, size_t N>
    auto operator-(const Polynomial<N,Real1>& poly, const Real2& r)->Polynomial<N,decltype(poly[0]-r)>
    {
      Polynomial<N,decltype(poly[0]-r)> retval(poly);
      retval[0] -= r;
      return retval;
    }
    template<class Real1, class Real2, size_t M, size_t N>
    auto operator-(const Polynomial<M,Real1>& poly1, const Polynomial<N,Real2>& poly2)->Polynomial<detail::max_order(M, N),decltype(poly1[0]-poly2[0])>
    {
      Polynomial<detail::max_order(M, N),decltype(poly1[0]-poly2[0])> retval(poly1);
      for (size_t i(0); i <= N; ++i)
	retval[i] -= poly2[i];
      return retval;
    }

    /*          Multiplication       */
    /* We assume that multiplication commutes. We do not assume that we
       have a closure. i.e. multiplication of two vectors (a dot product)
       yields a scalar and a vector and a scalar yields a vector. We must
       therefore determine the new type of the coefficients of the
       polynomial.
    */
    template<class Real1, class Real2, size_t N>
    auto operator*(const Real1& r, const Polynomial<N, Real2>& poly) -> Polynomial<N, decltype(poly[0] * r)>
    { return poly * r; }
    template<class Real1, class Real2, size_t N>
    auto operator*(const Polynomial<N, Real1>& poly, const Real2& r) -> Polynomial<N, decltype(poly[0] * r)>
    {
      Polynomial<N, decltype(Real1() * Real2())> retval;
      for (size_t i(0); i <= N; ++i)
	retval[i] = poly[i] * r;
      return retval;
    }
    template<class Real1, class Real2, size_t M, size_t N>
    auto operator*(const Polynomial<M, Real1>& poly1, const Polynomial<N, Real2>& poly2) -> Polynomial<M + N, decltype(poly1[0] * poly2[0])>
    {
      Polynomial<M + N, decltype(poly1[0] * poly2[0])> retval;
      for (size_t i(0); i <= N+M; ++i)
	for (size_t j(i>N?i-N:0); (j <= i) && (j <=M); ++j)
	  retval[i] += poly1[j] * poly2[i-j];
      return retval;
    }

    template<class Real, size_t N>
    inline std::ostream& operator<<(std::ostream& os, const Polynomial<N, Real>& poly) {
      os << poly[0];
      for (size_t i(1); i <= N; ++i) {
	if (poly[i] == 0) continue;
	if (poly[i] == 1)
	  os << "+x";
	else if (poly[i] == -1)
	  os << "-x";
	else if (poly[i] > 0)
	  os << "+" << poly[i] << "*x";
	else
	  os << poly[i] << "*x";
	if (i > 1) os << "^" << i;
      }
      return os;
    }

    template<class Real, size_t N>
    inline Polynomial<N-1, Real> derivative(const Polynomial<N, Real>& f) {
      Polynomial<N-1, Real> retval;
      for (size_t i(0); i < N; ++i) {
	retval[i] = f[i+1] * (i+1);
      }
      return retval;
    }
  }
}
