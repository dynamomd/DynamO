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
#include <magnet/math/operators.hpp>
#include <magnet/containers/stack_vector.hpp>
#include <magnet/math/precision.hpp>
#include <magnet/exception.hpp>
#include <magnet/math/vector.hpp>
#include <boost/math/tools/roots.hpp>
#include <stdexcept>
#include <ostream>
#include <array>
#include <tuple>

namespace magnet {
  namespace math {
    template<size_t Order, class Real = double> class Polynomial;
    namespace detail {
      constexpr size_t max_order(size_t N, size_t M) {
 	return N > M ? N : M;
      }
    
      /*! \relates Polynomial 
	
 	\brief Type trait which determines if an operation
 	(multiplication, addition) can be distributed over a
 	polynomials coefficients.

 	This general form allows all operations between fundamental
 	aritmatic types.
      */
      template<class OpType, class PolyType>
      struct distribute_poly {
 	static const bool value = std::is_arithmetic<OpType>::value && std::is_arithmetic<PolyType>::value;
      };

      /*! \relates Polynomial 
	
 	\brief Type trait enabling use of std::complex as a Polynomial
 	coefficient with arithmetic types.
      */
      template<class OpType, class T>
      struct distribute_poly<OpType, std::complex<T> > {
 	static const bool value = std::is_arithmetic<OpType>::value;
      };

      /*! \relates Polynomial 
	
 	\brief Type trait enabling use of std::complex as an operation
 	with Polynomials with arithmetic type coefficients.
      */
      template<class T, class PolyType>
      struct distribute_poly<std::complex<T>, PolyType> {
 	static const bool value = std::is_arithmetic<PolyType>::value;
      };

      /*! \relates Polynomial 
	
 	\brief Type trait enabling use of std::complex as an operation
 	with Polynomials.
      */
      template<class T>
      struct distribute_poly<std::complex<T>, std::complex<T> > {
 	static const bool value = true;
      };

      /*! \relates Polynomial 
	
 	\brief Type trait enabling use of NVector types as an
 	operation with Polynomials.
      */
      template<class T, size_t N, class PolyType>
      struct distribute_poly<NVector<T, N>, PolyType> {
 	static const bool value = std::is_arithmetic<PolyType>::value;
      };

      /*! \relates Polynomial 
	
 	\brief Type trait enabling use of NVector types as an
	operation with Polynomials.
      */
      template<class T, size_t N, class OpType>
      struct distribute_poly<OpType, NVector<T, N> > {
 	static const bool value = std::is_arithmetic<OpType>::value;
      };

      /*! \relates Polynomial 
	
 	\brief Type trait enabling use of NVector types as an
 	operation with Polynomials.
      */
      template<class T, size_t N>
      struct distribute_poly<NVector<T, N>, NVector<T, N> > {
 	static const bool value = true;
      };
    }

    /*! \brief Representation of Polynomial with basic algebra operations.

      This class allows basic computer algebra to be performed with
      Polynomial equations.
    
      For example, the polynomial \f$f(x)=x^2 + 2\,x + 3\f$ can be created like so:
      \code{.cpp}
      Polynomial<1> x{0,1};
      auto f = x*x + 2*x +3;    
      \endcode
      And evaluated at the point \f$x=3\f$ like so:
      \code{.cpp}
      double val = f(3);    
      \endcode
    
      The class also functions with NVector coefficients.

      \tparam Order The order of the Polynomial.

      \tparam Real The type of the coefficients of the
      Polynomial. These may be NVector types.
    */
    template<size_t Order, class Real>
    class Polynomial : public std::array<Real, Order+1>
    {
      typedef std::array<Real, Order+1> Base;
    public:
      using Base::operator[];
      /*! \brief Default constructor.  

 	This initialises all Polynomial orders to be equivalent to
 	zero.
      */
      Polynomial() { Base::fill(Real()); }
    
      /*! \brief List initializer for simple Polynomial construction. 
	
 	This allows a polynomial to be constructed using just a list
 	of coefficients. E.g.
 	\code{.cpp}
 	Polynomial<1> f{0.5,1,2};
 	//f = 2 * x*x + x + 0.5;
 	\endcode
	
      */
      Polynomial(std::initializer_list<Real> _list) {
 	if (_list.size() > Order+1)
 	  throw std::length_error("initializer list too long");
    
 	size_t i = 0;
 	for (auto it = _list.begin(); it != _list.end(); ++i, ++it)
 	  Base::operator[](i) = *it;

 	for (; i <= Order; ++i)
 	  Base::operator[](i) = Real();
      }

      template<class InputIt>
      Polynomial(InputIt first, InputIt last) {
	size_t i = 0;
	auto it = first;
 	for (; (it != last) && (i <= Order); ++i, ++it)
 	  Base::operator[](i) = *it;
	
	if (it != last)
	  M_throw() << "Polynomial type too short to hold this range of coefficients";

 	for (; i <= Order; ++i)
 	  Base::operator[](i) = Real();
      }

      /*! \brief Constructor for constructing higher-order Polynomial
          types from lower order Polynomial types. 
      */
      template<size_t N, class Real2,
	       typename = typename std::enable_if<N <= Order>::type >
      Polynomial(const Polynomial<N, Real2>& poly) {
 	size_t i(0);
 	for (; i <= N; ++i)
 	  Base::operator[](i) = poly[i];
 	for (; i <= Order; ++i)
 	  Base::operator[](i) = Real();
      }

      /*! \brief Unary negative operator to change the sign of a Polynomial. */
      Polynomial<Order> operator-() const {
	Polynomial<Order> retval;
	for (size_t i(0); i <= Order; ++i)
	  retval[i] = -Base::operator[](i);
	return retval;
      }
    };
    
    /*! \brief Change the order of a Polynomial.

      This can be dangerous, as if the order of a polynomial is
      lowered high order terms are simply truncated.
    */
    template<size_t NewOrder, size_t Order, class Real>
    Polynomial<NewOrder, Real> change_order(const Polynomial<Order, Real>& f) {
      //The default constructor blanks Polynomial coefficients to zero
      Polynomial<NewOrder, Real> retval;
      //Just copy the coefficients which overlap between the new and old polynomial orders.
      std::copy(f.begin(), f.begin() + std::min(Order, NewOrder) + 1, retval.begin());
      return retval;
    };

    /*! \relates Polynomial 
      \name Polynomial set properties
      \{
    */

    /*! \brief Returns the empty product of a Polynomial.
      
      The empty product is a term whose multiplicative action is null
      (can be ignored).
    */
    template<size_t Order, class Real> 
    constexpr Polynomial<Order, Real> empty_product(const Polynomial<Order, Real>&)
    { return Polynomial<Order, Real>{1}; }

    /*! \brief Returns the empty sum of a Polynomial.
      
      The empty sum is a term whose multiplicative action is null (can
      be ignored).
    */
    template<size_t Order, class Real>
    constexpr Polynomial<Order, Real> empty_sum(const Polynomial<Order, Real>&)
    { return Polynomial<Order, Real>{}; }

    /*! \} */

    /*! \relates Polynomial 
      \name Polynomial algebraic operations
      
      For all operations below we do not assume that we have a
      closure. For example, a vector multiplied by a vector is a
      scalar therefore the * operator may change the returned type of
      the polynomial.
      \{
    */

    /*! \brief Evaluates a Polynomial expression at a given point.

      This function also specially handles the cases where
      \f$x=+\infty\f$ or \f$-\infty\f$ and returns the correct sign of
      infinity (if the polynomial has one non-zero coefficients of
      x). This behaviour is crucial as it is used in the evaluation of
      Sturm chains.
     */
    template<class Real, size_t Order, class Real2> 
    Real eval(const Polynomial<Order, Real>& f, const Real2& x)
    {
      //Handle the case where this is actually a constant and not a
      //Polynomial. This is free to evaluate now as Order is a
      //template parameter.
      if (Order == 0)
	return f[0];

      //Special cases for infinite values of x
      if (std::isinf(x)) {
	//Look through the Polynomial to find the highest order term
	//with a non-zero coefficient.
	for(size_t i = Order; i > 0; --i)
	  if (f[i] != 0) {
	    //Determine if this is an odd or even function of x
	    if (Order % 2)
	      //This is an odd function of x.
	      return (1 - 2 * (std::signbit(f[i]) ^ std::signbit(x))) * std::numeric_limits<Real>::infinity();
	    else
	      //This is an even function of x, the sign of x doesn't
	      //matter!
	      return (1 - 2 * std::signbit(f[i])) * std::numeric_limits<Real>::infinity();
	  };
	//All terms in x have zero as their coefficient
	return f[0];
      }

      Real sum = Real();
      for(size_t i = Order; i > 0; --i)
	sum = sum * x + f[i];
      sum = sum * x + f[0];
      return sum;
    }

    namespace {
      constexpr size_t factorial(size_t i) {
	return (i <= 1) ? 1 : i * factorial(i - 1);
      }

      constexpr size_t falling_factorial(size_t i, size_t end) {
	return (i <= end) ? 1 : i * falling_factorial(i - 1, end);
      }
    }

    /*! \brief Fast evaluation of multiple derivatives of a
        polynomial.
      
      This function is provided to allow derivatives to be evaluated
      without symbolically taking the derivative (causing a copy of
      the coefficients).
    */
    template<size_t D, size_t Order, class Real, class Real2>
    std::array<Real, D+1> eval_derivatives(const Polynomial<Order, Real>& f, const Real2& x)
    {
      std::array<Real, D+1> retval;
      retval.fill(Real());
      retval[0] = f[Order];
      for (size_t i(Order); i>0; i--) {
	for (size_t j = std::min(D, Order-(i-1)); j>0; j--)
	  retval[j] = retval[j] * x + retval[j-1];
	retval[0] = retval[0] * x + f[i-1];
      }

      Real cnst(1.0);
      for (size_t i(2); i <= D; i++) {
	cnst *= i;
	retval[i] *= cnst;
      }

      return retval;
    }

    /*! \brief Perform Euclidean division of a polynomial.
      
      Given two polynomials \f$f(x)\f$ and \f$g(x)\f$, the Euclidean
      division is a determination of the quotient polynomial
      \f$q(x)\f$ and the remainder polynomial \f$r(x)\f$, which satisfy
      
      \f[
      f(x) = g(x)\,q(x) + r(x)
      \f]

      If g(x) only consists of factors of f(x), the remainder will be
      zero. The algorithm used here is based on polynomial long
      division.

      If the division is by a monomial \f$g(x)=(x-a)\f$ where \f$a\f$
      is a root of \f$f(x)\f$, then the deflate_polynomial function
      should be used as it is more numerically stable and efficient.

      As g(x) may contain leading order coefficients which are zero,
      we cannot lower the order of the quotient polynomial returned.
    */
    template<size_t Order1, class Real, size_t Order2>
    std::tuple<Polynomial<Order1, Real>, Polynomial<Order2 - 1, Real> >
      euclidean_division(const Polynomial<Order1, Real>& f, const Polynomial<Order2, Real>& g)
    {
      static_assert(Order2 < Order1, "Cannot perform division when the order of the denominator is larger than the numerator using this routine");
      static_assert(Order2 > 0, "Constant division fails with these loops");
      typedef std::tuple<Polynomial<Order1, Real>, Polynomial<Order2 - 1, Real> > RetType;
      //If the leading term of g is zero, drop to a lower order
      //euclidean division.
      if (g[Order2] == 0)
	return RetType(euclidean_division(f, change_order<Order2 - 1>(g)));

      //The quotient and remainder.
      Polynomial<Order1, Real> r(f);
      Polynomial<Order1, Real> q;

      //Loop from the highest order coefficient of f, down to where we
      //have a polynomial one order lower than g.
      for (size_t k(Order1); k>=Order2; --k) {
	//Calculate the term on the quotient
	q[k-Order2] = r[k] / g[Order2];
	//Subtract this factor of other terms from the remainder
	for (size_t j(0); j <= Order2; j++)
	  r[k+j-Order2] -= q[k-Order2] * g[j];
      }

      return RetType(q, change_order<Order2 - 1>(r));
    }

    /*!  \cond Specializations
      \brief Specialisation for division by a constant.
     */
    template<size_t Order1, class Real>
    std::tuple<Polynomial<Order1, Real>, Polynomial<0, Real> >
      euclidean_division(const Polynomial<Order1, Real>& f, const Polynomial<0, Real>& g)
    {
      typedef std::tuple<Polynomial<Order1, Real>, Polynomial<0, Real> > RetType;
      if (g[0] == 0)
	return RetType(Polynomial<Order1, Real>{std::numeric_limits<Real>::infinity()}, 
		       Polynomial<0, Real>{Real()});

      return RetType(f * (1.0 / g[0]), Polynomial<0, Real>{Real()});
    }

    /*! \brief Right-handed addition operation on a Polynomial.
      
      This operator is only enabled if the type of the Polynomial
      coefficients and the type being added is marked as compatitble
      for distribution over the Polnomial coefficients. This is tested
      using detail::distribute_poly.
    */
    template<class Real1, class Real2, size_t N>
    auto operator+(const Real1& r, const Polynomial<N,Real2>& poly) -> typename std::enable_if<detail::distribute_poly<Real1, Real2>::value, Polynomial<N, decltype(poly[0] + r)> >::type
    { return poly + r; }

    /*!\brief Left-handed addition operator for Polynomials 

      This operator is only enabled if the type of the Polynomial
      coefficients and the type being added is marked as compatitble
      for distribution over the Polnomial coefficients. This is tested
      using detail::distribute_poly.
    */
    template<class Real1, class Real2, size_t N>
    auto operator+(const Polynomial<N,Real1>& poly, const Real2& r) -> typename std::enable_if<detail::distribute_poly<Real1, Real2>::value, Polynomial<N, decltype(poly[0] + r)> >::type
    {
      Polynomial<N, decltype(poly[0] + r)> retval(poly);
      retval[0] += r;
      return retval;
    }

    /*!\brief Addition operator for two Polynomial types. 
     */
    template<size_t M, size_t N, class Real1, class Real2>
    auto operator+(const Polynomial<M, Real1>& poly1, const Polynomial<N, Real2>& poly2)->Polynomial<detail::max_order(M, N), decltype(poly1[0] + poly2[0])>
    {
      Polynomial<detail::max_order(M, N), decltype(poly1[0] + poly2[0])> retval(poly1);
      for (size_t i(0); i <= N; ++i)
	retval[i] += poly2[i];
      return retval;
    }

    /*! \brief Right-handed subtraction operator for Polynomial types.
     
      This will reorder and convert the operation to a unary negation
      operator with an addition if the left-handed addition form
      exists.
    */
    template<class Real1, class Real2, size_t N>
    auto operator-(const Real1& r, const Polynomial<N, Real2>& poly) -> typename std::enable_if<detail::distribute_poly<Real1, Real2>::value, decltype((-poly) + r)>::type {
      return (-poly) + r;
    }
    
    /*! \brief Left-handed subtraction from a Polynomial type.

      This will convert the operation to a unary negation operator
      with an addition if the left-handed form exists.
    */
    template<class Real1, class Real2, size_t N>
    auto operator-(const Polynomial<N,Real1>& poly, const Real2& r) -> typename std::enable_if<detail::distribute_poly<Real1, Real2>::value, decltype(poly + (-r))>::type {
      return poly + (-r);
    }

    /*! \brief Subtraction between two Polynomial types. 
     */
    template<class Real1, class Real2, size_t M, size_t N>
    auto operator-(const Polynomial<M,Real1>& poly1, const Polynomial<N,Real2>& poly2) -> Polynomial<detail::max_order(M, N),decltype(poly1[0]-poly2[0])>
    {
      Polynomial<detail::max_order(M, N),decltype(poly1[0]-poly2[0])> retval(poly1);
      for (size_t i(0); i <= N; ++i)
	retval[i] -= poly2[i];
      return retval;
    }

    /*! \brief Right-handed multiplication operation on a Polynomial.
      
      This operator is only enabled if the type of the Polynomial
      coefficients and the type being added is marked as compatitble
      for distribution over the Polnomial coefficients. This is tested
      using detail::distribute_poly.
    */
    template<class Real1, class Real2, size_t N>
    auto operator*(const Real1& r, const Polynomial<N, Real2>& poly) -> typename std::enable_if<detail::distribute_poly<Real1, Real2>::value, Polynomial<N, decltype(poly[0] * r)> >::type
    { return poly * r; }

    /*! \brief Left-handed multiplication on a Polynomial.

      This operator is only enabled if the type of the Polynomial
      coefficients and the type being added is marked as compatitble
      for distribution over the Polnomial coefficients. This is tested
      using detail::distribute_poly.
    */
    template<class Real1, class Real2, size_t N>
    auto operator*(const Polynomial<N, Real1>& poly, const Real2& r) -> typename std::enable_if<detail::distribute_poly<Real1, Real2>::value, Polynomial<N, decltype(poly[0] * r)> >::type
    {
      Polynomial<N, decltype(Real1() * Real2())> retval;
      for (size_t i(0); i <= N; ++i)
	retval[i] = poly[i] * r;
      return retval;
    }

    /*! \brief Multiplication between two Polynomial types.
     */
    template<class Real1, class Real2, size_t M, size_t N>
    auto operator*(const Polynomial<M, Real1>& poly1, const Polynomial<N, Real2>& poly2) -> Polynomial<M + N, decltype(poly1[0] * poly2[0])>
    {
      Polynomial<M + N, decltype(poly1[0] * poly2[0])> retval;
      for (size_t i(0); i <= N+M; ++i)
	for (size_t j(i>N?i-N:0); (j <= i) && (j <=M); ++j)
	  retval[i] += poly1[j] * poly2[i-j];
      return retval;
    }

    /*! \brief Division of a Polynomial by a constant. */
    template<class Real1, class Real2, size_t N>
    auto operator/(const Polynomial<N, Real1>& poly, const Real2& r) -> typename std::enable_if<detail::distribute_poly<Real1, Real2>::value, Polynomial<N, decltype(poly[0] / r)> >::type
    {
      Polynomial<N, decltype(Real1() / Real2())> retval;
      for (size_t i(0); i <= N; ++i)
	retval[i] = poly[i] / r;
      return retval;
    }

    /*! \brief Enable reordering of Polynomial types. */
    template<class R1, size_t N1, class R2, size_t N2> 
    struct Reorder<Polynomial<N1, R1>, Polynomial<N2, R2> > {
      static const bool value = true;
    };

    /*! \brief Enable reordering of Polynomial types with arithmetic types. */
    template<class R1, size_t N1, class R2> 
    struct Reorder<Polynomial<N1, R1>, R2 > {
      static const bool value = std::is_arithmetic<R2>::value;
    };

    /*! \brief Enable reordering of Polynomial types with arithmetic types. */
    template<class R1, size_t N1, class R2> 
    struct Reorder<R2,Polynomial<N1, R1> > {
      static const bool value = std::is_arithmetic<R2>::value;
    };

    /*! \endcond */

    /*! \} */

    /*! \relates Polynomial 
      \name Polynomial calculus operations
      \{
    */

    /*! \brief Derivatives of Polynomial classes (constants).*/
    template<class Real, size_t N>
    inline Polynomial<N-1, Real> derivative(const Polynomial<N, Real>& f) {
      Polynomial<N-1, Real> retval;
      for (size_t i(0); i < N; ++i) {
	retval[i] = f[i+1] * (i+1);
      }
      return retval;
    }
    
    /*! \cond Specializations
      \brief Specialisation for derivatives of 0th order Polynomial classes (constants).
    */
    template<class Real>
    inline NullSymbol derivative(const Polynomial<0, Real>& f) {
      return NullSymbol();
    }

    /*! \endcond \} */

    /*! \relates Polynomial 
      \name Polynomial input/output operations
      \{
    */
    /*! \brief Writes a human-readable representation of the Polynomial to the output stream. */
    template<class Real, size_t N>
    inline std::ostream& operator<<(std::ostream& os, const Polynomial<N, Real>& poly) {
      std::ostringstream oss;
      oss.precision(os.precision());
      size_t terms = 0;
      for (size_t i(N); i != 0; --i) {
	if (poly[i] == empty_sum(poly[i])) continue;
	if (terms != 0)
	  oss << " + ";
	++terms;
	if (poly[i] != empty_product(poly[i]))
	  oss << poly[i] << " * ";
	oss << "x";
	if (i > 1)
	  oss << "^" << i;
      }
      if ((poly[0] != empty_sum(poly[0])) || (terms == 0)) {
	if (terms != 0)
	  oss << " + ";
	++terms;
	oss << poly[0];
      }
      if (terms > 1)
	os << "(" << oss.str() << ")";
      else
	os << oss.str();
      return os;
    }

    /*! \} */
    

    /*! \relates Polynomial 
      \name Polynomial roots
      \{
    */

    /*! \brief Factors out a root of a polynomial and returns a
      lower-order polynomial with the remaining roots.
	
      Given a polynomial, we can rearrange it in factored form like so
      \f[
      \sum_{i=0}^N a_i\,x^i =(x - r_1)\sum_{i=0}^{N-1} b_i\, x^{i}
      \f]
      
      where \f$r_1\f$ is a root of the polynomial. Equating terms on the
      LHS with terms on the RHS with equal powers of \f$x\f$, we have:

      \f[
      b_i=\frac{b_{i-1} - a_i}{r_1}\qquad \textrm{for}\ i\in[1,\,N-1]
      \f]

      This formula, known as backward deflation, can be used to
      calculate all coefficients using the starting point \f$b_0=-a_0
      / r_1\f$. This approach is not always stable (for example if the
      root is zero, or if \f$b_{i-1}\f$ has the same sign as \f$a_i\f$
      we might have catastrophic cancellation).
      
      An alternative "forward" iterative form may be found by
      substituting \f$i\to i+1\f$, which gives:

      \f[
      b_{i} = a_{i+1} + r_1\,b_{i+1} \qquad \textrm{for}\ i\in[0,\,N-2]
      \f]

      Again this approach may be used given the starting point
      \f$b_{N-1}=a_N\f$. However, it might also suffer from
      catastrophic cancellation if \f$a_{i+1}\f$ has the opposite sign
      to \f$r_1\,b_{i+1}\f$.

      It should be noted that Numerical Recipies states that "Forward
      deflation is stable if the largest absolute root is always
      divided out... backward deflation is stable if the smallest
      absolute root is always divided out". Unfortunately we do not
      know a priori the magnitude of the root being divided out.
      
      As both approaches may suffer from catastrophic cancellation, we
      decide to switch between them. If we catch the special
      root-is-zero case, we only must avoid catastrophic
      cancellation. This arises if two non-zero terms are subtracted
      from each other (i.e., for the first approach this happens if
      \f$a_{i+1}\f$ and \f$r_1\,b_{i+1}\f$ are non-zero and have
      opposite sign). We could use this to monitor the bits of
      precision "lost" as we calculate from each end and select a
      point between the two methods where accuracy is highest, but
      this would require a more detailed analysis of the error. A
      simple approach is to solve from both ends of the polynomial at
      the same time and only actually accept whichever has the lowest 
      catastrophic cancellation accuracy in terms of bits.

      \param f The Polynomial to factor a root out of.
      \param root The root to remove.
    */
    template<size_t Order, class Real>
    inline Polynomial<Order-1, Real> deflate_polynomial(const Polynomial<Order, Real>& a, const double root) {
      Polynomial<Order-1, Real> b;
      
      //Check for the simple case where root==0. If this is the case,
      //then the deflated polynomial is actually just a divide by x.
      if (root == 0) {
	std::copy(a.begin()+1, a.end(), b.begin());
	return b;
      }
	
      //Calculate the highest and lowest order coefficients using
      //these stable approaches
      b[Order-1] = a[Order];
      b[0] = - a[0] / root;

      size_t i_t = Order-2;
      size_t i_b = 1;
      while (i_t >= i_b) {
	const Real d = root * b[i_t + 1];	
	if (subtraction_precision(b[i_b], a[i_b]) > addition_precision(a[i_t+1], d)) {
	  b[i_b] = (b[i_b-1] - a[i_b]) / root;
	  ++i_b;
	} else {
	  b[i_t] = a[i_t+1] + d;
	  --i_t;
	}
      }
      return b;
    }
    
    /*! \brief Returns a polynomial \f$g(x)=f(x+t)\f$.
      
      Given a polynomial \f$f(x)\f$:
      \f[
      f(x) = \sum_{i=0}^N a_i\,x^i
      \f]
      
      We wish to determine the coefficients of a polynomial
      \f$g(x)=f(x+t)\f$:

      \f[
      g(x) = \sum_{i=0}^N b_i\,x^i
      \f]

      We can define \f$g(x)\f$ by taking a Taylor expansion of
      \f$f(x)\f$ about the point \f$t\f$, we have:
      
      \f[
      g(x) = f(t+x) = \sum_{i=0}^n \frac{f^i(t)}{i!}x^i
      \f]
      
      where \f$f^i(x)\f$ is the \f$i\f$th derivative of \f$f(x)\f$ and
      \f$n\f$ is the order of the polynomial \f$f\f$. Each coefficient
      of \f$g\f$ is then given by:
      
      \f[
      b_i = \frac{f^i(t)}{i!}
      \f]

      Here, we then use a modified version of the \ref
      eval_derivatives function to actually calculate the derivatives
      while avoiding the factorial term.
     */
    template<size_t Order, class Real>
    inline Polynomial<Order, Real> shift_function(const Polynomial<Order, Real>& f, const double t) {
      //Check for the simple case where t == 0, nothing to be done
      if (t == 0) return f;

      Polynomial<Order, Real> retval;
      retval.fill(Real());
      retval[0] = f[Order];
      for (size_t i(Order); i>0; i--) {
	for (size_t j = Order-(i-1); j>0; j--)
	  retval[j] = retval[j] * t + retval[j-1];
	retval[0] = retval[0] * t + f[i-1];
      }
      return retval;
    }

    /*! \brief Returns a polynomial \f$g(x)=f(x+1)\f$.

      This is an optimised \ref shift_polynomial operation where the
      shift is unity. See \ref shift_polynomial for more
      implementation details
     */
    template<size_t Order, class Real>
    inline Polynomial<Order, Real> taylor_shift(const Polynomial<Order, Real>& f) {
      Polynomial<Order, Real> retval;
      retval.fill(Real());
      retval[0] = f[Order];
      for (size_t i(Order); i>0; i--) {
	for (size_t j = Order-(i-1); j>0; j--)
	  retval[j] += retval[j-1];
	retval[0] += f[i-1];
      }
      return retval;
    }

    /*! \brief Returns a polynomial \f$g(x)=(x+1)^{N}
        f\left(\frac{1}{x+1}\right)\f$.

	
	The creation of \f$p(x)\f$, given by
	\f[
	p(x) = \left(x+1\right)^d\,f\left(\frac{1}{x+1}\right)
	\f]

	where \f$d\f$ is the order of the polynomial \f$f(x)\f$, is
	often carried out while locating roots of a Polynomial. It has
	the useful property of generating a polynomial which has the
	same number of roots of \f$f(x)\f$ in the range
	\f$x\in[0,\,1]\f$, but now in the range
	\f$x\in[\infty,\,0]\f$. Therefore, small sections of a
	polynomial may be inspected for roots using scaling, shifting,
	and this transformation. It is a special case of a Mobius
	transformation of the Polynomial.

	Creation of \f$p(x)\f$ may be carried out in two steps. First,
	the following equation is generated:

	\f[
	p_1(x) = x^d\,f\left(\frac{1}{x}\right)
	\f]
	
	This operation may be performed simply by reversing the order
	of the coefficient array in the Polynomial. Then, a \ref
	taylor_shift is applied to complete the transformation:

	\f[
	p_2(x) = p_1\left(x+1\right)
	\f]

	This entire operation is performed using an optimised version
	of the \ref taylor_shift algorithm, which itself is an
	optimized \ref shift_polynomial function.
     */
    template<size_t Order, class Real>
    inline Polynomial<Order, Real> invert_taylor_shift(const Polynomial<Order, Real>& f) {
      Polynomial<Order, Real> retval;
      retval.fill(Real());
      retval[0] = f[0];
      for (size_t i(Order); i>0; i--) {
	for (size_t j = Order-(i-1); j>0; j--)
	  retval[j] += retval[j-1];
	retval[0] += f[Order - (i-1)];
      }
      return retval;
    }

    /*! \brief Returns a polynomial \f$g(x)=f\left(-x\right)\f$.
     */
    template<size_t Order, class Real>
    inline Polynomial<Order, Real> reflect_poly(const Polynomial<Order, Real>& f) {
      Polynomial<Order, Real> g(f);
      for (size_t i(1); i <= Order; i+=2)
	g[i] = -g[i];
      return g;
    }

    /*! \brief Returns a polynomial \f$g(x)=f\left(a\,x\right)\f$
        where \f$a\f$ is a scaling factor.
     */
    template<size_t Order, class Real>
    inline Polynomial<Order, Real> scale_poly(const Polynomial<Order, Real>& f, const Real& a) {
      Polynomial<Order, Real> g(f);
      Real factor = 1;
      for (size_t i(1); i <= Order; ++i)
	g[i] *= (factor *= a);
      return g;
    }

    /*!  \cond Specializations
      \brief Specialisation for no roots of a 0th order Polynomial.
      \param f The Polynomial to evaluate.
    */
    inline containers::StackVector<double, 0> solve_roots(const Polynomial<0, double>& f) {
      return containers::StackVector<double, 0>();
    }

    /*! \brief Calculate the single real root (if it exists) of a 1st
        order Polynomial.
      
      \param f The Polynomial to evaluate.
    */
    inline containers::StackVector<double, 1> solve_roots(const Polynomial<1, double>& f) {
      containers::StackVector<double, 1> roots;
      if (f[1] != 0)
	roots.push_back(-f[0] / f[1]);
      return roots;
    }

    /*! \brief Calcualte the real roots of a 2nd order Polynomial
        using radicals.
      
      \param f The Polynomial to evaluate.
    */
    inline containers::StackVector<double, 2> solve_roots(Polynomial<2, double> f) {
      //If this is actually a linear polynomial, drop down to that solver.
      if (f[2] == 0) 
	return solve_roots(change_order<1>(f));
      
      //Scale the constant of x^2 to 1
      f = f / f[2];

      if (f[0] == 0)
	//There is no constant term, so we actually have x^2 + f[1] * x = 0
	return containers::StackVector<double, 2>{-f[1]};
      
      static const double maxSqrt = std::sqrt(std::numeric_limits<double>::max());
      if ((f[1] > maxSqrt) || (f[1] < -maxSqrt)) {
	//arg contains f[1]*f[1], so it will overflow. In this case we
	//can approximate the equation as x^2 + a x = 0 to solve for
	//one root, and use root2 = f[0]/root1 to find the second
	//root. This should work even with large constant values
	return containers::StackVector<double, 2>{-f[1], -f[0] / f[1]};
      }

      const double arg = f[1] * f[1] - 4 * f[0];

      //Test if there are real roots   
      if (arg < 0)
	return containers::StackVector<double, 2>();

      //Test if there is a double root
      if (arg == 0)
	return containers::StackVector<double, 2>{-f[1] * 0.5};

      //Return both roots
      const double root1 = -(f[1] + std::copysign(std::sqrt(arg), f[1])) * 0.5;
      const double root2 = f[0] / root1;
      return containers::StackVector<double, 2>{root1, root2};
    }


    namespace {
      /*! \brief Deflate a Polynomial and solves for the remaining
	roots.

	This routine also polishes the root to try to improve accuracy;
	however, do not assume this function will behave well with
	inaccurate roots.
      */
      template<size_t Order, class Real>
      inline containers::StackVector<Real, Order> deflate_and_solve_polynomial(const Polynomial<Order, Real>& f, Real root) {
	containers::StackVector<Real, Order> roots = solve_roots(deflate_polynomial(f, root));
	roots.push_back(root);
	return roots;
      }

      /*! \brief Uses a quadratic scheme to polish up a root.
       */
      inline void cubicNewtonRootPolish(const Polynomial<3, double>& f, double& root)
      {
	//Stored for comparison later
	double error = eval(f, root);
	const size_t maxiterations = 4;
	for (size_t it = 0; (it < maxiterations) && (error != 0); ++it)
	  {
	    //Calculate the 1st and 2nd derivatives
	    double deriv = (3 * root + 2 * f[2]) * root + f[1];
	    double dderiv = 6 * root + 2 * f[2];
	    
	    //Try a quadratic scheme to improve the root
	    auto roots = solve_roots(Polynomial<2, double>{error, deriv, 0.5 * dderiv});
	    if (roots.size() == 2)
	      root += (std::abs(roots[0]) < std::abs(roots[1])) ? roots[0] : roots[1];
	    else
	      { //Switch to a linear scheme if the quadratic fails,
		//but if the derivative is zero then just accept this
		//is the closest we will get.
		if (deriv == 0) return;
		root -= error / deriv;
	      }
	    error = eval(f, root);
	  }
      }
    }
    
    /*! \brief Calculate the real roots of a 3rd order Polynomial
        using radicals.
      
      \param f The Polynomial to evaluate.
    */
    inline containers::StackVector<double, 3> solve_roots(const Polynomial<3, double>& f_original) {
      //Ensure this is actually a third order polynomial
      if (f_original[3] == 0)
	return solve_roots(change_order<2>(f_original));
      
      if (f_original[0] == 0)
	//If the constant is zero, one root is x=0.  We can divide
	//by x and solve the remaining quadratic
	return deflate_and_solve_polynomial(f_original, 0.0);

      //Convert to a cubic with a unity high-order coefficient
      auto f = f_original / f_original[3];
      
      if ((f[2] == 0) && (f[1] == 0))
	//Special case where f(x) = x^3 + f[0]
	return containers::StackVector<double, 3>{std::cbrt(-f[0])};

      static const double maxSqrt = std::sqrt(std::numeric_limits<double>::max());
      
      if ((f[2] > maxSqrt) || (f[2] < -maxSqrt))
	//The equation is limiting to x^3 + f[2] * x^2 == 0. Use
	//this to estimate the location of one root, polish it up,
	//then deflate the polynomial and solve the quadratic.
	return deflate_and_solve_polynomial(f, -f[2]);

      //NOT SURE THESE RANGE TESTS ARE BENEFICIAL
      if (f[1] > maxSqrt)
	//Special case, if f[1] is large (and f[2] is not) the root is
	//near -f[0] / f[1], the x^3 term is negligble, and all other terms
	//cancel.
	return deflate_and_solve_polynomial(f, -f[0] / f[1]);

      if (f[1] < -maxSqrt)
	//Special case, equation is approximated as x^3 + q x == 0
	return deflate_and_solve_polynomial(f, -std::sqrt(-f[1]));

      if ((f[0] > maxSqrt) || (f[0] < -maxSqrt))
	//Another special case where equation is approximated asf(x)= x^3 +f[0]
	return deflate_and_solve_polynomial(f, -std::cbrt(f[0]));

      const double v = f[0] + (2.0 * f[2] * f[2] / 9.0 - f[1]) * (f[2] / 3.0);

      if ((v > maxSqrt) || (v < -maxSqrt))
	return deflate_and_solve_polynomial(f, -f[2]);
      
      const double uo3 = f[1] / 3.0 - f[2] * f[2] / 9.0;
      const double u2o3 = uo3 + uo3;
      
      if ((u2o3 > maxSqrt) || (u2o3 < -maxSqrt))
	{
	  if (f[2]==0)
	    {
	      if (f[1] > 0)
		return deflate_and_solve_polynomial(f, -f[0] / f[1]);
	      
	      if (f[1] < 0)
		return deflate_and_solve_polynomial(f, -std::sqrt(-f[1]));
		
	      return deflate_and_solve_polynomial(f, 0.0);
	    }

	  return deflate_and_solve_polynomial(f, -f[1] / f[2]);
	}

      const double uo3sq4 = u2o3 * u2o3;
      if (uo3sq4 > maxSqrt)
	{
	  if (f[2] == 0)
	    {
	      if (f[1] > 0)
		return deflate_and_solve_polynomial(f, -f[0] / f[1]);

	      if (f[1] < 0)
		return deflate_and_solve_polynomial(f, -std::sqrt(-f[1]));

	      return deflate_and_solve_polynomial(f, 0.0);
	    }

	  return deflate_and_solve_polynomial(f, -f[1] / f[2]);
	}

      const double j = (uo3sq4 * uo3) + v * v;

      if (j > 0) 
	{//Only one root (but this test can be wrong due to a
	  //catastrophic cancellation in j 
	  //(i.e. (uo3sq4 * uo3) == v * v)

	  const double w = std::sqrt(j);
	  double root;
	  if (v < 0)
	    root = std::cbrt(0.5*(w-v)) - (uo3) * std::cbrt(2.0 / (w-v)) - f[2] / 3.0;
	  else
	    root = uo3 * std::cbrt(2.0 / (w+v)) - std::cbrt(0.5*(w+v)) - f[2] / 3.0;

	  return deflate_and_solve_polynomial(f, root);
	}
  
      if (uo3 >= 0)
	//Multiple root detected
	return containers::StackVector<double, 3>{std::cbrt(v) - f[2] / 3.0};

      const double muo3 = - uo3;
      double s = 0;
      if (muo3 > 0)
	{
	  s = std::sqrt(muo3);
	  if (f[2] > 0) s = -s;
	}
      
      const double scube = s * muo3;
      if (scube == 0)
	return containers::StackVector<double, 3>{ -f[2] / 3.0 };
      
      const double t = - v / (scube + scube);
      const double k = std::acos(t) / 3.0;
      const double cosk = std::cos(k);
      
      containers::StackVector<double, 3> roots{ (s + s) * cosk - f[2] / 3.0 };
      
      const double sinsqk = 1.0 - cosk * cosk;
      if (sinsqk < 0)
	return roots;

      double rt3sink = std::sqrt(3.0) * std::sqrt(sinsqk);
      roots.push_back(s * (-cosk + rt3sink) - f[2] / 3.0);
      roots.push_back(s * (-cosk - rt3sink) - f[2] / 3.0);

      cubicNewtonRootPolish(f, roots[0]);
      cubicNewtonRootPolish(f, roots[1]);
      cubicNewtonRootPolish(f, roots[2]);

      return roots;
    }
    /*! \endcond */

    namespace detail {
      /*! \brief Calculates the negative of the remainder of the
          division of \f$f(x)\f$ by \f$g(x)\f$.
       */
      template<size_t Order, class Real>
      Polynomial<Order-2, Real> 
	mrem(const Polynomial<Order, Real>& f, const Polynomial<Order-1, Real>&g)
      {
	Polynomial<Order-2, Real> rem;
	std::tie(std::ignore, rem) = euclidean_division(f, g);
	return -rem;
      }

      /*! \brief A collection of Polynomials which form a Sturm chain. 
       */
      template<size_t Order, class Real>
      struct SturmChain {
	/*! \brief Constructor if is the first Polynomial in the
            chain. 
	*/
	SturmChain(const Polynomial<Order, Real>& p_n):
	  _p_n(p_n), _p_nminus1(p_n, derivative(p_n))
	{}

	/*! \brief Constructor if is an intermediate Polynomial in the
            chain.
	*/
	SturmChain(const Polynomial<Order+1, Real>& p_nplus1, const Polynomial<Order, Real>& p_n):
	  _p_n(p_n), _p_nminus1(p_n, mrem(p_nplus1, p_n))
	{}

	Polynomial<Order, Real> _p_n;
	SturmChain<Order-1, Real> _p_nminus1;
	
	/*! Accessor function for the ith Polynomial in the Sturm
            chain.

	    This promotes the order of the Sturm chain polynomial to
	    the original order of the Polynomial as this is done at
	    runtime.
	*/
	Polynomial<Order, Real> get(size_t i) const {
	  if (i == 0)
	    return _p_n;
	  else
	    return _p_nminus1.get(i-1); 
	}
	
	/*! Count the number of sign changes in the Sturm chain
	  evaluated at \f$x\f$.
	  
	  This actually uses a helper function sign_change_helper to
	  carry out the calculation.
	*/
	template<class Real2>
	size_t sign_changes(const Real2& x) const {
	  return sign_change_helper(0, x);
	}

	template<class Real2>
	size_t roots(const Real2& a, const Real2& b) const {
#ifdef MAGNET_DEBUG
	  if (a > b)
	    M_throw() << "a<b!";
#endif
	  return sign_changes(a) - sign_changes(b);
	}
	
	
	/*! \brief Helper function for calculating the sign changes in
            the Sturm chain.

	    These functions use -1, 0, and +1 to denote the sign of an
	    evaluation of a polynomial. The sign of the previous
	    Polynomial in the Strum chain is given as last_sign. If
	    this is zero, then there has been no sign so far (all
	    previous polynomials were zero or this is the first
	    polynomial in the chain).
	 */
	template<class Real2>
	size_t sign_change_helper(const int last_sign, const Real2& x) const {
	  const Real currentx = eval(_p_n, x);
	  const int current_sign = (currentx != 0) * (1 - 2 * std::signbit(currentx));
	  const bool sign_change = (current_sign != 0) && (last_sign != 0) && (current_sign != last_sign);

	  const int next_sign = (current_sign != 0) ? current_sign : last_sign;
	  return _p_nminus1.sign_change_helper(next_sign, x) + sign_change;
	}

	void output_helper(std::ostream& os, const size_t max_order) const {
	  os << ",\n           p_" <<  max_order - Order << "=" << _p_n;
	  _p_nminus1.output_helper(os, max_order);
	}
      };

      /*! \brief Specialisation for a container holding the last Sturm
          chain Polynomial.
      */
      template<class Real>
      struct SturmChain<0, Real> {
	/*! \brief Constructor  is the first and last Polynomial in
            the chain.
	*/
	SturmChain(const Polynomial<0, Real>& p_n):
	  _p_n(p_n)
	{}

	/*! \brief Constructor if this is the last Polynomial in the
            chain.
	*/
	SturmChain(const Polynomial<1, Real>& p_nplus1, const Polynomial<0, Real>& p_n):
	  _p_n(p_n) {}

	Polynomial<0, Real> get(size_t i) const {
	  if (i == 0)
	    return _p_n;
	  return Polynomial<0,Real>{};
	}

	template<class Real2>
	size_t sign_changes(const Real2& x) const {
	  return 0;
	}

	Polynomial<0, Real> _p_n;

	template<class Real2>
	size_t sign_change_helper(const int last_sign, const Real2& x) const {
	  const Real currentx = eval(_p_n, x);
	  const int current_sign = (currentx != 0) * (1 - 2 * std::signbit(currentx));
	  const bool sign_change = (current_sign != 0) && (last_sign != 0) && (current_sign != last_sign);
	  return sign_change;
	}
	
	void output_helper(std::ostream& os, const size_t max_order) const {
	  os << ",\n           p_" <<  max_order << "=" << _p_n;
	}
      };
      
      template<size_t Order, class Real>
      std::ostream& operator<<(std::ostream& os, const SturmChain<Order, Real>& c) {
	os << "SturmChain{p_0=" << c._p_n;
	c._p_nminus1.output_helper(os, Order);
	os << "}";
	return os;
      }
    }

    /*! \brief Helper function to generate a SturmChain from a
        Polynomial.
      
      The actual calculation, storage, and evaluation of the Sturm
      chain is done by the detail::SturmChain type.

      The Sturm chain is a sequence of polynomials \f$p_0(x)\f$,
      \f$p_1(x)\f$, \f$p_2(x)\f$, \f$\ldots\f$, \f$p_n(x)\f$ generated
      from a single polynomial \f$f(x)\f$ of order \f$n\f$. The first
      two Sturm chain polynomials are given as
      
      \f{eqnarray*}{
      p_0(x) &=& f(x)\\
      p_1(x) &=& f'(x)
      \f}

      All higher polynomials are evaluated like so:

      \f{eqnarray*}{
      p_n(x) &=& -\mathrm{rem}(p_{n+2}, p_{n+1})
      \f}
      
      where \f$\mathrm{rem}(p_{n+2},\,p_{n+1})\f$ returns the
      remainder polynomial from a \ref euclidean_division of
      \f$p_{n+2}\f$ over \f$p_{n+1}\f$. This sequence terminates at
      \f$p_N\f$, where \f$N\f$ is the order of the original
      Polynomial, \f$f(x)\f$.
      
      The interesting property of this chain is that it allows a
      calculation of the count of distinct real roots of \f$f(x)\f$
      within a certain range. If we evaluate the Sturm chain of
      \f$f(x)\f$ at a point \f$\xi\f$ and count the number of changes
      in sign (ignoring zeros) in the sequence:
      
      \f[
      p_0(\xi),\,p_1(\xi),\,\ldots p_n(\xi)
      \f]
      
      and define this as \f$\sigma(\xi)\f$. Then, given two real
      numbers \f$a<b\f$, the number of distinct roots of \f$f(x)\f$ in
      \f$(a,b]\f$ is given by \f$\sigma(a)-\sigma(b)\f$.

      This gives a method to calculate the exact number of real and
      distinct roots in a region, which can then be used in a
      bisection routine to bound individual distinct roots. The Sturm
      sequence can also be easily evaluated at infinite bounds
      \f$(-\infty,+\infty)\f$ to determine the total number of real
      roots.

      Although this method allows the construction of an
      arbitrary-order Polynomial root finder through bisection,
      inexact methods for computing the number of roots in a region
      (such as bodans_01_test) are preferred as they are more
      computationally efficient.
    */
    template<size_t Order, class Real>
    detail::SturmChain<Order, Real> sturm_chain(const Polynomial<Order, Real>& f) {
      return detail::SturmChain<Order, Real>(f);
    }

    /*! \brief Calculates an upper bound estimate for the number of
      positive real roots of a Polynomial (including multiples).
      
      Descarte's rule of signs states that the number of positive real
      roots for a single-variable real-coefficient Polynomial is less
      than or equal to the number of sign changes between consecutive
      non-zero coefficients in the Polynomial. When the actual root
      count is less, it is less by an even number. Therefore, the
      values 0 or 1 are exact.
     */
    template<size_t Order, class Real>
    size_t descartes_rule_of_signs(const Polynomial<Order, Real>& f) {
      //Count the sign changes
      size_t sign_changes(0);
      int last_sign = 0;      
      for (size_t i(0); i <= Order; ++i) {
	const int current_sign = (f[i] != 0) * (1 - 2 * std::signbit(f[i]));
	sign_changes += (current_sign != 0) && (last_sign != 0) && (current_sign != last_sign);
	last_sign = (current_sign != 0) ? current_sign : last_sign;
      }
      return sign_changes;
    }

    /*! \brief Budan's upper bound estimate for the number of real
        roots in a Polynomial over the range \f$(0,\,1)\f$.

	Budan's test is actually just Descarte's test, but on a
	transformed Polynomial \f$p(x)\f$, which is related to
	\f$f(x)\f$ as follows:
	
	\f[
	p(x) = \left(x+1\right)^d\,f\left(\frac{1}{x+1}\right)
	\f]

	where \f$d\f$ is the order of the polynomial \f$f(x)\f$. The
	roots of \f$f(x)\f$ in the range \f$[0,1]\f$ are mapped over
	the range \f$[0,\infty]\f$ of \f$p(x)\f$. This allows
	Descarte's rule of signs to be applied to a limited range of
	the original polynomial.

	The actual transformation to \f$p(x)\f$ is carried out using
	the \ref invert_taylor_shift function, before this is passed
	to \ref descartes_rule_of_signs.

	\return An upper bound on the number of real roots in the
	interval \f$(0,\,1)\f$.
     */
    template<size_t Order, class Real>
    size_t budan_01_test(const Polynomial<Order, Real>& f) {
      return descartes_rule_of_signs(invert_taylor_shift(f));
    }

    /*! \brief Alesina-Galuzzi upper bound estimate for the numer of
        real roots in a Polynomial over a specified range
        \f$(a,\,b)\f$.

	This a generalisation of Budan's 01 test (see \ref
	budan_01_test) and is implemented that way. The polynomial is
	shifted so that \f$0\f$ corresponds to \f$a\f$, then scaled so
	that \f$x=1\f$ corresponds to \f$b\f$, before Budan's test is
	called on the transformed Polynomial.
    */
    template<size_t Order, class Real>
    size_t alesina_galuzzi_test(const Polynomial<Order, Real>& f, const Real& a, const Real& b) {
      return budan_01_test(scale_poly(shift_function(f, a), b - a));
    }
    
    /*! \brief Local-max Quadratic upper bound estimate for the value
        of the real roots of a Polynomial.

	This function is adapted from the thesis "Upper bounds on the
	values of the positive roots of polynomials" by Panagiotis
	S. Vigklas (2010). The main change is to generalise to
	arbitrary sign on the highest order coefficient, and to allow
	high-order coefficients with zero values.
     */
    template<class Real, size_t Order>
    Real LMQ_upper_bound(const Polynomial<Order, Real>& f) {
      std::array<size_t, Order+1> times_used;
      times_used.fill(1);
      Real ub = Real();

      size_t real_order = Order;
      while ((real_order > 0) && (f[real_order] == 0))
	--real_order;

      for (int m(real_order-1); m >= 0; --m)
	if (std::signbit(f[m]) != std::signbit(f[real_order])) {
	  Real tempub = std::numeric_limits<Real>::infinity();
	  for (int k(real_order); k > m; --k)
	    if (std::signbit(f[k]) != std::signbit(f[m]))
	      {
		Real temp = std::pow(-(1 << times_used[k]) * f[m] / f[k], 1.0 / (k - m));
		++times_used[k];
		tempub = std::min(temp, tempub);
	      }
	  ub = std::max(tempub, ub);
	}
      return ub;
    }

    /*! \brief Local-max Quadratic lower bound estimate for the real
      roots of a Polynomial.
      
      Given a Polynomial \f$f(x)\f$, this function performs the
      transformation:
      
      \f[
      g(x) = x^n\,f(\frac{1}{x})
      \f]
      
      Now the upper bound on the real roots of \f$g(x)\f$ are the
      inverse of the lower bound on the real roots of \f$f(x)\f$. The
      transformation is computationally equivalent to reversing the
      coefficient array of the polynomial \f$f(x)\f$.
    */
    template<class Real, size_t Order>
    Real LMQ_lower_bound(const Polynomial<Order, Real>& f) {
      return 1.0 / LMQ_upper_bound(Polynomial<Order, Real>(f.rbegin(), f.rend()));
    }

    /*! \cond Specializations

      \brief Specialisation of Local-max Quadratic upper-bound
      estimation for real roots of a Polynomial, where the Polynomial
      is a constant.
    */
    template<class Real>
    Real LMQ_upper_bound(const Polynomial<0, Real>& f) {
      return 0;
    }

    /*!
      \brief Specialisation of Local-max Quadratic
      lower-bound estimation for real roots of a Polynomial, where the
      Polynomial is a constant.
    */
    template<class Real>
    Real LMQ_lower_bound(const Polynomial<0, Real>& f) {
      return HUGE_VAL;
    }

    /*! \brief Calculate interval bounds on all of the positive real
      roots between \f$(0,1)\f$ of a squarefree Polynomial.
      
      This function uses the VCA algorithm to bound the roots. It
      assumes that the polynomial has a non-zero constant term and
      leading order coefficient term.
    */
    template<size_t Order, class Real>
    containers::StackVector<std::pair<Real,Real>, Order> 
    VCA_real_root_bounds_worker(const Polynomial<Order, Real>& f) {
      //Test how many roots are in the range (0,1)
      switch (budan_01_test(f)) {
      case 0:
	//No roots, return empty
	return containers::StackVector<std::pair<Real,Real>, Order>();
      case 1:
	//One root! the bound is (0,1)
	return containers::StackVector<std::pair<Real,Real>, Order>{std::make_pair(Real(0), Real(1))};
      default:
	//Possibly multiple roots, so divide the range and recursively
	//explore it.

	//Scale the polynomial so that all roots lie in the range
	//(0,2). We actually scale by 2^Order to make the division by 2 a multiply.
	//
	//We actually generate this function
	//p1(x) = 2^Order f(x/2)
	Polynomial<Order, Real> p1(f);
	for (size_t i(0); i <= Order; ++i)
	  p1[i] *= (1 << (Order-i)); //This gives (2^Order) / (2^i)

	//Perform a Taylor shift p2(x) = p1(x+1). This gives 
	//
	//p2(x) = 2^Order f(x/2 + 0.5) 
	//
	//in terms of the original function, f(x).
	const Polynomial<Order, Real> p2 = taylor_shift(p1);

	//Now that we have two polynomials, each of which is scaled so
	//the roots of interest lie in (0,1). Search them both and
	//combine the results.

	//Detect and scale the first range's roots
	auto retval = VCA_real_root_bounds_worker(p1);
	for (auto& root_bound : retval) {
	  root_bound.first /= 2;
	  root_bound.second /= 2;
	}

	//Detect, scale, and shift the second range's roots
	auto second_range = VCA_real_root_bounds_worker(p2);
	for (auto& root_bound : second_range)
	  retval.push_back(std::make_pair(root_bound.first / 2 + 0.5, root_bound.second / 2 + 0.5));
	
	return retval;
      }
    }

    /*! \endcond */

    /*! \brief Calculate bounds on all of the positive real roots of a
        Polynomial.

	This function uses the VCA algorithm to bound the roots and
	assumes the polynomial has non-zero constant and leading order
	coefficients. This function enforces these requirements and
	then simply scales the Polynomial so that all roots lie in the
	range (0,1) before passing it to \ref VCA_real_root_bounds_worker.
     */
    template<size_t Order, class Real>
    containers::StackVector<std::pair<Real,Real>, Order> 
    VCA_real_root_bounds(const Polynomial<Order, Real>& f) {
      //Calculate the upper bound on the polynomial real roots, and
      //return if no roots are detected
      const Real upper_bound = LMQ_upper_bound(f);
      if (upper_bound == 0)
	return containers::StackVector<std::pair<Real,Real>, Order>();
      
      //Scale the polynomial so that all roots lie in the region
      //(0,1), then solve for its roots
      auto bounds = VCA_real_root_bounds_worker(scale_poly(f, upper_bound) );

      //finally, we must undo the scaling on the roots found
      for (auto& bound : bounds) {
	bound.first *= upper_bound;
	bound.second *= upper_bound;
      }
      
      return bounds;
    }

    /*! \brief Class representing the current Mobius transformation
        applied to a Polynomial.
	
	This class is used in the implementation of the
	VAS_real_root_bounds_worker function. A mobius transformation is
	the following function:
	
	\[
	M(x) = \frac{a\,x+b}{c\,x+d}
	\]

	This allows a polynomial to be transformed and (provided the
	Mobius transformation has the same operations applied to it)
	we can map between the original polynomial x and the
	transformed x.
     */
    template<class Real>
    struct MobiusTransform: public std::array<Real, 4> {
      typedef std::array<Real, 4> Base;
      /*! \brief Constructor for the Mobius transformation.
       */
      MobiusTransform(Real a, Real b, Real c, Real d):
	Base({a,b,c,d})
      {}
      
      /*! \brief Evaluate the Mobius transformation at the transformed
          \f$x\f$.
       */
      Real eval(Real x) {
	//Catch the special case that there is no x term
	if (((*this)[0] == 0) && ((*this)[2] == 0))
	  return (*this)[1] / (*this)[3];

	//Special case of inf/inf
	if (std::isinf(x) && ((*this)[0] != 0) && ((*this)[2] != 0))
	  return (*this)[0] / (*this)[2];

	//We have to be careful not to do 0 * inf.
	Real numerator = (*this)[1];
	if ((*this)[0] != 0)
	  numerator += x * (*this)[0];
	
	Real denominator = (*this)[3];
	if ((*this)[2] != 0)
	  denominator += x * (*this)[2];
	
	return numerator / denominator;
      }

      /*!\brief Add the effect of a \ref shift_poly operation
         to the Mobius transform.
       */
      void shift(Real x) {
	(*this)[1] += (*this)[0] * x;
	(*this)[3] += (*this)[2] * x;
      }

      /*!\brief Add the effect of a \ref scale_poly operation to the
         Mobius transform.
       */
      void scale(Real x) {
	(*this)[0] *= x;
	(*this)[2] *= x;
      }
      
      /*!\brief Add the effect of a \ref scale_poly operation to the
         Mobius transform.
       */
      void invert_taylor_shift() {
	Base::operator=(std::array<Real, 4>{(*this)[1], (*this)[0] + (*this)[1], (*this)[3], (*this)[2] + (*this)[3]});
      }
    };

    /*! \brief Calculate bounds on all of the positive real roots in
      the range of a Polynomial.

      This function uses the VAS algorithm to bound the roots and
      assumes the polynomial has non-zero constant and leading order
      coefficients.
      
      The parameter M is an array of Mobius transformation
      coefficients {a,b,c,d}. By default it is a direct mapping to the
      original Polynomial.
     */
    template<size_t Order, class Real>
    containers::StackVector<std::pair<Real,Real>, Order> 
    VAS_real_root_bounds_worker(Polynomial<Order, Real> f, MobiusTransform<Real> M) {
      //This while loop is only used to allow restarting the method
      //without recursion, as this may recurse a large number of
      //times.
      while (true) {
	//Test how many positive roots exist
	const size_t sign_changes = descartes_rule_of_signs(f);
	
	if (sign_changes == 0)
	  //No roots, return empty
	  return containers::StackVector<std::pair<Real,Real>, Order>();
	
	if (sign_changes == 1)
	  //One root! the bounds are M(0) and M(\infty)
	  return containers::StackVector<std::pair<Real,Real>, Order>{std::make_pair(M.eval(0), M.eval(HUGE_VAL))};
	
	auto lb = LMQ_lower_bound(f);
	
	//Attempt to divide the polynomial range from [0, \infty] up
	//into [0, 1] and [1, \infty].
	
	//If there's a large jump in the lower bound, then this will
	//take too long to converge. Try rescaling the polynomial. The
	//factor 16 is empirically selected.
	if (lb >= 16) {
	  f = scale_poly(f, lb);
	  M.scale(lb);
	  lb = 1;
	}

	//Check if the lower bound suggests that this split is a waste
	//of time, if so, shift the polynomial and try again.
	if (lb >= 1) {
	  f = shift_function(f, lb);
	  M.shift(lb);
	  continue; //Start again
	}
      
	//Finally, take a lot of care with roots at the split
	//point. Unfortunately, this test is not as cheap as Descarte's
	//method, but if a root occurs at the split point, numerical
	//error makes it impossible to stably resolve this. Better to
	//just avoid it.
	const Real root_near_1_tol = 1e-8;
	if (alesina_galuzzi_test(f, 1.0 - root_near_1_tol, 1.0 + root_near_1_tol) > 0) {
	  //There is a strong chance of a root in the vicinity of
	  //1. Avoid splitting there by scaling the polynomial to place
	  //these roots at 0.5 and restart
	  const Real scale = 2.0;
	  f = scale_poly(f, scale);
	  M.scale(scale);
	  continue; //Start again
	}

	containers::StackVector<std::pair<Real,Real>, Order> retval;
	//Create the polynomial for [0, 1]
	Polynomial<Order, Real> p01 = invert_taylor_shift(f);
	auto M01 = M;
	M01.invert_taylor_shift();
	//Collect its roots
	auto first_range = VAS_real_root_bounds_worker(p01, M01);
	for (const auto& bound: first_range)
	  retval.push_back(bound);
      
	//Create the polynomial for [1, \infty]
	Polynomial<Order, Real> p1inf = taylor_shift(f);
	auto M1inf = M;
	M1inf.shift(1);
	//Collect its roots
	auto second_range = VAS_real_root_bounds_worker(p1inf, M1inf);
	for (const auto& bound: second_range)
	  retval.push_back(bound);

	return retval;
      }
    }
    
    template<class Real>
    containers::StackVector<std::pair<Real,Real>, 1> 
    VAS_real_root_bounds_worker(const Polynomial<0, Real>& f, MobiusTransform<Real> M) {
      return containers::StackVector<std::pair<Real,Real>, 1>();
    }

    /*! \endcond */
    
    /*! \brief Calculate bounds on all of the positive real roots of a
        Polynomial.
    
    	This function uses the VAS algorithm to bound the roots and
    	assumes the polynomial has non-zero constant and leading order
    	coefficients. This function enforces these conditions before
    	passing it to \ref VAS_real_root_bounds_worker.
     */
    template<size_t Order, class Real>
    containers::StackVector<std::pair<Real,Real>, Order> 
    VAS_real_root_bounds(const Polynomial<Order, Real>& f) {
      //Calculate the upper bound on the polynomial real roots, and
      //return if no roots are detected
      const Real upper_bound = LMQ_upper_bound(f);
      if (upper_bound == 0)
    	return containers::StackVector<std::pair<Real,Real>, Order>();
      
      auto bounds = VAS_real_root_bounds_worker(f, MobiusTransform<Real>(1,0,0,1));
      for (auto& bound: bounds) {
	//Sort the bounds into order
	if (bound.first > bound.second)
	  std::swap(bound.first, bound.second);
	//Replace infinite bounds with the upper bound estimate
	if (std::isinf(bound.second))
	  bound.second = upper_bound;
      }
      return bounds;
    }
    
    enum class  PolyRootBounder {
      VCA, VAS
    };

    /*! \brief Iterative solver for the real roots of a square-free
        Polynomial.

	This function uses an algorithm (VAS/VCA) to bound the roots,
	then a method to calculate them to full precision.
     */
    template<PolyRootBounder BoundMode, size_t Order, class Real>
    containers::StackVector<Real, Order>
    solve_real_roots_poly(const Polynomial<Order, Real>& f) {
      //Handle special cases 

      //The constant coefficient is zero: deflate the polynomial
      if (f[0] == Real())
	return deflate_and_solve_polynomial(f, Real());
      
      //The highest order coefficient is zero: drop to lower order
      //solvers
      if (f[Order] == Real())
	return solve_roots(change_order<Order-1>(f));

      containers::StackVector<std::pair<Real,Real>, Order> bounds;
      containers::StackVector<std::pair<Real,Real>, Order> neg_bounds;
      //Start by solving for the positive roots
      switch (BoundMode){
      case PolyRootBounder::VCA:
	bounds = VCA_real_root_bounds(f);
	neg_bounds = VCA_real_root_bounds(reflect_poly(f));
	break;
      case PolyRootBounder::VAS:
	bounds = VAS_real_root_bounds(f);
	neg_bounds = VAS_real_root_bounds(reflect_poly(f));
	break;
      }
      
      //We need to reflect the x values and put them into the bounds
      //list
      for (auto bound : neg_bounds)
	bounds.push_back(std::make_pair(-bound.second, -bound.first));
      
      //Now bisect to calculate the roots to full precision
      containers::StackVector<Real, Order> retval;
      for (const auto& bound : bounds) {
	boost::uintmax_t iter = 100;
	auto root = boost::math::tools::toms748_solve([&](Real x) { return eval(f, x); }, bound.first, bound.second, boost::math::tools::eps_tolerance<Real>(100), iter);
	retval.push_back((root.first + root.second) / 2);
      }
      return retval;
    }
    
    /*\brief Solve for the distinct real roots of a Polynomial.
      
      This is just a convenience function to select the preferred root
      solver combination.
     */
    template<size_t Order, class Real>
    containers::StackVector<Real, Order>
    solve_roots(const Polynomial<Order, Real>& f) {
      return solve_real_roots_poly<PolyRootBounder::VAS, Order, Real>(f);
    }

    /*! \brief Determine the next positive root of a function.

      This is a generic implementation for Polynomials.
     */
    template<size_t Order, class Real>
    Real next_root(const Polynomial<Order, Real>& f) {
      static_assert(Order, "Not implemented yet!");
    }
    /*! \cond Specializations

      \brief Trivial specialisation for the next positive root of a
      constant.
     */
    template<class Real>
    Real next_root(const Polynomial<0, Real>& f) {
      return std::numeric_limits<Real>::infinity();
    }
    
    /*! \brief Generic implementation of \ref next_root which utilises
      \ref solve_roots.
      
      This is mainly useful for functions with fast root detection
      algorithms. For example, the linear, quadratic, and cubic
      functions are solved via radicals, so this function is used to
      provide an implementation of next_root in these cases.
    */
    template<class Real, size_t Order>
    Real next_root_by_solve_roots(const Polynomial<Order, Real>& f) {
      auto roots = solve_roots(f);
      std::sort(roots.begin(), roots.end());
      for (const Real& root : roots)
	if (root >= 0)
	  return root;
      return std::numeric_limits<Real>::infinity();
    }

    /*! \brief Next positive root for linear functions.
      
      This simply redirects to \ref next_root_by_solve_roots.
    */
    template<class Real>
    Real next_root(const Polynomial<1, Real>& f) {
      return next_root_by_solve_root(f);
    }

    /*! \brief Next positive root for quadratic functions.
      
      This simply redirects to \ref next_root_by_solve_roots.
    */
    template<class Real>
    Real next_root(const Polynomial<2, Real>& f) {
      return next_root_by_solve_root(f);
    }

    /*! \brief Next positive root for cubic functions.
      
      This simply redirects to \ref next_root_by_solve_roots.
    */
    template<class Real>
    Real next_root(const Polynomial<3, Real>& f) {
      return next_root_by_solve_root(f);
    }

    /*! \endcond \} */
    
    /*! \relates Polynomial 
      \name Polynomial bounds
      \{
    */
    /*! \brief The maximum absolute value of an arbitrary order Polynomial in a specified range.
      \param f The Polynomial to evaluate.
      \param tmin The minimum bound.
      \param tmax The maximum bound.
    */
    template<class Real, size_t Order>
    inline auto minmax(const Polynomial<Order, Real>& f, const Real x_min, const Real x_max) -> std::pair<decltype(eval(f, x_min)), decltype(eval(f, x_max))>
    {
      auto roots = solve_roots(derivative(f));
      std::pair<decltype(eval(f, x_min)), decltype(eval(f, x_min))> retval = std::minmax(eval(f, x_min), eval(f, x_max));
      for (auto root : roots)
	if ((root > x_min) && (root < x_max))
	  retval = std::minmax({retval.first, retval.second, eval(f, root)});
      return retval;
    }

    /*! \cond Specializations
      \brief The maximum and minimum values of a 0th order Polynomial in a specified range. 
      \param f The Polynomial to evaluate.
      \param x_min The minimum bound.
      \param x_max The maximum bound.
    */
    template<class Real>
    inline auto minmax(const Polynomial<0, Real>& f, const Real x_min, const Real x_max) -> std::pair<decltype(eval(f, x_min)), decltype(eval(f, x_max))>
    { return std::pair<Real, Real>{eval(f, x_min), eval(f, x_max)}; }

    /*! \brief The maximum and minimum values of a 1st order Polynomial in a specified range. 
      \param f The Polynomial to evaluate.
      \param x_min The minimum bound.
      \param x_max The maximum bound.
    */
    template<class Real>
    inline auto minmax(const Polynomial<1, Real>& f, const Real x_min, const Real x_max) -> std::pair<decltype(eval(f, x_min)), decltype(eval(f, x_max))>
    { return std::pair<Real, Real>{eval(f, x_min), eval(f, x_max)}; }

    /*! \endcond \} */
  }
}
