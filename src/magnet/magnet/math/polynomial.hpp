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

#include <magnet/exception.hpp>
#include <magnet/containers/stack_vector.hpp>
#include <magnet/math/precision.hpp>
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

    /*! \brief Representation of Polynomial with basic algebra operations.

      This class allows basic computer algebra to be performed with
      Polynomial equations.
      
      For example, the polynomial $f(x)=x^2 + 2\,x + 3$ can be created like so:
      \code{.cpp}
      Polynomial<1> x{0,1};
      auto f = x*x + 2*x +3;    
      \endcode
      And evaluated at the point $x=3$ like so:
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

      /*! \brief Constructor for lowering the order of a Polynomial.
	
	This constructor should only be used if the highest
	coefficient in the passed Polynomial is zero. This requirement
	is only enforced in debug mode.
      */
      Polynomial(const Polynomial<Order+1, Real>& poly) {
#ifdef MAGNET_DEBUG
	if (poly[Order+1] != 0)
	  M_throw() << "Trying to reduce the order of a polynomial with non-zero highest order coefficients!";
#endif
	std::copy(poly.begin(), poly.end()-1, Base::begin());
      };
      
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
	auto it = _list.begin();
	for (; it != _list.end(); ++i, ++it)
	  Base::operator[](i) = *it;

	for (; i < Order+1; ++i)
	  Base::operator[](i) = Real();
      }

      /*! \brief Constructor for promoting from low-order to higher-order Polynomial types. */
      template<size_t N, class Real2>
	Polynomial(const Polynomial<N, Real2>& poly) {
	static_assert(N <= Order, "Can only promote to higher order polynomials");
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

      /*! \brief Evaluate the polynomial at x. */
      Real operator()(Real x) const {
	Real sum = Base::operator[](Order);
	for(int i = Order - 1; i >= 0; --i)
	  {
	    sum *= x;
	    sum += Base::operator[](i);
	  }
	return sum;
      }
    };

    /*! \relates Polynomial 
      \name Polynomial algebraic operations
      
       For all operations below we do not assume that we have a
       closure. For example, a vector multiplied by a vector is a
       scalar therefore the * operator may change the returned type of
       the polynomial.
       \{
    */
    template<class Real1, class Real2, size_t N>
    auto operator+(const Real1& r, const Polynomial<N,Real2>& poly)->Polynomial<N, decltype(poly[0] + r)>
    { return poly + r; }

    /*!\brief Addition operator for Polynomial and constants. */
    template<class Real1, class Real2, size_t N>
    auto operator+(const Polynomial<N,Real1>& poly, const Real2& r)->Polynomial<N, decltype(poly[0] + r)>
    {
      Polynomial<N, decltype(poly[0] + r)> retval(poly);
      retval[0] += r;
      return retval;
    }
    /*!\brief Addition operator for two Polynomial types. */
    template<size_t M, size_t N, class Real1, class Real2>
    auto operator+(const Polynomial<M, Real1>& poly1, const Polynomial<N, Real2>& poly2)->Polynomial<detail::max_order(M, N), decltype(poly1[0] + poly2[0])>
    {
      Polynomial<detail::max_order(M, N), decltype(poly1[0] + poly2[0])> retval(poly1);
      for (size_t i(0); i <= N; ++i)
	retval[i] += poly2[i];
      return retval;
    }

    /*! \brief Subtraction of a Polynomial from a constant. */
    template<class Real1, class Real2, size_t N>
    auto operator-(const Real1& r, const Polynomial<N, Real2>& poly)->Polynomial<N,decltype((-poly)[0]+r)>
    {
      Polynomial<N, decltype((-poly)[0]+r)> retval = -poly;
      retval[0] += r;
      return retval;  
    }
    /*! \brief Subtraction of a constant from a Polynomial. */
    template<class Real1, class Real2, size_t N>
    auto operator-(const Polynomial<N,Real1>& poly, const Real2& r)->Polynomial<N,decltype(poly[0]-r)>
    {
      Polynomial<N,decltype(poly[0]-r)> retval(poly);
      retval[0] -= r;
      return retval;
    }
    /*! \brief Subtraction between two Polynomial types. */
    template<class Real1, class Real2, size_t M, size_t N>
    auto operator-(const Polynomial<M,Real1>& poly1, const Polynomial<N,Real2>& poly2)->Polynomial<detail::max_order(M, N),decltype(poly1[0]-poly2[0])>
    {
      Polynomial<detail::max_order(M, N),decltype(poly1[0]-poly2[0])> retval(poly1);
      for (size_t i(0); i <= N; ++i)
	retval[i] -= poly2[i];
      return retval;
    }

    template<class Real1, class Real2, size_t N>
    auto operator*(const Real1& r, const Polynomial<N, Real2>& poly) -> Polynomial<N, decltype(poly[0] * r)>
    { return poly * r; }

    /*! \brief Multiplication between a constant and a Polynomial.
      
      We assume that multiplication commutes. 
    */
    template<class Real1, class Real2, size_t N>
    auto operator*(const Polynomial<N, Real1>& poly, const Real2& r) -> Polynomial<N, decltype(poly[0] * r)>
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
    auto operator/(const Polynomial<N, Real1>& poly, const Real2& r) -> Polynomial<N, decltype(Real1() / Real2())>
    {
      Polynomial<N, decltype(Real1() / Real2())> retval;
      for (size_t i(0); i <= N; ++i)
	retval[i] = poly[i] / r;
      return retval;
    }

    /*! \} */

    /*! \relates Polynomial 
      \name Polynomial calculus operations
      \{
    */

    /*! \relates Polynomial 
      \brief Derivatives of Polynomial classes (constants).*/
    template<class Real, size_t N>
    inline Polynomial<N-1, Real> derivative(const Polynomial<N, Real>& f) {
      Polynomial<N-1, Real> retval;
      for (size_t i(0); i < N; ++i) {
	retval[i] = f[i+1] * (i+1);
      }
      return retval;
    }
    
    /*! \relates Polynomial 
       \brief Specialisation for derivatives of 0th order Polynomial classes (constants).
    */
    template<class Real>
    inline Polynomial<0, Real> derivative(const Polynomial<0, Real>& f) {
      return Polynomial<0, Real>{0};
    }

    /*! \} */

    /*! \relates Polynomial 
      \name Polynomial input/output operations
      \{
    */
    /*! \brief Writes a human-readable representation of the Polynomial to the output stream. */
    template<class Real, size_t N>
    inline std::ostream& operator<<(std::ostream& os, const Polynomial<N, Real>& poly) {
      for (size_t i(N); i != 0; --i) {
	if (poly[i] == Real()) continue;
	os << "(" << poly[i] << ") * x";
	if (i > 1) 
	  os << "^" << i;
	os << " + ";
      }
      os << poly[0];
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

      This formula can be used to calculate all coefficients using the
      starting point \f$b_0=-a_0 / r_1\f$. This approach is not stable
      if the root is zero, or if \f$b_{i-1}\f$ has the same sign as
      \f$a_i\f$.
      
      An alternative iterative form may be found by substituting
      \f$i\to i+1\f$, which gives:

      \f[
      b_{i} = a_{i+1} + r_1\,b_{i+1} \qquad \textrm{for}\ i\in[0,\,N-2]
      \f]

      Again this approach may be used given the starting point
      \f$b_{N-1}=a_N\f$. However, it is not stable if \f$a_{i+1}\f$
      has the opposite sign to \f$r_1\,b_{i+1}\f$.
      
      As both approaches may suffer from catastrophic cancellation, we
      should switch between them. We prefer the second (forwards)
      approach as there is no divide and error does not accumulate as
      much for the higher order coefficients. This means we only
      switch to the backwards approach if the forwards approach has
      cancellation AND the backwards approach does not. An approach is
      dangerous if two non-zero terms are subtracted from each other
      (i.e., for the first approach this happens if \f$a_{i+1}\f$ and
      \f$r_1\,b_{i+1}\f$ are non-zero and have opposite sign).

      It should be noted that Numerical Recipies states that "Forward
      deflation is stable if the largest absolute root is always
      divided out... backward deflation is stable if the smallest
      absolute root is always divided out". Unfortunately we do not
      know the magnitude of the root being divided out.

      \param f The Polynomial to factor a root out of.
      \param root The root to remove.
     */
    template<size_t Order, class Real>
    inline Polynomial<Order-1, Real> deflate_polynomial(const Polynomial<Order, Real>& a, const double root) {
      Polynomial<Order-1, Real> b;
      
      //Check for the simple case where root==0
      if (root == 0) {
	for (size_t i(0); i < Order-1; ++i)
	  b[i] = a[i];
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

    /*! \brief A dummy function which returns no roots of a 0th order Polynomial.
      \param f The Polynomial to evaluate.
     */
    inline containers::StackVector<double, 0> solve_roots(const Polynomial<0, double>& f) {
      return containers::StackVector<double, 0>();
    }

    /*! \brief The root of a 1st order Polynomial.
      \param f The Polynomial to evaluate.
     */
    inline containers::StackVector<double, 1> solve_roots(const Polynomial<1, double>& f) {
      containers::StackVector<double, 1> roots;
      if (f[1] != 0)
	roots.push_back(-f[0] / f[1]);
      return roots;
    }

    /*! \brief The roots of a 2nd order Polynomial.
      \param f The Polynomial to evaluate.
     */
    inline containers::StackVector<double, 2> solve_roots(Polynomial<2, double> f) {
      //If this is actually a linear polynomial, drop down to that solver.
      if (f[2] == 0) 
	return solve_roots(Polynomial<1, double>(f));
      
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
      /*! \brief Uses a quadratic scheme to polish up a root.
      */
      inline void cubicNewtonRootPolish(const Polynomial<3, double>& f, double& root)
      {
	//Stored for comparison later
	double error = f(root);
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
	    error = f(root);
	  }
      }
    }


    /*! \brief The roots of a 3rd order Polynomial.
      \param f The Polynomial to evaluate.
     */
    inline containers::StackVector<double, 3> solve_roots(const Polynomial<3, double>& f_original) {
      //Ensure this is actually a third order polynomial
      if (f_original[3] == 0)
	return solve_roots(Polynomial<2, double>(f_original));
      
      if (f_original[0] == 0)
	{
	  //If the constant is zero, one root is x=0.  We can divide
	  //by x and solve the remaining quadratic
	  containers::StackVector<double, 3> roots = solve_roots(Polynomial<2, double>{f_original[1], f_original[2], f_original[3]});
	  roots.push_back(0);
	  return roots;
	}

      //Convert to a cubic with a unity high-order coefficient
      auto f = f_original / f_original[3];
      
      if ((f[2] == 0) && (f[1] == 0))
	//Special case where f(x) = x^3 + f[0]
	return containers::StackVector<double, 3>{std::cbrt(-f[0])};

      static const double maxSqrt = std::sqrt(std::numeric_limits<double>::max());
      
      if ((f[2] > maxSqrt) || (f[2] < -maxSqrt))
	{
	  //The equation is limiting to x^3 + f[2] * x^2 == 0. Use
	  //this to estimate the location of one root, polish it up,
	  //then deflate the polynomial and solve the quadratic.
	  double largeroot = -f[2];
	  cubicNewtonRootPolish(f, largeroot);
	  containers::StackVector<double, 3> roots = solve_roots(deflate_polynomial(f, largeroot));
	  roots.push_back(largeroot);
	  return roots;
	}

//NOT SURE THESE RANGE TESTS ARE BENEFICIAL
//      if (f[1] > maxSqrt)
//	//Special case, if f[1] is large (and f[2] is not) and the root
//	//is -f[0] / f[1], the x^3 term is negligble, and all other
//	//terms cancel.
//	return containers::StackVector<double, 3>{-f[0] / f[1]};
//
//      if (f[1] < -maxSqrt)
//	//Special case, equation is x^3 + q x == 0
//	return containers::StackVector<double, 3>{-std::sqrt(-f[1])};
//
//      if ((f[0] > maxSqrt) || (f[0] < -maxSqrt))
//	//Another special case where f(x)= x^3 +f[0]
//	return containers::StackVector<double, 3>{-std::cbrt(f[0])};

      const double v = f[0] + (2.0 * f[2] * f[2] / 9.0 - f[1]) * (f[2] / 3.0);

//      if ((v > maxSqrt) || (v < -maxSqrt))
//	return containers::StackVector<double, 3>{-f[2]};
      
      const double uo3 = f[1] / 3.0 - f[2] * f[2] / 9.0;
      const double u2o3 = uo3 + uo3;
      
//      if ((u2o3 > maxSqrt) || (u2o3 < -maxSqrt))
//	{
//	  if (f[2]==0)
//	    {
//	      if (f[1] > 0)
//		return containers::StackVector<double, 3>{-f[0] / f[1]};
//	      
//	      if (f[1] < 0)
//		return containers::StackVector<double, 3>{-std::sqrt(-f[1])};
//		
//	      return containers::StackVector<double, 3>{0};
//	    }
//
//	  return containers::StackVector<double, 3>{-f[1] / f[2]};
//	}

      const double uo3sq4 = u2o3 * u2o3;
//      if (uo3sq4 > maxSqrt)
//	{
//	  if (f[2] == 0)
//	    {
//	      if (f[1] > 0)
//		return containers::StackVector<double, 3>{-f[0] / f[1]};
//
//	      if (f[1] < 0)
//		return containers::StackVector<double, 3>{-std::sqrt(-f[1])};
//
//	      return containers::StackVector<double, 3>{0};
//	    }
//
//	  return containers::StackVector<double, 3>{-f[1] / f[2]};
//	}

      const double j = (uo3sq4 * uo3) + v * v;
  
      if (j > 0) 
	{//Only one root (but this test can be wrong due to a
	  //catastrophic cancellation in j 
	  //(i.e. (uo3sq4 * uo3) == v * v)
	  
	  const double w = std::sqrt(j);
	  double root1;
	  if (v < 0)
	    root1 = std::cbrt(0.5*(w-v)) - (uo3) * std::cbrt(2.0 / (w-v)) - f[2] / 3.0;
	  else
	    root1 = uo3 * std::cbrt(2.0 / (w+v)) - std::cbrt(0.5*(w+v)) - f[2] / 3.0;

	  //We now polish the root up before we use it in other calculations
	  cubicNewtonRootPolish(f, root1);
	 
	  //We double check that there are no more roots by deflating
	  //the polynomial using the calculated root and solving this.
	  containers::StackVector<double, 3> roots = solve_roots(deflate_polynomial(f, root1));
	  roots.push_back(root1);
	  return roots;
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
    
    /*! \} */

    /*! \relates Polynomial 
      \name Polynomial bounds
      \{
    */
    /*! \brief The maximum absolute value of a 0th order Polynomial in a specified range. 
      \param f The Polynomial to evaluate.
      \param tmin The minimum bound.
      \param tmax The maximum bound.
     */
    template<class Real>
    inline Real max_abs_val(const Polynomial<0, Real>& f, const double tmin, const double tmax) {
      return std::abs(f[0]);
    }

    /*! \brief The maximum absolute value of a 1st order Polynomial in a specified range. 
      \param f The Polynomial to evaluate.
      \param tmin The minimum bound.
      \param tmax The maximum bound.
     */
    template<class Real>
    inline Real max_abs_val(const Polynomial<1, Real>& f, const double tmin, const double tmax) {
      return std::max(std::abs(f(tmin)), std::abs(f(tmax)));
    }

    /*! \brief The maximum absolute value of an arbitrary order Polynomial in a specified range.
      \param f The Polynomial to evaluate.
      \param tmin The minimum bound.
      \param tmax The maximum bound.
     */
    template<class Real, size_t Order>
    inline Real max_abs_val(const Polynomial<Order, Real>& f, const double tmin, const double tmax) {
      auto roots = solve_roots(derivative(f));
      Real retval = std::max(std::abs(f(tmin)), std::abs(f(tmax)));
      for (auto root : roots)
	if ((root > tmin) && (root < tmax))
	  retval = std::max(std::abs(f(root)), retval);
      return retval;
    }
    /*! \} */
  }
}