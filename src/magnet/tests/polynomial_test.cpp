#define BOOST_TEST_MODULE Polynomial_test
#include <boost/test/included/unit_test.hpp>
#include <magnet/math/polynomial.hpp>
#include <magnet/math/vector.hpp>
#include <cmath>
#include <complex>

bool err(double val, double expected)
{
  return std::abs(val / expected - 1) < 0.0001;
}

template<class T1, class T2>
bool compare_expression(const T1& f, const T2& g) {
  std::ostringstream os;
  os << f;
  std::string f_str = os.str();
  os.str(""); os.clear();
  os << g;
  std::string g_str = os.str();
  if (!(f_str == g_str))
    std::cerr << f << " != " << g << std::endl;
  return f_str == g_str;
}

BOOST_AUTO_TEST_CASE( poly_addition )
{
  using namespace magnet::math;
  Polynomial<1> x{0, 2.5};
  Polynomial<0> C{0.3};
  auto poly1 = x+C;
  BOOST_CHECK_EQUAL(poly1[0], 0.3);
  BOOST_CHECK_EQUAL(poly1[1], 2.5);

  auto poly2 = x + 0.3;
  BOOST_CHECK_EQUAL(poly2[0], 0.3);
  BOOST_CHECK_EQUAL(poly2[1], 2.5);
}

BOOST_AUTO_TEST_CASE( poly_multiplication )
{
  using namespace magnet::math;
  Polynomial<1> x{0, 1};
  auto poly1 = -2.0;
  auto poly2 = 2.0 - x + x * x;
  auto poly3 = poly2 * poly1;
  BOOST_CHECK_EQUAL(poly3[0], -4);
  BOOST_CHECK_EQUAL(poly3[1], +2);
  BOOST_CHECK_EQUAL(poly3[2], -2);
}

BOOST_AUTO_TEST_CASE( poly_division )
{
  using namespace magnet::math;
  Polynomial<1> x{0, 1};
  auto poly1 = 2.0 - x + x * x;
  auto poly2 = poly1 / 0.5;
  BOOST_CHECK_EQUAL(poly2[0], 4);
  BOOST_CHECK_EQUAL(poly2[1], -2);
  BOOST_CHECK_EQUAL(poly2[2], 2);
}

BOOST_AUTO_TEST_CASE( poly_vector )
{
  using namespace magnet::math;
  Polynomial<1, Vector> x{Vector(), Vector{1,2,3}};
  Polynomial<0, Vector> C{Vector{3,2,1}};
  auto poly1 = x+C;
  BOOST_CHECK_EQUAL(poly1[0], (Vector{3,2,1}));
  BOOST_CHECK_EQUAL(poly1[1], (Vector{1,2,3}));
  
  auto poly2 = poly1 * poly1;
  BOOST_CHECK_EQUAL(poly2[0], 14);
  BOOST_CHECK_EQUAL(poly2[1], 20);
  BOOST_CHECK_EQUAL(poly2[2], 14);
}

BOOST_AUTO_TEST_CASE( poly_lower_order )
{
  using namespace magnet::math;
  Polynomial<1> x{0, 1};
  Polynomial<2> poly2 = 2.0 - x + x * x;
  Polynomial<3> poly3 = poly2 + 0 * x * x * x;
  //Try to cast down one level as the highest order coefficient is zero
  BOOST_CHECK(poly3[3] == 0);
  Polynomial<2> poly4(poly3);

  BOOST_CHECK_EQUAL(poly4[0], 2);
  BOOST_CHECK_EQUAL(poly4[1], -1);
  BOOST_CHECK_EQUAL(poly4[2], 1);
  BOOST_CHECK_EQUAL(eval(poly3, 123), eval(poly4, 123));
}

BOOST_AUTO_TEST_CASE( poly_derivative )
{
  using namespace magnet::math;
  Polynomial<1> x{0, 1};
  auto poly1 = x + x*x + x*x*x + x*x*x*x;
  auto poly2 = derivative(poly1);
  BOOST_CHECK_EQUAL(poly2[0], 1);
  BOOST_CHECK_EQUAL(poly2[1], 2);
  BOOST_CHECK_EQUAL(poly2[2], 3);  
  BOOST_CHECK_EQUAL(poly2[3], 4);  

  auto poly3 = 2.0 - x + 2 * x * x;
  auto poly4 = derivative(poly3);
  BOOST_CHECK_EQUAL(poly4[0], -1);
  BOOST_CHECK_EQUAL(poly4[1], 4);
  BOOST_CHECK_EQUAL(eval(poly4, 0), -1);
  BOOST_CHECK_EQUAL(eval(poly4, 1), 3);
}

BOOST_AUTO_TEST_CASE( poly_zero_derivative)
{
  using namespace magnet::math;
  const Polynomial<1> x{0, 1};
  const auto poly1 = derivative(x);
  BOOST_CHECK_EQUAL(poly1[0], 1);

  const auto poly2 = derivative(poly1);
  BOOST_CHECK(compare_expression(poly2, NullSymbol()));
}

BOOST_AUTO_TEST_CASE( poly_quadratic_roots)
{
  using namespace magnet::math;
  Polynomial<1> x{0, 1};
  
  {//Quadratic with no roots
    auto poly = x * x - 3 * x + 4;
    auto roots = solve_roots(poly);
    BOOST_CHECK(roots.size() == 0);
  }
  
  {//Quadratic with one root
    auto poly = -4 * x * x + 12 * x - 9;
    auto roots = solve_roots(poly);
    BOOST_CHECK(roots.size() == 1);
    if (roots.size() == 1)
      BOOST_CHECK_CLOSE(roots[0], 1.5, 1e-10);
  }
  
  {//linear function with one root
    auto poly =  0 * x * x + 12 * x - 9;
    auto roots = solve_roots(poly);
    BOOST_CHECK(roots.size() == 1);
    if (roots.size() == 1)
      BOOST_CHECK_CLOSE(roots[0], 0.75, 1e-10);
  }

  {//constant function, with no roots
    auto poly =  0 * x * x + 0 * x - 9;
    auto roots = solve_roots(poly);
    BOOST_CHECK(roots.size() == 0);
  }
}

BOOST_AUTO_TEST_CASE( poly_quadratic_special_cases)
{
  using namespace magnet::math;
  Polynomial<1> x{0, 1};
  
  {//Quadratic with catastrophic cancellation of error
    auto poly = x * x + 712345.12 * x + 1.25;  
    auto roots = solve_roots(poly);
    BOOST_CHECK(roots.size() == 2);

    std::sort(roots.begin(), roots.end());
    BOOST_CHECK_CLOSE(roots[0], -712345.1199985961, 1e-10);
    BOOST_CHECK_CLOSE(roots[1], -1.754767408250742e-6, 1e-10);
  }

  const double maxsqrt = std::sqrt(std::numeric_limits<double>::max());
  const double largeterm = maxsqrt * 100;
  {//Large linear coefficient
    auto poly = x * x + largeterm * x + 1.25;
    auto roots = solve_roots(poly);
    BOOST_CHECK(roots.size() == 2);

    std::sort(roots.begin(), roots.end());
    //Mathematica value
    BOOST_CHECK_CLOSE(roots[0], -1.3407807929942599e156, 1e-10);
    BOOST_CHECK_CLOSE(roots[1], -9.322925914000258e-157, 1e-10);
  }

  {//Large (+ve) constant coefficient
    auto poly = x * x + x + largeterm;
    auto roots = solve_roots(poly);
    BOOST_CHECK(roots.size() == 0);
  }
  {//Large (-ve) constant coefficient
    auto poly = x * x + x - largeterm;
    auto roots = solve_roots(poly);
    BOOST_CHECK(roots.size() == 2);

    std::sort(roots.begin(), roots.end());
    //Mathematica value
    BOOST_CHECK_CLOSE(roots[0], -1.157920892373162e78, 1e-10);
    BOOST_CHECK_CLOSE(roots[1], 1.157920892373162e78, 1e-10);
  }
}

double cubic_rootvals[] = {-1e6, -1e3, -100, -1, 0, 1, +100, 1e3, 1e6};

BOOST_AUTO_TEST_CASE( poly_cubic_triple_roots )
{
  using namespace magnet::math;
  const Polynomial<1> x{0, 1};

  for (double root1 : cubic_rootvals)
    for (double root2 : cubic_rootvals)
      if (root2 != root1)
	for (double root3 : cubic_rootvals)
	  if ((root1 != root2) && (root2 != root3) && (root1 != root3))
	    {
	      auto f = (x - root1) * (x - root2) * (x - root3);
	      //Don't test the case where there is only one root (x^3=c)
	      if ((f[2] == 0) && (f[1] == 0)) continue;

	      auto roots = solve_roots(f);
	      decltype(roots) actual_roots = {root1, root2, root3};
	      std::sort(actual_roots.begin(), actual_roots.end());
	      std::sort(roots.begin(), roots.end());

	      BOOST_CHECK_MESSAGE(roots.size() == 3, f << " roots=[" << roots[0] << "," << roots[1] << "," << roots[2] << "] actual_roots=[" << root1 << "," << root2 << "," << root3 << "]");

	      if (roots.size() == 3)
		for (size_t i = 0; i < 3; ++i)
		  {
		    double root_error = std::abs((roots[i] - actual_roots[i]) / (actual_roots[i] + (actual_roots[i] == 0)));
		    BOOST_CHECK_MESSAGE(root_error < 0.001, "root_error=" << root_error << " " << f << " roots=[" << roots[0] << "," << roots[1] << "," << roots[2] << "] actual_roots=[" << root1 << "," << root2 << "," << root3 << "]");
		  }
	    }
}

BOOST_AUTO_TEST_CASE( poly_cubic_single_roots )
{
  using namespace magnet::math;
  const Polynomial<1> x{0, 1};

  for (double root1 : cubic_rootvals)
    for (double root2real : cubic_rootvals)
      for (double root2im : cubic_rootvals)
	{
	  //Skip real second root cases, and the symmetric cases
	  if (root2im <= 0) continue;
	  
	  std::complex<double>
	    root1val(root1, 0),
	    root2val(root2real, root2im),
	    root3val(root2real, -root2im);
	  
	  auto poly_c = (x - root1val) * (x - root2val) * (x - root3val);
	  
	  Polynomial<3, double> f = poly_c[0].real() + poly_c[1].real() * x + poly_c[2].real() * x * x + poly_c[3].real() * x * x * x;

	  auto roots = solve_roots(f);
	  BOOST_CHECK_MESSAGE(roots.size() == 1, "rootcount=" << roots.size() << " " << f << " roots=[" << roots[0] << "," << roots[1] << "," << roots[2] << "] actual_roots=[" << root1 << "," << root2real << " +- " << root2im << "i]");

	  if (roots.size() == 1)
	    {
	      double root_error = std::abs((roots[0] - root1) / (root1 + (root1 == 0)));
	      BOOST_CHECK_MESSAGE(root_error < 0.001, "root error is " << root_error);
	    }
	}
}

BOOST_AUTO_TEST_CASE( poly_cubic_special_cases )
{
  using namespace magnet::math;
  const Polynomial<1> x{0, 1};

  {//Zero constant term with three roots
    auto poly = (x * x + 712345.12 * x + 1.25) * x;
    auto roots = solve_roots(poly);
    BOOST_CHECK(roots.size() == 3);

    std::sort(roots.begin(), roots.end());    
    BOOST_CHECK_CLOSE(roots[0], -712345.1199985961, 1e-10);
    BOOST_CHECK_CLOSE(roots[1], -1.754767408250742e-6, 1e-10);
    BOOST_CHECK_CLOSE(roots[2], 0, 1e-10);
  }

  {//Zero constant term with one root
    auto poly = (x * x - 3 * x + 4) * x;
    auto roots = solve_roots(poly);
    BOOST_CHECK(roots.size() == 1);
    BOOST_CHECK_CLOSE(roots[0], 0 , 1e-10);
  }

  {//Special case where f(x) = a * x^3 + d
    auto poly = x * x * x + 1e3;
    auto roots = solve_roots(poly);
    BOOST_CHECK(roots.size() == 1);
    BOOST_CHECK_CLOSE(roots[0], -10, 1e-10);

    poly = x * x * x - 1e3;
    roots = solve_roots(poly);
    BOOST_CHECK(roots.size() == 1);
    BOOST_CHECK_CLOSE(roots[0], 10, 1e-10);
  }

  const double maxsqrt = std::sqrt(std::numeric_limits<double>::max());
  const double largeterm = maxsqrt * 100;

  {//Large x^2 term
    auto poly = x * x * x - largeterm * x * x + 1.25;
    auto roots = solve_roots(poly);
    std::sort(roots.begin(), roots.end());
    BOOST_CHECK_EQUAL(roots.size(), 3);
    BOOST_CHECK_CLOSE(roots[0], -9.655529977168658e-79, 1e-10);
    BOOST_CHECK_CLOSE(roots[1], +9.655529977168658e-79, 1e-10);
    BOOST_CHECK_CLOSE(roots[2], 1.3407807929942599e156, 1e-10);
  }

  {//Large x term
    auto poly = x * x * x - x * x - largeterm * x + 1.25;
    auto roots = solve_roots(poly);
    std::sort(roots.begin(), roots.end());
    BOOST_CHECK_EQUAL(roots.size(), 3);
    BOOST_CHECK_CLOSE(roots[0], -1.1579208923731622e78, 1e-10);
    BOOST_CHECK_CLOSE(roots[1], 9.322925914000258e-157, 1e-10);
    BOOST_CHECK_CLOSE(roots[2], 1.1579208923731622e78, 1e-10);
  }

  const double smallerterm = maxsqrt * 1e-1;
  {
    //Large v term
    auto poly = x * x * x  - smallerterm * x * x - smallerterm * x + 2;
    auto roots = solve_roots(poly);
    std::sort(roots.begin(), roots.end());
    BOOST_CHECK_EQUAL(roots.size(), 3);
    BOOST_CHECK_CLOSE(roots[0], -1.0, 1e-10);
    BOOST_CHECK_CLOSE(roots[1], 1.491668146240041472864517142264024641421371730393e-153, 1e-10);
    BOOST_CHECK_CLOSE(roots[2], 1.340780792994259598314974448015366224371799690462e153, 1e-10);
  }
}

BOOST_AUTO_TEST_CASE( poly_bounds)
{
  using namespace magnet::math;
  const Polynomial<1> x{0, 1};

  { //Check constant
    auto f1 = Polynomial<0>{23};
    auto bounds = minmax(f1, -4.0, 10.0);
    BOOST_CHECK_CLOSE(bounds.first, 23, 1e-10);
    BOOST_CHECK_CLOSE(bounds.second, 23, 1e-10);
  }

  { //Check linear
    auto f1 = 2*x +12;
    auto bounds = minmax(f1, -4.0, 10.0);
    BOOST_CHECK_CLOSE(bounds.first, 4, 1e-10);
    BOOST_CHECK_CLOSE(bounds.second, 32, 1e-10);
  }

  { //Check quadratic
    auto f1 = x * x + 2*x +12;
    auto bounds = minmax(f1, -4.0, 10.0);
    BOOST_CHECK_CLOSE(bounds.first, 11, 1e-10);
    BOOST_CHECK_CLOSE(bounds.second, 132, 1e-10);
  }

  {//Check cubic
    auto f1 = 4 * (x*x*x) - x * x - 2*x +12;

    auto roots = solve_roots(f1);
    BOOST_CHECK(roots.size() == 1);
    std::sort(roots.begin(), roots.end());
    BOOST_CHECK_CLOSE(roots[0], -1.472711896724616002268033950475380144341, 1e-10);

    auto droots = solve_roots(derivative(f1));
    BOOST_CHECK(droots.size() == 2);
    std::sort(droots.begin(), droots.end());
    BOOST_CHECK_CLOSE(droots[0], -1.0/3, 1e-10);
    BOOST_CHECK_CLOSE(droots[1], 0.5, 1e-10);

    auto bounds = minmax(f1, -4.0, 10.0);    
    BOOST_CHECK_CLOSE(bounds.first, -252, 1e-10);
    BOOST_CHECK_CLOSE(bounds.second, 3892, 1e-10);
  }

  {//Check quartic
    auto f1 = 10 * (x*x*x*x) + x*x*x - 30 * x * x -23;

    //Quartic roots aren't available yet
    //auto roots = solve_roots(f1);
    //BOOST_CHECK(roots.size() == 1);
    //std::sort(roots.begin(), roots.end());
    //BOOST_CHECK_CLOSE(roots[0], -1.949403904489790210996459054473124835057, 1e-10);
    //BOOST_CHECK_CLOSE(roots[1], +1.864235880634589025006445510389799368569, 1e-10);

    auto droots = solve_roots(derivative(f1));
    BOOST_CHECK_EQUAL(droots.size(),3);
    std::sort(droots.begin(), droots.end());
    BOOST_CHECK_CLOSE(droots[0], -1.262818836058599076329128653113014315066, 1e-10);
    BOOST_CHECK_CLOSE(droots[1], 0, 1e-10);
    BOOST_CHECK_CLOSE(droots[2], +1.187818836058599076329128653113014315066, 1e-10);

    auto bounds = minmax(f1, -4.0, 10.0);
    BOOST_CHECK_CLOSE(bounds.first, -47.42412909307610601944478081683796164898, 1e-10);
    BOOST_CHECK_CLOSE(bounds.second, 97977.0, 1e-10);
  }

  {//Check PowerOp quartic
    auto f1 = pow<2>(30 * x * x + x - 23);
    
    //Quartic roots aren't available yet
    //auto roots = solve_roots(f1);
    //BOOST_CHECK(roots.size() == 2);
    //std::sort(roots.begin(), roots.end());
    //NOTE THESE ROOTS ARE DOUBLE ROOTS (roots.size() may equal 2,3, or 4)
    //BOOST_CHECK_CLOSE(roots[0], -0.8924203103613100773375343963347855860436, 1e-10);
    //BOOST_CHECK_CLOSE(roots[1], 0.8590869770279767440042010630014522527103, 1e-10);

    auto droots = solve_roots(derivative(f1));
    BOOST_CHECK(droots.size() == 3);
    std::sort(droots.begin(), droots.end());
    BOOST_CHECK_CLOSE(droots[0], -0.8924203103613100773375343963347855860436, 1e-10);
    BOOST_CHECK_CLOSE(droots[1], -0.01666666666666666666666666666666666666667, 1e-10);
    BOOST_CHECK_CLOSE(droots[2], +0.8590869770279767440042010630014522527103, 1e-10);

    auto bounds = minmax(f1, -4.0, 10.0);
    BOOST_CHECK_CLOSE(bounds.first, 0, 1e-10);
    BOOST_CHECK_CLOSE(bounds.second, 8.922169e6, 1e-10);
  }
}
