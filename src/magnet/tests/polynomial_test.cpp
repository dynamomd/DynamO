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
  BOOST_CHECK_EQUAL(poly3(123), poly4(123));
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
  BOOST_CHECK_EQUAL(poly4(0), -1);
  BOOST_CHECK_EQUAL(poly4(1), 3);
}

BOOST_AUTO_TEST_CASE( poly_zero_derivative)
{
  using namespace magnet::math;
  const Polynomial<1> x{0, 1};
  const auto poly1 = derivative(x);
  BOOST_CHECK_EQUAL(poly1[0], 1);

  const auto poly2 = derivative(poly1);
  BOOST_CHECK_EQUAL(poly2[0], 0);

  const auto poly3 = derivative(poly2);
  BOOST_CHECK_EQUAL(poly3[0], 0);

  const auto poly4 = derivative(poly3);
  BOOST_CHECK_EQUAL(poly4[0], 0);
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
    
    BOOST_CHECK_CLOSE(roots[0], 1.5, 1e-10);
  }
  
  {//linear function with one root
    auto poly =  0 * x * x + 12 * x - 9;
    auto roots = solve_roots(poly);
    BOOST_CHECK(roots.size() == 1);
    
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

double cubic_rootvals[] = {-1e7, -1e6, -1e3, -100, -1, 0, 1, +100, 1e3, 1e6, 1e7};

BOOST_AUTO_TEST_CASE( poly_cubic_triple_roots )
{
  using namespace magnet::math;
  const Polynomial<1> x{0, 1};

  for (double root1 : cubic_rootvals)
    for (double root2 : cubic_rootvals)
      if (root2 != root1)
	for (double root3 : cubic_rootvals)
	  if ((root3 != root2) && (root3 != root1))
	    {
	      auto f = (x - root1) * (x - root2) * (x - root3);
	      //Don't test the case where there is only one root (x^3=c)
	      if ((f[2] == 0) && (f[1] == 0)) continue;

	      auto roots = solve_roots(f);
	      decltype(roots) actual_roots = {root1, root2, root3};
	      std::sort(actual_roots.begin(), actual_roots.end());
	      std::sort(roots.begin(), roots.end());

	      BOOST_CHECK_EQUAL(roots.size(), 3);

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
	  
	  double a = (-root1val - root2val  - root3val).real(),
	    b = (root1val * root2val + root1val * root3val + root2val * root3val).real(),
	    c = - (root1val * root2val * root3val).real();
	  
	  auto f = ((x + a) * x + b) * x + c;

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
}
