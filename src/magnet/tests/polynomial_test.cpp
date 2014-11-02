#define BOOST_TEST_MODULE Polynomial_test
#include <boost/test/included/unit_test.hpp>
#include <magnet/math/polynomial.hpp>
#include <cmath>

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
