#define BOOST_TEST_MODULE Polynomial_test
#include <boost/test/included/unit_test.hpp>
#include <magnet/math/polynomial.hpp>
#include <cmath>

bool err(double val, double expected)
{
  return std::abs(val / expected - 1) < 0.0001;
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

BOOST_AUTO_TEST_CASE( poly_derivative )
{
  using namespace magnet::math;
  Polynomial<1> x{0, 1};
  auto poly2 = 2.0 - x + x * x;
  auto poly3 = derivative(poly2);
  
  BOOST_CHECK_EQUAL(poly3[0], -1);
  BOOST_CHECK_EQUAL(poly3[1], 1);
  BOOST_CHECK_EQUAL(poly3(0), -1);
  BOOST_CHECK_EQUAL(poly3(1), 0);
}
