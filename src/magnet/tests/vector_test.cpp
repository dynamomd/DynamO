#define BOOST_TEST_MODULE Vector_test
#include <boost/test/unit_test.hpp>
#include <magnet/math/vector.hpp>
#include <cmath>

bool err(double val, double expected)
{
  return std::abs(val / expected - 1) < 0.0001;
}

BOOST_AUTO_TEST_CASE( vector_addition )
{
  magnet::math::Vector A(1,2,3);
  magnet::math::Vector B(4,5,6);

  magnet::math::Vector C = A+B;
  BOOST_CHECK_EQUAL(C[0], 5);
  BOOST_CHECK_EQUAL(C[1], 7);
  BOOST_CHECK_EQUAL(C[2], 9);
}

BOOST_AUTO_TEST_CASE( vector_norm )
{
  magnet::math::Vector B(1,1,1);

  BOOST_CHECK(err(B.nrm(), std::sqrt(3)));
}
