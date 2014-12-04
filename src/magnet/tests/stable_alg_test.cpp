#define BOOST_TEST_MODULE Stable_intersection_test
#include <boost/test/included/unit_test.hpp>
#include <magnet/intersection/stable_poly.hpp>
#include <magnet/math/polynomial.hpp>
#include <cmath>
#include <complex>

using namespace magnet::math;
const Polynomial<1> x{0, 1};

double rootvals[] = {-1e7, -1e6, -1e3, -100, -1, 0, 1, +100, 1e3, 1e6, 1e7 };

next_event(double t, std::initializer_list<)

BOOST_AUTO_TEST_CASE( Linear_function )
{
  for (root : rootvals) {
    auto poly = (x - root);
    
    //nextEvent()
  }
}
