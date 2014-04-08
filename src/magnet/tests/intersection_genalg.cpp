#define BOOST_TEST_MODULE Intersection_GeneralAlgorithm_Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <magnet/intersection/generic_algorithm.hpp>
#include <magnet/intersection/polynomial.hpp>
#include <magnet/intersection/parabola_sphere.hpp>
#include <iostream>
#include <random>

std::mt19937 RNG;
std::normal_distribution<double> normal_dist(0.0, 1.0);
std::uniform_real_distribution<double> angle_dist(0, M_PI);
std::uniform_real_distribution<double> dist01(0, 1);
using namespace magnet::math;

Vector random_vec() {
  return Vector(normal_dist(RNG), normal_dist(RNG), normal_dist(RNG));
}

Vector random_unit_vec() {
  Vector vec = random_vec();
  return vec / vec.nrm();
}

const size_t testcount = 1000;
const double errlvl = 1e-8;


using namespace magnet::intersection;

template<size_t Order>
class PolyGeneral {
public:
  PolyGeneral(const magnet::intersection::detail::PolynomialFunction<Order>& f,
	      double t_min, double t_max): _f(f), _t_min(t_min), _t_max(t_max) {}

  template<size_t Deriv=0>
  double eval(double dt) const {
    return _f.template eval<Deriv>(dt);
  }

  template<size_t Deriv=0>
  double max() const {
    double accum(0);
    for (size_t i(Deriv); i < Order + 1; ++i) {
      accum += std::max(std::pow(_t_min,i) * _f[i], std::pow(_t_max,i) * _f[i])  / factorial(i);
    }
    return accum;
  }

private:
  static constexpr uint64_t factorial(uint64_t n) { 
    return n == 0 ? 1  :  n * factorial(n-1); 
  }

  const magnet::intersection::detail::PolynomialFunction<Order>& _f;
  const double _t_min;
  const double _t_max;  
};

BOOST_AUTO_TEST_CASE( GravitySphere_Test )
{
  RNG.seed(5489u);
  const size_t tests = 100000;

  for (size_t i(0); i < tests; ++i){
    const Vector aij = random_unit_vec();
    const Vector rij = random_unit_vec() * 1.5;
    const Vector vij = random_vec();
    const double r = 1;
  
    magnet::intersection::detail::PolynomialFunction<4> f_radical{rij.nrm2() - r * r, 2 * (vij | rij), 2 * (vij.nrm2() + (aij | rij)), 6 * (aij | vij), 6 * aij.nrm2()};
    double radical_root = magnet::intersection::detail::nextEvent(f_radical, r * r);

    PolyGeneral<4> f_numerical(f_radical, 0, 10);

    const double t_max = 10;
    auto numerical_root = magnet::intersection::nextEvent(f_numerical, 0, t_max);
    while (!numerical_root.first && !std::isinf(numerical_root.second))
      numerical_root = magnet::intersection::nextEvent(f_numerical, numerical_root.second, t_max);
    
    if (std::isinf(radical_root) != std::isinf(numerical_root.second))
      { BOOST_CHECK(false); }
    else if (!std::isinf(radical_root))
      { BOOST_CHECK_CLOSE(radical_root, numerical_root.second, 1e-12); }
  }
}
