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
private:
  PolyGeneral(const magnet::intersection::detail::PolynomialFunction<Order>& f,
	      double t_min, double t_max): _f(f), _t_min(t_min), _t_max(t_max) {}

  template<size_t Deriv=0>
  double eval(double dt) const {
    return _f.template eval<Deriv>();
  }

  template<size_t Deriv=0>
  double max() const {
    double accum(0);
    for (size_t i(Deriv); i < Order + 1; ++i) {
      accum += std::max(std::pow(_t_min,i) * _f[i], std::pow(_t_max,i) * _f[i])  / factorial(i);
    }
    return accum;
  }

public:
  static constexpr uint64_t factorial(uint64_t n) { 
    return n == 0 ? 1  :  n * factorial(n-1); 
  }

  const magnet::intersection::detail::PolynomialFunction<Order>& _f;
  const double _t_min;
  const double _t_max;  
};


//template<bool inverse = false>
//inline double parabola_sphere(const math::Vector& R, const math::Vector& V, const math::Vector& A, const double& r)
//{
//  magnet::intersection::detail::PolynomialFunction<4> f{R.nrm2() - r * r, 2 * (V | R), 2 * (V.nrm2() + (A | R)), 6 * (A | V), 6 * A.nrm2()};
//  if (inverse) f.flipSign();
//  return detail::nextEvent(f, r * r);
//}



BOOST_AUTO_TEST_CASE( TimeToEvent_Test )
{
  RNG.seed();
}
