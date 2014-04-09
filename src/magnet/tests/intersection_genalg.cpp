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

BOOST_AUTO_TEST_CASE( OffCentreSphere_Test )
{
  RNG.seed(5489u);

  const long double systime1 = 0.29228171769001683291721097046611533;
  const Vector rij1(0.33930816635469108, 1.971007348602491, 0);
  const Vector vij(1.1608942531073687, -4.0757606085691398, 0);
  const Vector angvi(-0,-0,-1.0326096458374654);
  const Vector angvj(0,0,3.0759235803301794);
  const Vector relativeposi1(0.19838653763498912, -0.45895836596057499, 2.2204460492503128e-16);
  const Vector relativeposj1(0.32578919839301484, 0.37929065136177137, 0);
  double diameteri=1, diameterj=1, maxdist = 2;

  const long double systime2 = 0.32486537488833596663004646409866893;
  const double dt = (systime2 - systime1)*0.99999;
  const Vector rij2 = rij1 + dt * vij;
  const Vector relativeposi2 = Rodrigues(angvi * dt) * relativeposi1;
  const Vector relativeposj2 = Rodrigues(angvj * dt) * relativeposj1;

//  std::cout << "\nrij2=" << rij2.toString() 
//	    << "\nrposi2=" << relativeposi2.toString()
//	    << "\nrposj2=" << relativeposj2.toString()
//	    << std::endl;

  magnet::intersection::detail::OffcentreSpheresOverlapFunction f1(rij1, vij, angvi, angvj, relativeposi1, relativeposj1, diameteri, diameterj, 2);
  auto result1 = magnet::intersection::nextEvent(f1, 0, 0.49421681707429921);
  //0.032812502395565935

  std::cout << "\nresult1="<< result1.first << "," << result1.second - dt << std::endl;

  magnet::intersection::detail::OffcentreSpheresOverlapFunction f2(rij2, vij, angvi, angvj, relativeposi2, relativeposj2, diameteri, diameterj, 2);
  auto result2 = magnet::intersection::nextEvent(f2, 0, 0.81815864721356835);

  std::cout << "\nresult2="<< result2.first << "," << result2.second << std::endl;
}
