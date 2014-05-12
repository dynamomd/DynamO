#define BOOST_TEST_MODULE Intersection_OffcenterSpheres_Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <magnet/intersection/generic_algorithm.hpp>
#include <magnet/intersection/polynomial.hpp>
#include <magnet/intersection/parabola_sphere.hpp>
#include <magnet/intersection/offcentre_spheres.hpp>
#include <magnet/math/matrix.hpp>
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

using namespace magnet::intersection;

BOOST_AUTO_TEST_CASE( OffCentreSphere_Test )
{
  RNG.seed(5489u);
  using namespace magnet::math;

  const Vector rij1(0.33930816635469108, 1.971007348602491, 0);
  const Vector vij(1.1608942531073687, -4.0757606085691398, 0);
  const Vector angvi(-0,-0,-1.0326096458374654);
  const Vector angvj(0,0,3.0759235803301794);
  const Vector relativeposi1(0.19838653763498912, -0.45895836596057499, 2.2204460492503128e-16);
  const Vector relativeposj1(0.32578919839301484, 0.37929065136177137, 0);
  const double diameteri=1, diameterj=1, maxdist = 2;

  magnet::intersection::detail::OffcentreSpheresOverlapFunction f1(rij1, vij, angvi, angvj, relativeposi1, relativeposj1, diameteri, diameterj, maxdist, 0, 0, 0, 0.49421681707429921);
  auto result1 = f1.nextEvent();
  //Check against verified result
  BOOST_CHECK(result1.first);
  BOOST_CHECK_CLOSE(result1.second, 0.032812502395565935, 1e-10);

  //Check that time shifts do not yield different roots
  const size_t tests = 1000000;
  for (size_t i(0); i < tests; ++i){
    const double dt = result1.second * dist01(RNG);
    const Vector rij2 = rij1 + dt * vij;
    const Vector relativeposi2 = Rodrigues(angvi * dt) * relativeposi1;
    const Vector relativeposj2 = Rodrigues(angvj * dt) * relativeposj1;
    magnet::intersection::detail::OffcentreSpheresOverlapFunction f2(rij2, vij, angvi, angvj, relativeposi2, relativeposj2, diameteri, diameterj, maxdist, 0, 0, 0, 0.81815864721356835);
    auto result2 = f2.nextEvent();
    
    BOOST_CHECK(result2.first);
    BOOST_CHECK_CLOSE(result2.second + dt, result1.second, 1e-10);
  }

  //Check that exceeding the predicted event time by a small amount
  //results in instant collisions.
  for (size_t i(0); i < tests; ++i){
    const double dt = result1.second * (1 + 0.01 * dist01(RNG));
    const Vector rij2 = rij1 + dt * vij;
    const Vector relativeposi2 = Rodrigues(angvi * dt) * relativeposi1;
    const Vector relativeposj2 = Rodrigues(angvj * dt) * relativeposj1;
    magnet::intersection::detail::OffcentreSpheresOverlapFunction f2(rij2, vij, angvi, angvj, relativeposi2, relativeposj2, diameteri, diameterj, maxdist, 0, 0, 0, 0.81815864721356835);
    auto result2 = f2.nextEvent();
    
    BOOST_CHECK(result2.first);
    BOOST_CHECK(result2.second == 0);
  }
}

BOOST_AUTO_TEST_CASE( OffCentreSphereGrowing_Test )
{

  const Vector r1(-58.502136713115014, -36.366193091936353, 0);
  const Vector v1(0.84416336042235662, 1.2211750596392101, 0);
  const Vector angv1(0, 0, 3.2360018350663511);
  const Quaternion orientation1(0.70710678118654757, Vector(0.51568947764362061, -0.48380198702325511, 0));
  const Vector relativepos1 = - 0.5 * (orientation1 * Quaternion::initialDirector());
  const double diam1 = 1.0;

  const Vector r2(-58.502136713115014, -33.203915431767982, 0);
  const Vector v2(1.1393960536121741, -0.73973475637700326, 0);
  const Vector angv2(0, 0, 1.2537004072984947);
  const Quaternion orientation2(0.70710678118654757, Vector(-0.61761216742478009, 0.34431847270930083, 0));
  const Vector relativepos2 = - 0.5 * (orientation2 * Quaternion::initialDirector());
  const double diam2 = 1.0;

  const double invgamma = 1.0;

  magnet::intersection::detail::OffcentreSpheresOverlapFunction f1(r1 - r2, v1 - v2, angv1, angv2, relativepos1, relativepos2, diam1, diam2, Vector(r1-r2).nrm() + 1.0, 0, invgamma, 0, 0.363046);
  auto result1 = f1.nextEvent();
  std::cout << "Result1.first = " << result1.first << "Result1.second = " << result1.second << std::endl;
  std::cout << "f = " << f1.eval(result1.second).front() << "Result1.second = " << result1.second << std::endl;
  
}
