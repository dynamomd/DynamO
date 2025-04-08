#define BOOST_TEST_MODULE Plane_Intersection_Tests
#include <boost/test/included/unit_test.hpp>
#include <magnet/intersection/ray_plane.hpp>
#include <random>

std::mt19937 RNG;
std::normal_distribution<double> normal_dist(0.0, 1.0);
std::uniform_real_distribution<double> angle_dist(0, M_PI);
std::uniform_real_distribution<double> dist01(0, 1);
using namespace magnet::math;

Vector random_vec() {
  return Vector{normal_dist(RNG), normal_dist(RNG), normal_dist(RNG)};
}

Vector random_unit_vec() {
  Vector vec = random_vec();
  return vec / vec.nrm();
}

const size_t testcount = 1000;
const double errlvl = 1e-8;

BOOST_AUTO_TEST_CASE(TimeToEvent_Test) {
  RNG.seed();

  for (size_t i(0); i < testcount; ++i) {
    // Generate a wall somewhere
    Vector n = random_unit_vec();
    Vector wallpos = random_vec();
    // Generate a particle somewhere
    Vector velocity = random_vec();
    Vector position = random_vec();
    double diam = std::abs(normal_dist(RNG));
    // Place it in contact with the wall
    position -= n * (position | n);
    position += diam * n + wallpos;
    if ((velocity | n) > 0)
      velocity = -velocity;

    // Now generate the position at some time away from the wall
    double deltat = std::abs(normal_dist(RNG));
    position += -deltat * velocity;

    // Test the collision is detected
    double calc_deltat =
        magnet::intersection::ray_plane(position - wallpos, velocity, n, diam);
    BOOST_CHECK_CLOSE(deltat, calc_deltat, errlvl);
  }
}

BOOST_AUTO_TEST_CASE(Overlapped_Approaching_Test) {
  RNG.seed();

  for (size_t i(0); i < testcount; ++i) {
    // Generate a wall somewhere
    Vector n = random_unit_vec();
    Vector wallpos = random_vec();
    // Generate a particle somewhere
    Vector velocity = random_vec();
    Vector position = random_vec();
    double diam = std::abs(normal_dist(RNG));
    // Place it in contact with the wall
    position -= n * (position | n);
    position += diam * n + wallpos;
    if ((velocity | n) > 0)
      velocity = -velocity;

    // Now generate the position 10% into the wall
    double deltat = 0.1 * diam / (-velocity | n);
    position += deltat * velocity;

    // Test the collision is detected as immediate
    double calc_deltat =
        magnet::intersection::ray_plane(position - wallpos, velocity, n, diam);
    BOOST_CHECK_EQUAL(calc_deltat, 0);
  }
}

BOOST_AUTO_TEST_CASE(Overlapped_Receeding_Test) {
  RNG.seed();

  for (size_t i(0); i < testcount; ++i) {
    // Generate a wall somewhere
    Vector n = random_unit_vec();
    Vector wallpos = random_vec();
    // Generate a particle somewhere
    Vector velocity = random_vec();
    Vector position = random_vec();
    double diam = std::abs(normal_dist(RNG));
    // Place it in contact with the wall
    position -= n * (position | n);
    position += diam * n + wallpos;
    if ((velocity | n) > 0)
      velocity = -velocity;

    // Now generate the position into the wall (but exiting the wall)
    double deltat = 1.01 * diam / (-velocity | n);
    position += deltat * velocity;

    // Test the collision is never detected
    double calc_deltat =
        magnet::intersection::ray_plane(position - wallpos, velocity, n, diam);
    BOOST_CHECK_EQUAL(calc_deltat, HUGE_VAL);
  }
}
