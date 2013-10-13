#define BOOST_TEST_MODULE Quaternion_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <magnet/math/quaternion.hpp>
#include <magnet/math/matrix.hpp>
#include <iostream>
#include <random>

std::mt19937 RNG;
std::normal_distribution<double> normal_dist(0, 1);
std::uniform_real_distribution<double> angle_dist(0, std::atan(1)*4);
std::uniform_real_distribution<double> dist01(0, 1);
using namespace magnet::math;

Vector random_unit_vec() {
  Vector vec(normal_dist(RNG), normal_dist(RNG), normal_dist(RNG));
  return vec/vec.nrm();
}

const size_t testcount = 1000;
const double errlvl = 1e-12;


BOOST_AUTO_TEST_CASE( Quaternion_fromToVector )
{
  RNG.seed();
  for (size_t i(0); i < testcount; ++i)
    {
      Vector start = random_unit_vec();
      Vector end = random_unit_vec();
      Vector err = end - (Quaternion::fromToVector(end, start) * start);
      
      BOOST_CHECK(std::abs(err[0]) < errlvl);
      BOOST_CHECK(std::abs(err[1]) < errlvl);
      BOOST_CHECK(std::abs(err[2]) < errlvl);
    }
}

BOOST_AUTO_TEST_CASE( Quaternion_fromAngleAxis )
{
  RNG.seed();
  for (size_t i(0); i < testcount; ++i)
    {
      double angle = angle_dist(RNG);
      Vector axis = random_unit_vec();
      Vector start = random_unit_vec();
      Vector end = Rodrigues(axis * angle) * start;
      
      Vector err = end - (Quaternion::fromAngleAxis(angle, axis) * start);
      
      BOOST_CHECK(std::abs(err[0]) < errlvl);
      BOOST_CHECK(std::abs(err[1]) < errlvl);
      BOOST_CHECK(std::abs(err[2]) < errlvl);
    }
}

BOOST_AUTO_TEST_CASE( Quaternion_toMatrix )
{
  RNG.seed();
  for (size_t i(0); i < testcount; ++i)
    {
      double angle = angle_dist(RNG);
      Vector axis = random_unit_vec();
      Vector start = random_unit_vec();
      Vector end = Rodrigues(axis * angle) * start;
      
      Vector err = end - (Quaternion::fromAngleAxis(angle, axis).toMatrix() * start);

      BOOST_CHECK(std::abs(err[0]) < errlvl);
      BOOST_CHECK(std::abs(err[1]) < errlvl);
      BOOST_CHECK(std::abs(err[2]) < errlvl);
    }
}

BOOST_AUTO_TEST_CASE( Quaternion_multiply )
{
  RNG.seed();
  for (size_t i(0); i < testcount; ++i)
    {
      Vector start = random_unit_vec();
      
      double angle1 = angle_dist(RNG);
      Vector axis1 = random_unit_vec();

      double angle2 = angle_dist(RNG);
      Vector axis2 = random_unit_vec();

      double angle3 = angle_dist(RNG);
      Vector axis3 = random_unit_vec();

      Vector end = Rodrigues(axis3 * angle3) * Rodrigues(axis2 * angle2) * Rodrigues(axis1 * angle1) * start;
      Vector err = end - (Quaternion::fromAngleAxis(angle3, axis3) * Quaternion::fromAngleAxis(angle2, axis2) * Quaternion::fromAngleAxis(angle1, axis1) * start);

      BOOST_CHECK(std::abs(err[0]) < errlvl);
      BOOST_CHECK(std::abs(err[1]) < errlvl);
      BOOST_CHECK(std::abs(err[2]) < errlvl);
    }
}

BOOST_AUTO_TEST_CASE( Quaternion_inverse )
{
  RNG.seed();
  for (size_t i(0); i < testcount; ++i)
    {
      Vector start = random_unit_vec();

      double angle1 = angle_dist(RNG);
      Vector axis1 = random_unit_vec();

      double angle2 = angle_dist(RNG);
      Vector axis2 = random_unit_vec();

      Vector end = Rodrigues(axis1 * angle1) * start;
      Vector err = end - (Quaternion::fromAngleAxis(angle1, axis1) * Quaternion::fromAngleAxis(angle2, axis2) * Quaternion::fromAngleAxis(angle2, axis2).inverse() * start);

      BOOST_CHECK(std::abs(err[0]) < errlvl);
      BOOST_CHECK(std::abs(err[1]) < errlvl);
      BOOST_CHECK(std::abs(err[2]) < errlvl);
    }
}

BOOST_AUTO_TEST_CASE( GLSL_rotation_formula )
{
  RNG.seed();
  for (size_t i(0); i < testcount; ++i)
    {
      Vector start = random_unit_vec();
      Vector end = random_unit_vec();
      Quaternion q = Quaternion::fromToVector(end, start);

      Vector result = start + 2.0 * (q.imaginary() ^ ((q.imaginary() ^ start) + q.real() * start));
      Vector err = end - result;

      BOOST_CHECK(std::abs(err[0]) < errlvl);
      BOOST_CHECK(std::abs(err[1]) < errlvl);
      BOOST_CHECK(std::abs(err[2]) < errlvl);
    }
}
