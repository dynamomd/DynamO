#define BOOST_TEST_MODULE Plane_Intersection_Tests
#include <boost/test/included/unit_test.hpp>
#include <magnet/intersection/ray_plane.hpp>
#include <iostream>
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

std::pair<double, double> ray_wall_bounds(const Vector& relPos, const Vector& relVel, 
				   const Vector& normal, const double dist)
{
  using namespace std;
  //Objects must move towards each other to intersect
  const double r = (relPos | normal);
  const double v = (relVel | normal);
  const double val1 = max(-(r - dist) / v, 0.0);
  const double val2 = max(-(r + dist) / v, 0.0);
  return pair<double, double>(min(val1, val2), max(val1, val2));
}

double ray_prism(const Vector V[3], const Vector& pos, 
		 const Vector& vel, const double dist)
{
//  using namespace std;
//  const Vector M = (V[0] + V[1] + V[2]) / 3;
//  Vector Vrel[3];
//  double Vdist[3];
//
//  for (size_t i(0); i<3; ++i) {
//    Vrel[i] = V[i] - M;
//    Vdist[i] = Vrel[i].nrm();
//    Vrel[i] /= Vdist[i];
//  }
//
//  double tmin = HUGE_VAL;
//  double tmax = - HUGE_VAL;
//  Vector relPos = pos - M;
//
//  //Test against the bounds given by the three walls
//  for (size_t i(0); i<3; ++i) {
//    auto bounds = ray_wall(relPos, vel, Vrel[i], Vdist[i]);
//    tmin = min(bounds.first, tmin);
//    tmax = max(bounds.second, tmax);
//  }
//  
//  //Test against the plane the triangle is in
//  Vector N = Vrel[0] ^ Vrel[1];
//  N /= N.nrm();
//  auto bounds = ray_wall(relPos, vel, N, dist);
//  tmin = min(bounds.first, tmin);
//  tmax = max(bounds.second, tmax);
//
  return HUGE_VAL;
}


BOOST_AUTO_TEST_CASE(TimeToEvent_Test)
{
//  const Vector V[] = {random_vec(), random_vec(), random_vec()};
//
//  RNG.seed();
//  auto out = ray_wall(Vector(1.5141, 0, 0), Vector(-0.92378417, 0, 0), Vector(-1, 0, 0), 0.5);
//  std::cout << "rootpair = " << out.first << ", " << out.second << std::endl;
//
//  const size_t testcount = 1000;
//  for (size_t i(0); i < testcount; ++i)
//    {
//      //Generate a wall somewhere
//      Vector n = random_unit_vec();
//      Vector wallpos = random_vec();
//      //Generate a particle somewhere
//      Vector velocity = random_vec();
//      Vector position = random_vec();
//      double diam = std::abs(normal_dist(RNG));
//      //Place it in contact with the wall
//      position -= n * (position | n);
//      position += diam * n + wallpos;
//      if ((velocity | n) > 0)
//	velocity = -velocity;
//
//      //Now generate the position at some time away from the wall
//      double deltat = std::abs(normal_dist(RNG));
//      position += -deltat * velocity;
//
//      //Test the collision is detected
//      double calc_deltat = magnet::intersection::ray_plane(position - wallpos, velocity, n, diam);
//      const double errlvl = 1e-8;
//      BOOST_CHECK_CLOSE(deltat, calc_deltat, errlvl);
//    }
}
