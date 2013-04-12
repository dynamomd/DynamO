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

int main(int argc, char *argv[])
{
  RNG.seed();
  double errlvl = 1e-12;
  size_t testcount = 100000;
  std::cout << "Quaternion Testing:" << std::endl;

  std::cout << "Quaternion::fromToVector()" << std::endl;
  for (size_t i(0); i < testcount; ++i)
    {
      Vector start = random_unit_vec();
      Vector end = random_unit_vec();
      Vector err = end - (Quaternion::fromToVector(end, start) * start);

      if ((std::abs(err[0]) > errlvl) || (std::abs(err[1]) > errlvl) || (std::abs(err[2]) > errlvl))
	{
	  std::cout << "test " << i << ": error " << err.toString() << std::endl;
	  return EXIT_FAILURE;
	}
    }

  std::cout << "Quaternion::fromAngleAxis()" << std::endl;
  for (size_t i(0); i < testcount; ++i)
    {
      double angle = angle_dist(RNG);
      Vector axis = random_unit_vec();
      Vector start = random_unit_vec();
      Vector end = Rodrigues(axis * angle) * start;
      
      Vector err = end - (Quaternion::fromAngleAxis(angle, axis) * start);
      
      if ((std::abs(err[0]) > errlvl) || (std::abs(err[1]) > errlvl) || (std::abs(err[2]) > errlvl))
	{
	  std::cout << "test " << i << ": error " << err.toString() << std::endl;
	  return EXIT_FAILURE;
	}
    }

  std::cout << "Quaternion::toMatrix()" << std::endl;
  for (size_t i(0); i < testcount; ++i)
    {
      double angle = angle_dist(RNG);
      Vector axis = random_unit_vec();
      Vector start = random_unit_vec();
      Vector end = Rodrigues(axis * angle) * start;
      
      Vector err = end - (Quaternion::fromAngleAxis(angle, axis).toMatrix() * start);
      
      if ((std::abs(err[0]) > errlvl) || (std::abs(err[1]) > errlvl) || (std::abs(err[2]) > errlvl))
	{
	  std::cout << "test " << i << ": error " << err.toString() << std::endl;
	  return EXIT_FAILURE;
	}
    }

  std::cout << "Quaternion::operator*(Quaternion&)" << std::endl;
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
      
      Vector err = end - (Quaternion::fromAngleAxis(angle3, axis3) 
			  * Quaternion::fromAngleAxis(angle2, axis2)
			  * Quaternion::fromAngleAxis(angle1, axis1)
			  * start);
      
      if ((std::abs(err[0]) > errlvl) || (std::abs(err[1]) > errlvl) || (std::abs(err[2]) > errlvl))
	{
	  std::cout << "test " << i << ": error " << err.toString() << std::endl;
	  return EXIT_FAILURE;
	}
    }

  std::cout << "Quaternion::inverse()" << std::endl;
  for (size_t i(0); i < testcount; ++i)
    {
      Vector start = random_unit_vec();

      double angle1 = angle_dist(RNG);
      Vector axis1 = random_unit_vec();

      double angle2 = angle_dist(RNG);
      Vector axis2 = random_unit_vec();

      Vector end = Rodrigues(axis1 * angle1) * start;
      
      Vector err = end - (Quaternion::fromAngleAxis(angle1, axis1) 
			  * Quaternion::fromAngleAxis(angle2, axis2)
			  * Quaternion::fromAngleAxis(angle2, axis2).inverse()
			  * start);
      
      if ((std::abs(err[0]) > errlvl) || (std::abs(err[1]) > errlvl) || (std::abs(err[2]) > errlvl))
	{
	  std::cout << "test " << i << ": error " << err.toString() << std::endl;
	  return EXIT_FAILURE;
	}
    }

  std::cout << "GLSL rotation()" << std::endl;
  for (size_t i(0); i < testcount; ++i)
    {
      Vector start = random_unit_vec();
      Vector end = random_unit_vec();
      Quaternion q = Quaternion::fromToVector(end, start);

      Vector result = start + 2.0 * (q.imaginary() ^ ((q.imaginary() ^ start) + q.real() * start));
      Vector err = end - result;

      if ((std::abs(err[0]) > errlvl) || (std::abs(err[1]) > errlvl) || (std::abs(err[2]) > errlvl))
	{
	  std::cout << "test " << i << ": error " << err.toString() << std::endl;
	  return EXIT_FAILURE;
	}
    }

  return EXIT_SUCCESS;
}
