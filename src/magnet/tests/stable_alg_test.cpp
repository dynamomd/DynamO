#define BOOST_TEST_MODULE Stable_intersection_test
#include <boost/test/included/unit_test.hpp>
#include <magnet/intersection/stable_poly.hpp>
#include <magnet/math/polynomial.hpp>
#include <magnet/intersection/stable_poly.hpp>
#include <cmath>
#include <complex>
#include <random>

using namespace magnet::math;
const Polynomial<1> x{0, 1};
std::mt19937 RNG;

const double rootvals[] = {-1e7, -1e6, -1e3, -3.14159265, -1, 0, 1, 3.14159265, +100, 1e3, 1e6, 1e7 };
const size_t tests = 1000;

template<class F, class DF, class R>
void test_solution(const F& f, const DF& df, double solution, const R& roots) {
  if (solution == HUGE_VAL) {
    //No root! so check there really are no roots
    for (double root : roots)
      if (root > 0) {
	BOOST_CHECK(eval(df, root) > 0);
      }
  } else if (solution == 0) {
    //Immediate collision! Check that the particles are currently
    //approaching and overlapping
    BOOST_CHECK(eval(f, 0) <= 0);
    BOOST_CHECK(eval(df, 0) < 0);
  } else {
    //There is a root. First determine what the next root actually is
    double nextroot = HUGE_VAL;
    for (double root : roots)
      if (root > 0)
	nextroot = std::min(nextroot, root);

    if (eval(f, 0) > 0) {
      //Particles started out not overlapping, therefore the solution
      //must be an actual root
      if (nextroot == HUGE_VAL) {
	//Check if this is a phantom root
	BOOST_CHECK(nextroot == HUGE_VAL);
      } else {	
	BOOST_CHECK_CLOSE(solution, nextroot, 1e-10);

	if (std::abs((solution - nextroot) / nextroot) > 1e-10) {
	  std::cout << "f(x)=" << f << std::endl;
	  std::cout << "f'(x)=" << df << std::endl;
	  std::cout << "f("<< solution <<")=" << eval(f, solution) << std::endl;
	  std::cout << roots << std::endl;
	}
      }
    } else /*if (f(0) <= 0)*/ {
      //The particles started out overlapping, check if the detected
      //root is still overlapping
      if (eval(f, solution) < 0)
	//This is an approaching root at a turning point
	BOOST_CHECK_SMALL(eval(df, solution), 1e-6);
      else {
	//The detected root must be after the nextroot (which is the
	//exit root)
	BOOST_CHECK(solution > nextroot);
	BOOST_CHECK_SMALL(eval(f, solution), 1e-8);
	double err = HUGE_VAL;
	for (double root : roots)
	  err = std::min(err, std::abs(root - solution));
	BOOST_CHECK_SMALL(err, 1e-8);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE( Linear_function )
{
  RNG.seed(1);

  for (double root : rootvals) {
    auto poly = (x - root);
    
    std::uniform_real_distribution<double> shift_dist(-10, 10);
    for (size_t i(0); i < tests; ++i) {
      auto s_poly = shift_polynomial(poly, shift_dist(RNG));
      auto roots = solve_roots(s_poly);
      test_solution(s_poly, derivative(s_poly), magnet::intersection::nextEvent(s_poly), roots);
    }
  }
}

BOOST_AUTO_TEST_CASE( Quadratic_function )
{
  RNG.seed(1);

  for (double root1 : rootvals)
    for (double root2 : rootvals) {
      auto poly = (x - root1) * (x - root2);
      
      std::uniform_real_distribution<double> shift_dist(-10, 10);
      for (size_t i(0); i < tests; ++i) {
	auto s_poly = shift_polynomial(poly, shift_dist(RNG));
	auto roots = solve_roots(s_poly);
	test_solution(s_poly, derivative(s_poly), magnet::intersection::nextEvent(s_poly), roots);
      }
    }
}
