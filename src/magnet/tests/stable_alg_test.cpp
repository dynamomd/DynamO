#define BOOST_TEST_MODULE Stable_intersection_test
#include <boost/test/included/unit_test.hpp>
#include <magnet/intersection/stable_poly.hpp>
#include <magnet/math/polynomial.hpp>
#include <cmath>
#include <complex>
#include <random>
#include <fenv.h>

using namespace magnet::math;
const Polynomial<1, double, 't'> t{0, 1};
std::mt19937 RNG;

const double rootvals[] = {-1e7, -1e3, -3.14159265, -1, 0, 1, 3.14159265, +100, 1e3, 1e7 };
const size_t tests = 1000;

template<class F, class R>
void test_solution(const F& f, double tol, R actual_roots) {
  //Counter of how many tests have been done (to give an idea of where failures are)
  static size_t counter = 0;
  ++counter;

  std::sort(actual_roots.begin(), actual_roots.end());

//  if (counter == 222718)
//    std::cout << "PING" << std::endl;

  double nextroot = HUGE_VAL;
  const double solution = magnet::intersection::nextEvent(f);
  Variable<'t'> t;
  const auto roots = solve_roots(f);
  const auto df = derivative(f, t);
  const auto droots = solve_roots(df);

  //Multiplicity of each root
  std::map<double, size_t> root_counters;
  for (auto root : actual_roots)
    root_counters[root] = 0;
  for (auto root : actual_roots)
    ++root_counters[root];


  try {

    if (solution == 0) {
      //Immediate collision! Check that the particles are currently
      //approaching and overlapping
      if (eval(f, t==0.0) > 0)
	throw std::runtime_error("Not sufficiently overlapped during an immediate collision");
      if (eval(df, t==0.0) > 0)
	throw std::runtime_error("Not sufficiently approaching during an immediate collision");
      return;
    }

    auto it = root_counters.begin();

    if (f[0] < 0) {
      //Particles have started out overlapping
    
      //Check if the event corresponds to a root in the derivative (turning point)
      for (auto droot : droots)
	if (!std::isinf(solution) && ::boost::test_tools::check_is_close(solution, droot, ::boost::test_tools::percent_tolerance(tol))) {
	  //Its a turning point, check there are no roots before it
	  for (auto root : actual_roots)
	    if ((root > 0) && (root < droot) && !::boost::test_tools::check_is_close(root, droot, ::boost::test_tools::percent_tolerance(tol)))
	      throw std::runtime_error("Root between turnaround and current time");
	  return;
	}

      for (; it != root_counters.end(); ++it) {
	//Check that this root is in the future
	if (it->first < 0) continue;

	nextroot = it->first;
	bool is_odd_root = it->second % 2;
	if (!is_odd_root) {
	  //Its an even root, we must detect this!
	  if (!std::isinf(solution) && !::boost::test_tools::check_is_close(solution, nextroot, ::boost::test_tools::percent_tolerance(tol)))
	    throw std::runtime_error("Missed a turnback root?");
	  else
	    return;
	} else {
	  //Odd root. This means the particle passes through
	  //this root, to outside. Carry on with the normal detection!
	  ++it;
	  break;
	}
      }
    }
  
    //We expect the next, non-negative, odd multiplicity root to be an
    //event. even-multiplicity roots will be accepted though (as they
    //may be numerically distinct)
  
    for (; it != root_counters.end(); ++it) {
      if (it->first < 0) continue;
      nextroot = it->first;
      bool is_odd_root = it->second % 2;
      
      if (is_odd_root) {
	//Its an odd root, we must detect this!
	if (!::boost::test_tools::check_is_close(solution, nextroot, ::boost::test_tools::percent_tolerance(tol))
	    && (std::abs(eval(f, t==solution)) > 4 * precision(f, solution)))
	  throw std::runtime_error("Missed a root?");
	else 
	  return;
      } else if (!std::isinf(solution) && ::boost::test_tools::check_is_close(solution, nextroot, ::boost::test_tools::percent_tolerance(tol))) 
	//Detected an even root as a crossing, not a huge error.
	return;
      else
	continue;
    }

    if (solution != HUGE_VAL)
      throw std::runtime_error("Detected an extra root?");

  } catch (std::exception& e) {
    BOOST_ERROR(e.what());
    std::cout.precision(50);
    std::cout << "TEST " << counter << std::endl;
    std::cout << "f(x)=" << f << std::endl;
    std::cout << "f'(x)=" << df << std::endl;
    std::cout << "f''(x)=" << derivative(df, t) << std::endl;
    std::cout << "f(0)=" << eval(f, t==0) << std::endl;
    std::cout << "f'(0)=" << eval(df, t==0) << std::endl;
    std::cout << "solution = " << solution << std::endl;
    std::cout << "f(solution)=" << eval(f, t==solution) << std::endl;
    std::cout << "f'(solution)=" << eval(df, t==solution) << std::endl;
    std::cout << "f(nextroot="<<nextroot<<")=" << eval(f, t==nextroot) << std::endl;
    std::cout << "f'(nextroot="<<nextroot<<")=" << eval(df, t==nextroot) << std::endl;
    std::cout << "actual_roots = " << actual_roots << std::endl;
    std::cout << "roots = " << roots << std::endl;
    std::cout << "f' roots = " << droots << std::endl;
    std::cout << "d|f|(nextroot) = " << precision(f, nextroot) << std::endl;
    std::cout << "d|f'|(nextroot) = " << precision(df, nextroot) << std::endl;
    std::cout << "d|f|(next_event) = " << precision(f, solution) << std::endl;
    std::cout << "d|f'|(next_event) = " << precision(df, solution) << std::endl;
    throw;
  }
}

BOOST_AUTO_TEST_CASE( Linear_function )
{
  RNG.seed(1);

  for (double sign : {-1.0, +1.0})
    for (double root : rootvals) {
      auto poly = (t - root) * sign;
      
      std::uniform_real_distribution<double> shift_dist(-10, 10);
      for (size_t i(0); i < tests; ++i) {
	double shift = shift_dist(RNG);
	auto s_poly = shift_function(poly, shift);
	auto roots = solve_roots(s_poly);
	test_solution(s_poly, 1e-10, magnet::containers::StackVector<double,1>{root - shift});
      }
    }
}

BOOST_AUTO_TEST_CASE( Quadratic_function )
{
  RNG.seed(1);
  for (double sign : {-1.0, +1.0})
    for (double root1 : rootvals)
      for (double root2 : rootvals) {
	auto poly = (t - root1) * (t - root2) * sign;
	std::uniform_real_distribution<double> shift_dist(-10, 10);
	for (size_t i(0); i < tests; ++i) {
	  double shift = shift_dist(RNG);
	  auto s_poly = shift_function(poly, shift);
	  auto roots = solve_roots(s_poly);
	  test_solution(s_poly, 1e-2, magnet::containers::StackVector<double,2>{root1 - shift, root2 - shift});
	}
      }
}

BOOST_AUTO_TEST_CASE( Cubic_function )
{
  RNG.seed(1);
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  const double rootvals[] = {-1e7, -1e3, -3.14159265, -1, 0, 1, 3.14159265, +100, 1e3, 1e7 };

  std::cout.precision(50);  
  auto poly = (t*t*t + -7.9258638382968493729663350677583366632461547851562 * t*t + 19.41096616822921561151815694756805896759033203125 * t + -13.674032423594180585268986760638654232025146484375);
  auto roots = solve_roots(poly);
  
  for (double sign : {-1.0, +1.0})
    for (double root1 : rootvals)
      for (double root2 : rootvals) 
	for (double root3 : rootvals) 
	  {
	    auto poly = (t - root1) * (t - root2) * (t - root3) * sign;
	    std::uniform_real_distribution<double> shift_dist(-10, 10);
	    for (size_t i(0); i < tests; ++i) {
	      double shift = shift_dist(RNG);
	      auto s_poly = shift_function(poly, shift);
	      test_solution(s_poly, 1e-1, magnet::containers::StackVector<double,3>{root1 - shift, root2 - shift, root3 - shift});
	    }
	  }
}

//BOOST_AUTO_TEST_CASE( Quartic_function )
//{
//  RNG.seed(1);
//  for (double sign : {-1.0, +1.0})
//    for (double root1 : rootvals)
//      for (double root2 : rootvals) 
//	for (double root3 : rootvals) 
//	  for (double root4 : rootvals) 
//	    {
//	      auto poly = (x - root1) * (x - root2) * (x - root3) * (x - root4) * sign;
//	      std::uniform_real_distribution<double> shift_dist(-10, 10);
//	      for (size_t i(0); i < tests; ++i) {
//		double shift = shift_dist(RNG);
//		auto s_poly = shift_polynomial(poly, shift);
//		test_solution(s_poly, 1e-4);
//	      }
//	    }
//}
