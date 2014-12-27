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

const std::array<double, 10> rootvals{-1e7, -1e3, -3.14159265, -1, 0, 1, 3.14159265, +100, 1e3, 1e7 };
const size_t tests = 1000;

template<class F, class R>
void test_solution(const F& f, double tol, R actual_roots) {
  //Counter of how many tests have been done (to give an idea of where failures are)
  static size_t counter = 0;
  ++counter;

  Variable<'t'> t;

  const auto df = derivative(f, t);
  double solution;
  double nextroot;
  decltype(solve_real_roots(f)) roots;
  decltype(solve_real_roots(df)) droots;
  
  nextroot = HUGE_VAL;
  solution = magnet::intersection::nextEvent(f);
  try {
  roots = solve_real_roots(f);
  droots = solve_real_roots(df);

  std::sort(actual_roots.begin(), actual_roots.end());

  //Multiplicity of each root
  std::map<double, size_t> root_counters;
  for (auto root : actual_roots)
    root_counters[root] = 0;
  for (auto root : actual_roots)
    ++root_counters[root];


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
    
      //Check if the event corresponds to a root in the derivative (turning point) before the next root
      double next_droot = HUGE_VAL;
      for (auto droot : droots)
	if (droot > 0)
	  next_droot = droot;

      //Find the first positive root
      while (it->first < 0) ++it;

      if (!std::isinf(solution) && !std::isinf(next_droot) 
	  && !::boost::test_tools::check_is_close(it->first, next_droot, ::boost::test_tools::percent_tolerance(tol))
	  && (next_droot < it->first)) {
	//It should be a turning point root
	nextroot = next_droot;
	if (!::boost::test_tools::check_is_close(solution, next_droot, ::boost::test_tools::percent_tolerance(tol)))
	  throw std::runtime_error("Turning point root missed?");
	else
	  return;
      }
      
      for (; it != root_counters.end(); ++it) {
	nextroot = it->first;

	bool is_odd_root = it->second % 2;
	if (!is_odd_root) {
	  //Its an even root, we must detect this!
	  if (!std::isinf(solution) && !::boost::test_tools::check_is_close(solution, nextroot, ::boost::test_tools::percent_tolerance(tol)))
	    throw std::runtime_error("Missed a turnback root?");
	  else
	    return;
	} else {
	  //Odd root. This means the particle eventually passes
	  //through this root, to outside. 
	  
	  if (it->second > 1) {
	    //distinct roots should always pass through the
	    //transition. However, odd multiplicity roots may
	    //numerically turn around for an instant. We have to
	    //accept this as a recoverable error.
	    if (!std::isinf(solution) && ::boost::test_tools::check_is_close(solution, it->first, ::boost::test_tools::percent_tolerance(tol)))
	      return;
	  }

	  //We got through and back into non-overlapped space. Carry
	  //on with the normal event detection.
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
	if (!std::isinf(solution) && !::boost::test_tools::check_is_close(solution, nextroot, ::boost::test_tools::percent_tolerance(tol))
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
  feenableexcept(FE_INVALID | FE_OVERFLOW);

  for (double sign : {-1.0, +1.0})
    for (double root : rootvals) {
      auto poly = (t - root) * sign;
      
      std::uniform_real_distribution<double> shift_dist(-10, 10);
      for (size_t i(0); i < tests; ++i) {
	double shift = shift_dist(RNG);
	auto s_poly = shift_function(poly, shift);
	test_solution(s_poly, 1e-10, magnet::containers::StackVector<double,1>{root - shift});
      }
    }
}

BOOST_AUTO_TEST_CASE( Quadratic_function )
{
  RNG.seed(1);
  feenableexcept(FE_INVALID | FE_OVERFLOW);

  for (double sign : {-1.0, +1.0})
    for (auto it1 = rootvals.begin(); it1 != rootvals.end(); ++it1)
      for (auto it2 = it1; it2 != rootvals.end(); ++it2) {
	auto poly = (t - *it1) * (t - *it2) * sign;
	std::uniform_real_distribution<double> shift_dist(-10, 10);
	for (size_t i(0); i < tests; ++i) {
	  double shift = shift_dist(RNG);
	  auto s_poly = shift_function(poly, shift);
	  test_solution(s_poly, 1e-2, magnet::containers::StackVector<double,2>{*it1 - shift, *it2 - shift});
	}
      }
}

BOOST_AUTO_TEST_CASE( Cubic_function )
{
  RNG.seed(1);

  const std::array<double, 10> rootvals{-1e7, -1e3, -3.14159265, -1, 0, 1, 3.14159265, +100, 1e3, 1e7 };
  
  for (auto it1 = rootvals.begin(); it1 != rootvals.end(); ++it1)
    for (auto it2 = it1; it2 != rootvals.end(); ++it2)
      for (auto it3 = it2; it3 != rootvals.end(); ++it3)
	for (double sign : {-1.0, +1.0})
	  {
	    auto poly = (t - *it1) * (t - *it2) * (t - *it3) * sign;
	    std::uniform_real_distribution<double> shift_dist(-10, 10);
	    for (size_t i(0); i < tests; ++i) {
	      double shift = shift_dist(RNG);
	      auto s_poly = shift_function(poly, shift);
	      test_solution(s_poly, 1e-1, magnet::containers::StackVector<double,3>{*it1 - shift, *it2 - shift, *it3 - shift});
	    }
	  }
}


BOOST_AUTO_TEST_CASE( Quartic_function )
{
  RNG.seed(1);
  const std::array<double, 10> rootvals{-1e7, -1e3, -3.14159265, -1, 0, 1, 3.14159265, +100, 1e3, 1e7 };
  
  for (auto it1 = rootvals.begin(); it1 != rootvals.end(); ++it1)
    for (auto it2 = it1+1; it2 != rootvals.end(); ++it2)
      for (auto it3 = it2+1; it3 != rootvals.end(); ++it3)
	for (auto it4 = it3+1; it4 != rootvals.end(); ++it4)
	  for (double sign : {-1.0, +1.0})
	    {
	      auto poly = (t - *it1) * (t - *it2) * (t - *it3) * (t - *it4) * sign;
	      std::uniform_real_distribution<double> shift_dist(-10, 10);
	      for (size_t i(0); i < tests; ++i) {
		double shift = shift_dist(RNG);
		auto s_poly = shift_function(poly, shift);
		test_solution(s_poly, 1e-1, magnet::containers::StackVector<double,4>{*it1 - shift, *it2 - shift, *it3 - shift, *it4 - shift});
	      }
	    }
}

BOOST_AUTO_TEST_CASE( Quintic_function )
{
  RNG.seed(1);
  const std::array<double, 10> rootvals{-1e7, -1e3, -3.14159265, -1, 0, 1, 3.14159265, +100, 1e3, 1e7 };
  
  for (auto it1 = rootvals.begin(); it1 != rootvals.end(); ++it1)
    for (auto it2 = it1+1; it2 != rootvals.end(); ++it2)
      for (auto it3 = it2+1; it3 != rootvals.end(); ++it3)
	for (auto it4 = it3+1; it4 != rootvals.end(); ++it4)
	  for (auto it5 = it4+1; it5 != rootvals.end(); ++it5)
	    for (double sign : {-1.0, +1.0})
	      {
		auto poly = (t - *it1) * (t - *it2) * (t - *it3) * (t - *it4) * (t - *it5) * sign;
	      std::uniform_real_distribution<double> shift_dist(-10, 10);
	      for (size_t i(0); i < tests; ++i) {
		double shift = shift_dist(RNG);
		auto s_poly = shift_function(poly, shift);
		test_solution(s_poly, 1e-1, magnet::containers::StackVector<double,5>{*it1 - shift, *it2 - shift, *it3 - shift, *it4 - shift, *it5 - shift});
	      }
	    }
}
