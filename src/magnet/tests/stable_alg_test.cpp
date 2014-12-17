#define BOOST_TEST_MODULE Stable_intersection_test
#include <boost/test/included/unit_test.hpp>
#include <magnet/intersection/stable_poly.hpp>
#include <magnet/math/polynomial.hpp>
#include <cmath>
#include <complex>
#include <random>

using namespace magnet::math;
const Polynomial<1, double, 't'> t{0, 1};
std::mt19937 RNG;

const double rootvals[] = {-1e7, -1e3, -3.14159265, -1, 0, 1, 3.14159265, +100, 1e3, 1e7 };
const size_t tests = 1000;

template<class F>
void test_solution(const F& f, double solution, double tol) {
  Variable<'t'> t;
  const auto roots = solve_roots(f);
  const auto df = derivative(f, t);
  const auto droots = solve_roots(df);
  if (solution == HUGE_VAL) {
    //No root! so check there really are no roots
    for (double root : roots)
      if (root > 0) {
	if ((eval(df, t == root) < 0) || ((eval(df, t== root) == 0) && (eval(derivative(df, t), t == root) < 0))) {
	  std::cout.precision(30);
	  std::cout << "f(x)=" << f << std::endl;
	  std::cout << "f'(x)=" << df << std::endl;
	  std::cout << "f''(x)=" << derivative(df, t) << std::endl;
	  std::cout << "f(0)=" << eval(f, t == 0) << std::endl;
	  std::cout << "f'(0)=" << eval(df, t == 0) << std::endl;
	  std::cout << "f("<< solution <<")=" << eval(f, t == solution) << std::endl;
	  std::cout << "f("<< root <<")=" << eval(f, root) << std::endl;
	  std::cout << "f'("<< root <<")=" << eval(df, root) << std::endl;
	  std::cout << "roots = " << roots << std::endl;
	  std::cout << "f' roots = " << droots << std::endl;
	  BOOST_ERROR("Did not detect a root!");
	}
      }
  } else if (solution == 0) {
    //Immediate collision! Check that the particles are currently
    //approaching and overlapping
    BOOST_CHECK(eval(f, t==0.0) <= 0);
    BOOST_CHECK(eval(df, t==0.0) < 0);
  } else {
    //There is a root. First determine what the next root actually is
    double nextroot = HUGE_VAL;
    for (double root : roots)
      if (root > 0)
	nextroot = std::min(nextroot, root);

    double nextdroot = HUGE_VAL;
    for (double root : droots)
      if (root > 0)
	nextdroot = std::min(nextdroot, root);

    if (eval(f, t==0) >= 0) {
      //Particles started out not overlapping, therefore the solution
      //must be an actual root
      if (nextroot == HUGE_VAL) {
	//Check if this is a phantom root
	BOOST_CHECK(nextroot == HUGE_VAL);
      } else {
	if (!::boost::test_tools::check_is_close(solution, nextroot, ::boost::test_tools::percent_tolerance(tol))) {
	  std::cout.precision(30);
	  std::cout << magnet::intersection::nextEvent(f) << std::endl;
	  std::cout << "f(x)=" << f << std::endl;
	  std::cout << "f'(x)=" << df << std::endl;
	  std::cout << "f''(x)=" << derivative(df, t) << std::endl;
	  std::cout << "f("<< solution <<")=" << eval(f, t == solution) << std::endl;
	  std::cout << "f("<< nextroot <<")=" << eval(f, t == nextroot) << std::endl;
	  std::cout << "f'("<< nextroot <<")=" << eval(df, t == nextroot) << std::endl;
	  std::cout << "roots = " << roots << std::endl;
	  std::cout << "f' roots = " << droots << std::endl;
	  BOOST_REQUIRE_CLOSE(solution, nextroot, tol);
	}
      }
    } else /*if (f(0) <= 0)*/ {
      //The particles started out overlapping, 
      if (::boost::test_tools::check_is_close(solution, nextdroot, ::boost::test_tools::percent_tolerance(tol))) {
	//Its at the next turning point, the root must be overlapping
	if (eval(f, t==solution) > 4 * precision(f, solution)) {
	  std::cout.precision(30);
	  std::cout << magnet::intersection::nextEvent(f) << std::endl;
	  std::cout << "f(x)=" << f << std::endl;
	  std::cout << "f'(x)=" << df << std::endl;
	  std::cout << "f''(x)=" << derivative(df, t) << std::endl;
	  std::cout << "f(0)=" << eval(f, t==0) << std::endl;
	  std::cout << "f'(0)=" << eval(df, t==0) << std::endl;
	  std::cout << "f("<< solution <<")=" << eval(f, t==solution) << std::endl;
	  std::cout << "f("<< nextroot <<")=" << eval(f, t==nextroot) << std::endl;
	  std::cout << "f'("<< nextroot <<")=" << eval(df, t==nextroot) << std::endl;
	  std::cout << "roots = " << roots << std::endl;
	  std::cout << "f' roots = " << droots << std::endl;
	  BOOST_ERROR("Turning point event is not within the overlapped zone");
	}
	if (std::abs(eval(df, t==solution)) > 4 * precision(f, solution)){
	  std::cout.precision(30);
	  std::cout << magnet::intersection::nextEvent(f) << std::endl;
	  std::cout << "f(x)=" << f << std::endl;
	  std::cout << "f'(x)=" << df << std::endl;
	  std::cout << "f''(x)=" << derivative(df, t) << std::endl;
	  std::cout << "f(0)=" << eval(f, t==0) << std::endl;
	  std::cout << "f'(0)=" << eval(df, t==0) << std::endl;
	  std::cout << "f("<< solution <<")=" << eval(f, t==solution) << std::endl;
	  std::cout << "f'("<< solution <<")=" << eval(df, t==solution) << std::endl;
	  std::cout << "f("<< nextroot <<")=" << eval(f, t==nextroot) << std::endl;
	  std::cout << "f'("<< nextroot <<")=" << eval(df, t==nextroot) << std::endl;
	  std::cout << "roots = " << roots << std::endl;
	  std::cout << "f' roots = " << droots << std::endl;
	  BOOST_ERROR("Particles are approaching at turning point root!");
	}
      } else {
	//The detected root must be after the nextroot (which is the
	//exit root)
	BOOST_CHECK(solution >= nextroot);
	if (std::abs(eval(f, t==solution)) > std::abs(tol * solution)) {
	  std::cout.precision(30);
	  std::cout << magnet::intersection::nextEvent(f) << std::endl;
	  std::cout << "f(x)=" << f << std::endl;
	  std::cout << "f'(x)=" << df << std::endl;
	  std::cout << "f''(x)=" << derivative(df, t) << std::endl;
	  std::cout << "f(0)=" << eval(f, t==0) << std::endl;
	  std::cout << "f'(0)=" << eval(df, t==0) << std::endl;
	  std::cout << "f("<< solution <<")=" << eval(f, t==solution) << std::endl;
	  std::cout << "f("<< nextroot <<")=" << eval(f, t==nextroot) << std::endl;
	  std::cout << "f'("<< nextroot <<")=" << eval(df, t==nextroot) << std::endl;
	  std::cout << "roots = " << roots << std::endl;
	  std::cout << "f' roots = " << droots << std::endl;
	  BOOST_ERROR("This is not a root!");
	}
	double err = HUGE_VAL;
	for (double root : roots)
	  err = std::min(err, std::abs(root - solution));
	BOOST_CHECK_SMALL(err, tol);
      }
    }
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
	auto s_poly = shift_function(poly, shift_dist(RNG));
	auto roots = solve_roots(s_poly);
	test_solution(s_poly, magnet::intersection::nextEvent(s_poly), 1e-10);
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
	  test_solution(s_poly, magnet::intersection::nextEvent(s_poly), 1e-8);
	}
      }
}

BOOST_AUTO_TEST_CASE( Cubic_function )
{
  RNG.seed(1);
  
const double rootvals[] = {-1e7, -1e3, -3.14159265, -1, 0, 1, 3.14159265, +100, 1e3, 1e7 };
  double sign = -1;
  double root1 = -1e7;
  double root2 = 3.14159265;
  double root3 = 3.14159265;
//  for (double sign : {-1.0, +1.0})
//    for (double root1 : rootvals)
//      for (double root2 : rootvals) 
//	for (double root3 : rootvals) 
	  {
	    auto poly = (t - root1) * (t - root2) * (t - root3) * sign;
	    std::uniform_real_distribution<double> shift_dist(-10, 10);
	    for (size_t i(0); i < tests; ++i) {
	      double shift = shift_dist(RNG);
	      auto s_poly = shift_function(poly, shift);
	      std::cout << "roots = " << root1-shift << ", " << root2-shift << ", " << root3-shift << std::endl;
	      std::cout << "shift = " << shift << std::endl;
	      std::cout << "Next_event = " << magnet::intersection::nextEvent(s_poly) << std::endl;
	      std::cout << "f(0) = " << eval(s_poly, Variable<'t'>()==0.0) << std::endl;
	      std::cout << "f'(0) = " << eval(derivative(s_poly, Variable<'t'>()), Variable<'t'>()==0.0) << std::endl;
	      test_solution(s_poly, magnet::intersection::nextEvent(s_poly), 1e-4);
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
//		test_solution(s_poly, magnet::intersection::nextEvent(s_poly, 1.0), 1e-4);
//	      }
//	    }
//}
