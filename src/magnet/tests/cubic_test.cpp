#define BOOST_TEST_MODULE Cubic_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <magnet/math/quartic.hpp>
#include <vector>
#include "quartic_original.hpp"
#include <complex>

double rootvals[] = {1e8, -1e6, -1e3, -100, -1, 0, 1, +100, 1e3, 1e6, 1e8 };
BOOST_AUTO_TEST_SUITE( cubic_root_solver )

BOOST_AUTO_TEST_CASE( triple_roots )
{
  for (double root1 : rootvals)
    for (double root2 : rootvals)
      for (double root3 : rootvals)
	{
	  double a = - root1 - root2  - root3,
	    b = root1 * root2 + root1 * root3 + root2 * root3,
	    c = - root1 * root2 * root3;
	  
	  //Don't test the case where there is only one root (x^3=c)
	  if ((a == 0) && (b == 0)) continue;

	  double roots[3] = {0, 0, 0};
	  size_t rootcount = magnet::math::cubicSolve(a, b, c, roots[0], roots[1], roots[2]);
	  double actual_roots[] = {root1, root2, root3};
	  std::sort(actual_roots, actual_roots + 3);
	  std::sort(roots, roots + rootcount);

	  BOOST_CHECK_EQUAL(rootcount, 3);

	  if (rootcount == 3)
	    for (size_t i = 0; i < 3; ++i)
	      {
		double root_error = std::abs((roots[i] - actual_roots[i]) / (actual_roots[i] + (actual_roots[i] == 0)));
		BOOST_CHECK_MESSAGE(root_error < 0.001, "root_error=" << root_error << " [a,b,c]=[" << a << "," << b << "," << c << "] roots=[" << roots[0] << "," << roots[1] << "," << roots[2] << "] actual_roots=[" << root1 << "," << root2 << "," << root3 << "]");
	      }
	}
}
  
BOOST_AUTO_TEST_CASE( single_roots )
{
  for (double root1 : rootvals)
    for (double root2real : rootvals)
      for (double root2im : rootvals)
	{
	  //Skip three real root cases
	  if (root2im <= 0) continue;

	  std::complex<double>
	    root1val(root1, 0),
	    root2val(root2real, root2im),
	    root3val(root2real, -root2im);
	  
	  double a = (-root1val - root2val  - root3val).real(),
	    b = (root1val * root2val + root1val * root3val + root2val * root3val).real(),
	    c = - (root1val * root2val * root3val).real();
	  
	  double roots[] = {0,0,0};
	  size_t rootcount = magnet::math::cubicSolve(a, b, c, roots[0], roots[1], roots[2]);

	  BOOST_CHECK_MESSAGE(rootcount == 1, "rootcount=" << rootcount << " [a,b,c]=[" << a << "," << b << "," << c << "] roots=[" << roots[0] << "," << roots[1] << "," << roots[2] << "] actual_roots=[" << root1 << "," << root2real << " +- " << root2im << "i]");

	  if (rootcount == 1)
	    {
	      double root_error = std::abs((roots[0] - root1) / (root1 + (root1 == 0)));
	      BOOST_CHECK_MESSAGE(root_error < 0.001, "root error is " << root_error);
	    }
	}
}

BOOST_AUTO_TEST_SUITE_END()
