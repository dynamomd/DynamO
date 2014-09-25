#define BOOST_TEST_MODULE Cubic_Quartic_test
#include <boost/test/unit_test.hpp>
#include <magnet/math/quartic.hpp>
#include <complex>

double cubic_rootvals[] = {-1e7, -1e6, -1e3, -100, -1, 0, 1, +100, 1e3, 1e6, 1e7 };

BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES( cubic_triple_roots, 0)
BOOST_AUTO_TEST_CASE( cubic_triple_roots )
{
  for (double root1 : cubic_rootvals)
    for (double root2 : cubic_rootvals)
      if (root2 != root1)
	for (double root3 : cubic_rootvals)
	  if ((root3 != root2) && (root3 != root1))
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

BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES( cubic_single_roots, 0 )
BOOST_AUTO_TEST_CASE( cubic_single_roots )
{
  for (double root1 : cubic_rootvals)
    for (double root2real : cubic_rootvals)
      for (double root2im : cubic_rootvals)
	{
	  //Skip real second root cases, and the symmetric cases
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


double quartic_rootvals[] = {-1e3, -100, -10, -1, 0, 1, 10, +100, 1e3};

BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES( quartic_four_roots, 0 )
BOOST_AUTO_TEST_CASE( quartic_four_roots )
{
  for (double root1 : quartic_rootvals)
    for (double root2 : quartic_rootvals)
      if (root1 != root2)
	for (double root3 : quartic_rootvals)
	  if ((root3 != root2) && (root3 != root1))
	    for (double root4 : quartic_rootvals)
	      if ((root4 != root3) && (root4 != root2) && (root4 != root1))
		{
		  double a = - root1 - root2 - root3 - root4,
		    b = root1 * root2 + root1 * root3 + root1 * root4 + root2 * root3 + root2 * root4 + root3 * root4,
		    c = -root1 * root2 * root3 - root1 * root2 * root4 - root1 * root3 * root4 - root2 * root3 * root4,
		    d = root1 * root2 * root3 * root4;
	    
		  double actual_roots[] = {root1, root2, root3, root4};
		  double roots[4] = {0, 0, 0, 0};
		  size_t rootcount = magnet::math::quarticSolve(a, b, c, d, roots[0], roots[1], roots[2], roots[3]);
		  //size_t rootcount = quartic(a, b, c, d, roots);
		  std::sort(actual_roots, actual_roots + 4);
		  std::sort(roots, roots + rootcount);

		  BOOST_CHECK_MESSAGE(rootcount == 4, "rootcount=" << rootcount << " [a,b,c,d]=[" << a << "," << b << "," << c << "," << d << "] roots=[" << roots[0] << "," << roots[1] << "," << roots[2] << "," << roots[3] << "] actual_roots=[" << root1 << "," << root2 << "," << root3 << "," << root4 << "]");
	    
		  if (rootcount == 4)
		    for (size_t i = 0; i < 4; ++i)
		      {
			double root_error = std::abs((roots[i] - actual_roots[i]) / (actual_roots[i] + (actual_roots[i] == 0)));
			BOOST_CHECK_MESSAGE(root_error < 0.0002, "root_error=" << root_error << " [a,b,c,d]=[" << a << "," << b << "," << c << "," << d << "] roots=[" << roots[0] << "," << roots[1] << "," << roots[2] << "," << roots[3] << "] actual_roots=[" << root1 << "," << root2 << "," << root3 << "," << root4 << "]");
		      }
		}
}

BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES( quartic_two_roots, 0 )
BOOST_AUTO_TEST_CASE( quartic_two_roots )
{
  for (double root1 : quartic_rootvals)
    for (double root2 : quartic_rootvals)
      if (root2 != root1)
	for (double root3real : quartic_rootvals)
	  for (double root3im : quartic_rootvals)
	    if (root3im > 0)
	      {
		std::complex<double> root1val(root1, 0), root2val(root2, 0), root3val(root3real, root3im), root4val(root3real, -root3im);
		
		double a = (- root1val - root2val - root3val - root4val).real(),
		  b = (+ root1val * root2val + root1val * root3val + root1val * root4val + root2val * root3val + root2val * root4val + root3val * root4val).real(),
		  c = (- root1val * root2val * root3val - root1val * root2val * root4val - root1val * root3val * root4val - root2val * root3val * root4val).real(),
		  d = (root1val * root2val * root3val * root4val).real();

		double actual_roots[] = {root1, root2};
		double roots[4] = {0, 0, 0, 0};
		size_t rootcount = magnet::math::quarticSolve(a, b, c, d, roots[0], roots[1], roots[2], roots[3]);
		std::sort(actual_roots, actual_roots + 2);
		std::sort(roots, roots + rootcount);
		
		BOOST_CHECK_MESSAGE(rootcount == 2, "rootcount=" << rootcount << " [a,b,c,d]=[" << a << "," << b << "," << c << "," << d << "] roots=[" << roots[0] << "," << roots[1] << "," << roots[2] << "," << roots[3] << "]" << " actual_roots=[" << root1 << "," << root2 << "," << root3real << "+-" << root3im << "i]");
		
		if (rootcount == 2)
		  for (size_t i = 0; i < 2; ++i)
		    {
		      double root_error = std::abs((roots[i] - actual_roots[i]) / (actual_roots[i] + (actual_roots[i] == 0)));
		      BOOST_CHECK_MESSAGE(root_error < 0.0002, "root_error=" << root_error << " [a,b,c,d]=[" << a << "," << b << "," << c << "," << d << "] roots=[" << roots[0] << "," << roots[1] << "," << roots[2] << "," << roots[3] << "]" << " actual_roots=[" << root1 << "," << root2 << "," << root3real << "+-" << root3im << "i]");
		    }
	      }
}

BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES( quartic_no_roots, 0 )
BOOST_AUTO_TEST_CASE( quartic_no_roots)
{
  for (double root1real : quartic_rootvals)
    for (double root1im : quartic_rootvals)
      for (double root2real : quartic_rootvals)
	for (double root2im : quartic_rootvals)
	  {
	    if (root1im <= 0) continue;
	    if (root2im <= 0) continue;
	    if ((root1real == root2real) && (root1im == root2im)) continue;
	    
	    std::complex<double> root1val(root1real, root1im), root2val(root1real, -root1im), root3val(root2real, root2im), root4val(root2real, -root2im);
		
	    double a = (- root1val - root2val - root3val - root4val).real(),
	      b = (+ root1val * root2val + root1val * root3val + root1val * root4val + root2val * root3val + root2val * root4val + root3val * root4val).real(),
	      c = (- root1val * root2val * root3val - root1val * root2val * root4val - root1val * root3val * root4val - root2val * root3val * root4val).real(),
	      d = (root1val * root2val * root3val * root4val).real();

	    double roots[4] = {0, 0, 0, 0};
	    size_t rootcount = magnet::math::quarticSolve(a, b, c, d, roots[0], roots[1], roots[2], roots[3]);
	    std::sort(roots, roots + rootcount);
		
	    BOOST_CHECK_MESSAGE(rootcount == 0, "rootcount=" << rootcount << " [a,b,c,d]=[" << a << "," << b << "," << c << "," << d << "] roots=[" << roots[0] << "," << roots[1] << "," << roots[2] << "," << roots[3] << "]" << " actual_roots=[" << root1real << ",+-" << root1im << "i," << root2real << "+-" << root2im << "i]");
	  }
}
