#define BOOST_TEST_MODULE Polynomial_test
#include <boost/test/included/unit_test.hpp>
#include <vector>
#include <magnet/math/symbolic.hpp>
#include <magnet/math/vector.hpp>
#include <cmath>
#include <complex>

bool err(double val, double expected)
{
  return std::abs(val / expected - 1) < 0.0001;
}

template<class T1, class T2>
bool compare_expression(const T1& f, const T2& g) {
  using namespace magnet::math;
  std::ostringstream os;
  os << f;
  std::string f_str = os.str();
  os.str(""); os.clear();
  os << g;
  std::string g_str = os.str();
  if (!(f_str == g_str))
    std::cerr << f << " != " << g << std::endl;
  return f_str == g_str;
}

BOOST_AUTO_TEST_CASE( poly_variables )
{
  using namespace magnet::math;
  Polynomial<1> x{0, 1};
  Polynomial<1,double,'y'> y{0, 1};

  BOOST_CHECK(compare_expression(x * x * x, "x³"));
  BOOST_CHECK(compare_expression(y * y * y, "y³"));
  BOOST_CHECK(compare_expression(substitution(y * y * y, Variable<'y'>()==Variable<'x'>()), "x³"));
}

BOOST_AUTO_TEST_CASE( poly_addition )
{
  using namespace magnet::math;
  Polynomial<1> x{0, 2.5};
  Polynomial<0> C{0.3};
  auto poly1 = x+C;
  BOOST_CHECK_EQUAL(poly1[0], 0.3);
  BOOST_CHECK_EQUAL(poly1[1], 2.5);

  auto poly2 = x + 0.3;
  BOOST_CHECK_EQUAL(poly2[0], 0.3);
  BOOST_CHECK_EQUAL(poly2[1], 2.5);
}

BOOST_AUTO_TEST_CASE( poly_multiplication )
{
  using namespace magnet::math;
  Polynomial<1> x{0, 1};
  auto poly1 = -2.0;
  auto poly2 = 2.0 - x + x * x;
  auto poly3 = poly2 * poly1;
  BOOST_CHECK_EQUAL(poly3[0], -4);
  BOOST_CHECK_EQUAL(poly3[1], +2);
  BOOST_CHECK_EQUAL(poly3[2], -2);
}

BOOST_AUTO_TEST_CASE( poly_division )
{
  using namespace magnet::math;
  Polynomial<1> x{0, 1};
  auto poly1 = 2.0 - x + x * x;
  auto poly2 = poly1 / 0.5;
  BOOST_CHECK_EQUAL(poly2[0], 4);
  BOOST_CHECK_EQUAL(poly2[1], -2);
  BOOST_CHECK_EQUAL(poly2[2], 2);
}

BOOST_AUTO_TEST_CASE( poly_vector )
{
  using namespace magnet::math;
  Polynomial<1, Vector> x{Vector(), Vector{1,2,3}};
  Polynomial<0, Vector> C{Vector{3,2,1}};
  auto poly1 = x+C;
  BOOST_CHECK_EQUAL(poly1[0], (Vector{3,2,1}));
  BOOST_CHECK_EQUAL(poly1[1], (Vector{1,2,3}));
  
  auto poly2 = poly1 * poly1;
  BOOST_CHECK_EQUAL(poly2[0], 14);
  BOOST_CHECK_EQUAL(poly2[1], 20);
  BOOST_CHECK_EQUAL(poly2[2], 14);
}

BOOST_AUTO_TEST_CASE( poly_lower_order )
{
  using namespace magnet::math;
  Polynomial<1> x{0, 1};
  Polynomial<2> poly2 = 2.0 - x + x * x;
  Polynomial<3> poly3 = poly2 + 0 * x * x * x;
  //Try to cast down one level as the highest order coefficient is zero
  BOOST_CHECK(poly3[3] == 0);
  Polynomial<2> poly4(change_order<2>(poly3));

  BOOST_CHECK_EQUAL(poly4[0], 2);
  BOOST_CHECK_EQUAL(poly4[1], -1);
  BOOST_CHECK_EQUAL(poly4[2], 1);
  BOOST_CHECK_EQUAL(eval(poly3, 123), eval(poly4, 123));
}

BOOST_AUTO_TEST_CASE( poly_eval_limits )
{
  using namespace magnet::math;
  Polynomial<1> x{0, 1};
  
  {//Check even positive polynomials
    auto f = x*x-x+3;
    BOOST_CHECK_EQUAL(eval(f, 0), 3);
    BOOST_CHECK_EQUAL(eval(f, +HUGE_VAL), +HUGE_VAL);
    BOOST_CHECK_EQUAL(eval(f, -HUGE_VAL), +HUGE_VAL);
  }

  {//Check even negative polynomials
    auto f = -x*x+x+3;
    BOOST_CHECK_EQUAL(eval(f, 0), 3);
    BOOST_CHECK_EQUAL(eval(f, +HUGE_VAL), -HUGE_VAL);
    BOOST_CHECK_EQUAL(eval(f, -HUGE_VAL), -HUGE_VAL);
  }

  {//Check odd positive polynomials
    auto f = x*x*x + x + 3;
    BOOST_CHECK_EQUAL(eval(f, +HUGE_VAL), +HUGE_VAL);
    BOOST_CHECK_EQUAL(eval(f, -HUGE_VAL), -HUGE_VAL);
  }

  {//Check odd negative polynomials
    auto f = -x*x*x + x + 3;
    BOOST_CHECK_EQUAL(eval(f, +HUGE_VAL), -HUGE_VAL);
    BOOST_CHECK_EQUAL(eval(f, -HUGE_VAL), +HUGE_VAL);
  }
}

BOOST_AUTO_TEST_CASE( poly_derivative )
{
  using namespace magnet::math;
  Polynomial<1> x{0, 1};
  auto poly1 = x + x*x + x*x*x + x*x*x*x;
  auto poly2 = derivative(poly1, Variable<'x'>());
  BOOST_CHECK_EQUAL(poly2[0], 1);
  BOOST_CHECK_EQUAL(poly2[1], 2);
  BOOST_CHECK_EQUAL(poly2[2], 3);  
  BOOST_CHECK_EQUAL(poly2[3], 4);  

  BOOST_CHECK_CLOSE(eval(poly2, 3.14159), eval_derivatives<1>(poly1, 3.14159)[1], 1e-10);

  auto poly3 = 2.0 - x + 2 * x * x;
  auto poly4 = derivative(poly3, Variable<'x'>());
  BOOST_CHECK_EQUAL(poly4[0], -1);
  BOOST_CHECK_EQUAL(poly4[1], 4);
  BOOST_CHECK_EQUAL(eval(poly4, 0), -1);
  BOOST_CHECK_EQUAL(eval(poly4, 1), 3);
}

BOOST_AUTO_TEST_CASE( poly_zero_derivative)
{
  using namespace magnet::math;
  const Polynomial<1> x{0, 1};
  const auto poly1 = derivative(x, Variable<'x'>());
  BOOST_CHECK_EQUAL(poly1[0], 1);

  const auto poly2 = derivative(poly1, Variable<'x'>());
  BOOST_CHECK(compare_expression(poly2, NullSymbol()));
}

BOOST_AUTO_TEST_CASE( poly_deflation)
{
  using namespace magnet::math;
  Polynomial<1> x{0, 1};

  const double roots[] = {-1e3, 4e3, 0, 3.14159265, -3.14159265};
  for (double root1: roots)
    for (double root2: roots)
      for (double root3: roots) {
	auto poly = (x - root1) * (x - root2) * (x-root3);
	
	auto deflated = deflate_polynomial(poly, root1);
	auto exact = (x - root2) * (x-root3);
	for (size_t i(0); i < 3; ++i)
	  if (exact[i] != 0)
	    BOOST_CHECK_CLOSE(deflated[i], exact[i], 1e-10);
	  else
	    BOOST_CHECK_SMALL(deflated[i], 1e-10);
	
	deflated = deflate_polynomial(poly, root2);
	exact = (x - root1) * (x-root3);
	for (size_t i(0); i < 3; ++i)
	  if (exact[i] != 0)
	    BOOST_CHECK_CLOSE(deflated[i], exact[i], 1e-10);
	  else
	    BOOST_CHECK_SMALL(deflated[i], 1e-10);

	deflated = deflate_polynomial(poly, root3);
	exact = (x - root1) * (x-root2);
	for (size_t i(0); i < 3; ++i)
	  if (exact[i] != 0)
	    BOOST_CHECK_CLOSE(deflated[i], exact[i], 1e-10);
	  else
	    BOOST_CHECK_SMALL(deflated[i], 1e-10);
      }
}

BOOST_AUTO_TEST_CASE( poly_shift)
{
  using namespace magnet::math;
  Polynomial<1> x{0, 1};

  const double roots[] = {-1e3, 4e3, 0, 3.14159265, -3.14159265};
  for (double root1: roots)
    for (double root2: roots)
      for (double root3: roots) {
	auto f = (x - root1) * (x - root2) * (x-root3);

	for (double shift : {-1.0, 2.0, 1e3, 3.14159265, -1e5}) {
	  auto g = shift_function(f, shift);
	  
	  BOOST_CHECK_CLOSE(eval(g, 0), eval(f,shift), 1e-10);
	  BOOST_CHECK_CLOSE(eval(g, 1e3), eval(f, 1e3 + shift), 1e-10);
	}
      }
}

BOOST_AUTO_TEST_CASE( poly_quadratic_roots)
{
  using namespace magnet::math;
  Polynomial<1> x{0, 1};
  
  {//Quadratic with no roots
    auto poly = x * x - 3 * x + 4;
    auto roots = solve_real_roots(poly);
    BOOST_CHECK(roots.size() == 0);
    BOOST_CHECK_EQUAL(next_root(poly), HUGE_VAL);
  }
  
  {//Quadratic with one root
    auto poly = -4 * x * x + 12 * x - 9;
    auto roots = solve_real_roots(poly);
    BOOST_CHECK(roots.size() == 1);
    if (roots.size() == 1)
      BOOST_CHECK_CLOSE(roots[0], 1.5, 1e-10);

    BOOST_CHECK_CLOSE(next_root(poly), 1.5, 1e-10);
  }
  
  {//quadratic but linear function with one root
    auto poly =  0 * x * x + 12 * x - 9;
    auto roots = solve_real_roots(poly);
    BOOST_CHECK(roots.size() == 1);
    if (roots.size() == 1)
      BOOST_CHECK_CLOSE(roots[0], 0.75, 1e-10);
    BOOST_CHECK_CLOSE(next_root(poly), 0.75, 1e-10);
  }

  {//constant function, with no roots
    auto poly =  0 * x * x + 0 * x - 9;
    auto roots = solve_real_roots(poly);
    BOOST_CHECK(roots.size() == 0);
    BOOST_CHECK_EQUAL(next_root(poly), HUGE_VAL);
  }
}

BOOST_AUTO_TEST_CASE( poly_quadratic_special_cases)
{
  using namespace magnet::math;
  Polynomial<1> x{0, 1};
  
  {//Quadratic with catastrophic cancellation of error
    auto poly = x * x + 712345.12 * x + 1.25;  
    auto roots = solve_real_roots(poly);
    BOOST_CHECK(roots.size() == 2);
    BOOST_CHECK_CLOSE(roots[0], -712345.1199985961, 1e-10);
    BOOST_CHECK_CLOSE(roots[1], -1.754767408250742e-6, 1e-10);
    BOOST_CHECK_EQUAL(next_root(poly), HUGE_VAL);
  }

  const double maxsqrt = std::sqrt(std::numeric_limits<double>::max());
  const double largeterm = maxsqrt * 100;
  {//Large linear coefficient
    auto poly = x * x + largeterm * x + 1.25;
    auto roots = solve_real_roots(poly);
    BOOST_CHECK(roots.size() == 2);
    //Mathematica value
    BOOST_CHECK_CLOSE(roots[0], -1.3407807929942599e156, 1e-10);
    BOOST_CHECK_CLOSE(roots[1], -9.322925914000258e-157, 1e-10);
    BOOST_CHECK_EQUAL(next_root(poly), HUGE_VAL);
  }

  {//Large (+ve) constant coefficient
    auto poly = x * x + x + largeterm;
    auto roots = solve_real_roots(poly);
    BOOST_CHECK(roots.size() == 0);
    BOOST_CHECK_EQUAL(next_root(poly), HUGE_VAL);
  }
  {//Large (-ve) constant coefficient
    auto poly = x * x + x - largeterm;
    auto roots = solve_real_roots(poly);
    BOOST_CHECK(roots.size() == 2);

    //Mathematica value
    BOOST_CHECK_CLOSE(roots[0], -1.157920892373162e78, 1e-10);
    BOOST_CHECK_CLOSE(roots[1], 1.157920892373162e78, 1e-10);
    BOOST_CHECK_CLOSE(next_root(poly), 1.157920892373162e78, 1e-10);
  }
}

double cubic_rootvals[] = {-1e6, -1e3, -100, -1, 0, 1, +100, 1e3, 1e6};

BOOST_AUTO_TEST_CASE( poly_cubic_triple_roots )
{
  using namespace magnet::math;
  const Polynomial<1> x{0, 1};

  for (double root1 : cubic_rootvals)
    for (double root2 : cubic_rootvals)
      if (root2 != root1)
	for (double root3 : cubic_rootvals)
	  if ((root1 != root2) && (root2 != root3) && (root1 != root3))
	    {
	      auto f = (x - root1) * (x - root2) * (x - root3);
	      //Don't test the case where there is only one root (x^3=c)
	      if ((f[2] == 0) && (f[1] == 0)) continue;

	      auto roots = solve_real_roots(f);
	      decltype(roots) actual_roots = {root1, root2, root3};
	      std::sort(actual_roots.begin(), actual_roots.end());

	      BOOST_CHECK_MESSAGE(roots.size() == 3, f << " roots=[" << roots[0] << "," << roots[1] << "," << roots[2] << "] actual_roots=[" << root1 << "," << root2 << "," << root3 << "]");

	      if (roots.size() == 3) {
		for (size_t i = 0; i < 3; ++i)
		  {
		    double root_error = std::abs((roots[i] - actual_roots[i]) / (actual_roots[i] + (actual_roots[i] == 0)));
		    BOOST_CHECK_MESSAGE(root_error < 0.001, "root_error=" << root_error << " " << f << " roots=[" << roots[0] << "," << roots[1] << "," << roots[2] << "] actual_roots=[" << root1 << "," << root2 << "," << root3 << "]");
		  }
		
		//Test the next_root implementation
		std::vector<double> pos_roots;
		for (double root : actual_roots)
		  if (root >= 0)
		    pos_roots.push_back(root);
		std::sort(pos_roots.begin(), pos_roots.end());
		if (pos_roots.size() != 0)
		  BOOST_CHECK_CLOSE(pos_roots[0], next_root(f), 1e-11);
		else
		  BOOST_CHECK_EQUAL(next_root(f), HUGE_VAL);
	      }
	    }
}

BOOST_AUTO_TEST_CASE( poly_cubic_single_roots )
{
  using namespace magnet::math;
  const Polynomial<1> x{0, 1};

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
	  
	  auto poly_c = (x - root1val) * (x - root2val) * (x - root3val);
	  
	  Polynomial<3, double> f = poly_c[0].real() + poly_c[1].real() * x + poly_c[2].real() * x * x + poly_c[3].real() * x * x * x;

	  auto roots = solve_real_roots(f);
	  BOOST_CHECK_MESSAGE(roots.size() == 1, "rootcount=" << roots.size() << " " << f << " roots=[" << roots[0] << "," << roots[1] << "," << roots[2] << "] actual_roots=[" << root1 << "," << root2real << " +- " << root2im << "i]");

	  if (roots.size() == 1)
	    {
	      double root_error = std::abs((roots[0] - root1) / (root1 + (root1 == 0)));
	      BOOST_CHECK_MESSAGE(root_error < 0.001, "root error is " << root_error);
	      if (roots[0] >= 0)
		BOOST_CHECK_CLOSE(roots[0], next_root(f), 1e-11);
	    }
	}
}

BOOST_AUTO_TEST_CASE( poly_cubic_special_cases )
{
  using namespace magnet::math;
  const Polynomial<1> x{0, 1};

  {//Zero constant term with three roots
    auto poly = (x * x + 712345.12 * x + 1.25) * x;
    auto roots = solve_real_roots(poly);
    BOOST_CHECK(roots.size() == 3);
    BOOST_CHECK_CLOSE(roots[0], -712345.1199985961, 1e-10);
    BOOST_CHECK_CLOSE(roots[1], -1.754767408250742e-6, 1e-10);
    BOOST_CHECK_CLOSE(roots[2], 0, 1e-10);
    BOOST_CHECK_CLOSE(next_root(poly), 0, 1e-11);
  }

  {//Zero constant term with one root
    auto poly = (x * x - 3 * x + 4) * x;
    auto roots = solve_real_roots(poly);
    BOOST_CHECK(roots.size() == 1);
    BOOST_CHECK_CLOSE(roots[0], 0, 1e-10);
    BOOST_CHECK_CLOSE(next_root(poly), 0, 1e-11);
  }

  {//Special case where f(x) = a * x^3 + d
    auto poly = x * x * x + 1e3;
    auto roots = solve_real_roots(poly);
    BOOST_CHECK(roots.size() == 1);
    BOOST_CHECK_CLOSE(roots[0], -10, 1e-10);

    poly = x * x * x - 1e3;
    roots = solve_real_roots(poly);
    BOOST_CHECK(roots.size() == 1);
    BOOST_CHECK_CLOSE(roots[0], 10, 1e-10);
    BOOST_CHECK_CLOSE(next_root(poly), 10, 1e-11);
  }

  const double maxsqrt = std::sqrt(std::numeric_limits<double>::max());
  const double largeterm = maxsqrt * 100;

  {//Large x^2 term
    auto poly = x * x * x - largeterm * x * x + 1.25;
    auto roots = solve_real_roots(poly);
    BOOST_CHECK_EQUAL(roots.size(), 3);
    BOOST_CHECK_CLOSE(roots[0], -9.655529977168658e-79, 1e-10);
    BOOST_CHECK_CLOSE(roots[1], +9.655529977168658e-79, 1e-10);
    BOOST_CHECK_CLOSE(roots[2], 1.3407807929942599e156, 1e-10);
    BOOST_CHECK_CLOSE(next_root(poly), +9.655529977168658e-79, 1e-11);
  }

  {//Large x term
    auto poly = x * x * x - x * x - largeterm * x + 1.25;
    auto roots = solve_real_roots(poly);
    BOOST_CHECK_EQUAL(roots.size(), 3);
    BOOST_CHECK_CLOSE(roots[0], -1.1579208923731622e78, 1e-10);
    BOOST_CHECK_CLOSE(roots[1], 9.322925914000258e-157, 1e-10);
    BOOST_CHECK_CLOSE(roots[2], 1.1579208923731622e78, 1e-10);
    BOOST_CHECK_CLOSE(next_root(poly), 9.322925914000258e-157, 1e-11);
  }

  const double smallerterm = maxsqrt * 1e-1;
  {
    //Large v term
    auto poly = x * x * x  - smallerterm * x * x - smallerterm * x + 2;
    auto roots = solve_real_roots(poly);
    BOOST_CHECK_EQUAL(roots.size(), 3);
    BOOST_CHECK_CLOSE(roots[0], -1.0, 1e-10);
    BOOST_CHECK_CLOSE(roots[1], 1.491668146240041472864517142264024641421371730393e-153, 1e-10);
    BOOST_CHECK_CLOSE(roots[2], 1.340780792994259598314974448015366224371799690462e153, 1e-10);
    BOOST_CHECK_CLOSE(next_root(poly), 1.491668146240041472864517142264024641421371730393e-153, 1e-11);
  }
}

BOOST_AUTO_TEST_CASE( poly_root_tests)
{
  using namespace magnet::math;
  const Polynomial<1> x{0, 1};

  {//Check cubic
    auto f1 = 4 * (x*x*x) - x * x - 2*x +12;
    
    auto roots = solve_real_roots(f1);
    BOOST_CHECK(roots.size() == 1);
    BOOST_CHECK_CLOSE(roots[0], -1.472711896724616002268033950475380144341, 1e-10);
    
    auto droots = solve_real_roots(derivative(f1, Variable<'x'>()));
    BOOST_CHECK(droots.size() == 2);
    BOOST_CHECK_CLOSE(droots[0], -1.0/3, 1e-10);
    BOOST_CHECK_CLOSE(droots[1], 0.5, 1e-10);
  }

  {//Check quartic
    auto f1 = 10 * (x*x*x*x) + x*x*x - 30 * x * x -23;

    auto roots = solve_real_roots<>(f1);
    BOOST_CHECK_EQUAL(roots.size(), 2);
    BOOST_CHECK_CLOSE(roots[0], -1.949403904489790210996459054473124835057, 1e-10);
    BOOST_CHECK_CLOSE(roots[1], +1.864235880634589025006445510389799368569, 1e-10);

    auto droots = solve_real_roots(derivative(f1, Variable<'x'>()));
    BOOST_CHECK_EQUAL(droots.size(),3);
    BOOST_CHECK_CLOSE(droots[0], -1.262818836058599076329128653113014315066, 1e-10);
    BOOST_CHECK_CLOSE(droots[1], 0, 1e-10);
    BOOST_CHECK_CLOSE(droots[2], +1.187818836058599076329128653113014315066, 1e-10);
  }

  {//Check PowerOp quartic
    auto f1 = simplify(pow<2>(30 * x * x + x - 23));
    
    //auto roots = solve_real_roots(f1);
    //FIX THIS UNIT TEST!
    //BOOST_CHECK(roots.size() == 2);
    ////NOTE THESE ROOTS ARE DOUBLE ROOTS (roots.size() may equal 2,3, or 4)
    //BOOST_CHECK_CLOSE(roots[0], -0.8924203103613100773375343963347855860436, 1e-7);
    //BOOST_CHECK_CLOSE(roots[1], -0.8924203103613100773375343963347855860436, 1e-7);

    auto droots = solve_real_roots(derivative(f1, Variable<'x'>()));
    BOOST_CHECK(droots.size() == 3);
    BOOST_CHECK_CLOSE(droots[0], -0.8924203103613100773375343963347855860436, 1e-10);
    BOOST_CHECK_CLOSE(droots[1], -0.01666666666666666666666666666666666666667, 1e-10);
    BOOST_CHECK_CLOSE(droots[2], +0.8590869770279767440042010630014522527103, 1e-10);
  }
}

BOOST_AUTO_TEST_CASE( poly_Euclidean_division )
{
  using namespace magnet::math;
  const Polynomial<1> x{0, 1};

  { //Standard division with remainder
    auto q = x * x * x + 3 * x - 2;
    auto g = x * x - 2 * x;
    auto r = 4 * x - 2;
    auto f = q * g + r;
    auto euclid = euclidean_division(f, g);  
    BOOST_CHECK(compare_expression(q, std::get<0>(euclid)));
    BOOST_CHECK(compare_expression(r, std::get<1>(euclid)));
  }
  {//Standard division without remainder
    auto q = x * x * x + 3 * x - 2;
    auto g = x * x - 2 * x;
    auto r = 0;
    auto f = q * g + r;
    auto euclid = euclidean_division(f, g);
    
    BOOST_CHECK(compare_expression(q, std::get<0>(euclid)));
    BOOST_CHECK(compare_expression(r, std::get<1>(euclid)));
  }

  {//Division with a zero leading order coefficient
    auto q = x * x * x + 3 * x - 2;
    auto g = 0 * x*x*x + x*x - 2 * x;
    auto r = 0;
    auto f = q * g + r;
    auto euclid = euclidean_division(f, g);
    
    BOOST_CHECK(compare_expression(q, std::get<0>(euclid)));
    BOOST_CHECK(compare_expression(r, std::get<1>(euclid)));
  }

  {//Division by a constant
    auto q = x * x * x + 3 * x - 2;
    auto g = Polynomial<0>{0.5};
    auto r = 0;
    auto f = q * g + r;
    auto euclid = euclidean_division(f, g);
    
    BOOST_CHECK(compare_expression(q, std::get<0>(euclid)));
    BOOST_CHECK(compare_expression(r, std::get<1>(euclid)));
  }

  {//Division by a high-order Polynomial which is actually a constant
    auto q = x * x * x + 3 * x - 2;
    auto g = Polynomial<3>{0.25};
    auto r = 0;
    auto f = q * g + r;
    auto euclid = euclidean_division(f, g);
    
    BOOST_CHECK(compare_expression(q, std::get<0>(euclid)));
    BOOST_CHECK(compare_expression(r, std::get<1>(euclid)));
  }
}

BOOST_AUTO_TEST_CASE( poly_Sturm_chains )
{
  using namespace magnet::math;
  const Polynomial<1> x{0, 1};

  { //Example from wikipedia (x^4+x^3-x-1)
    auto f = x*x*x*x + x*x*x - x - 1;
    auto chain = sturm_chain(f);

    BOOST_CHECK(compare_expression(chain.get(0), f));
    BOOST_CHECK(compare_expression(chain.get(1), 4*x*x*x + 3*x*x - 1));
    BOOST_CHECK(compare_expression(chain.get(2), (3.0/16) * x*x + (3.0/4)*x + (15.0/16)));
    BOOST_CHECK(compare_expression(chain.get(3), -32*x -64));
    BOOST_CHECK(compare_expression(chain.get(4), -3.0/16));
    BOOST_CHECK(compare_expression(chain.get(5), 0));
    BOOST_CHECK(compare_expression(chain.get(6), 0));
    
    //This polynomial has roots at -1 and +1
    BOOST_CHECK_EQUAL(chain.sign_changes(-HUGE_VAL), 3);
    BOOST_CHECK_EQUAL(chain.sign_changes(0), 2);
    BOOST_CHECK_EQUAL(chain.sign_changes(HUGE_VAL), 1);
    
    BOOST_CHECK_EQUAL(chain.roots(0.5, 3.0), 1);
    BOOST_CHECK_EQUAL(chain.roots(-2.141, -0.314159265), 1);
    BOOST_CHECK_EQUAL(chain.roots(-HUGE_VAL, HUGE_VAL), 2);
 }
}

BOOST_AUTO_TEST_CASE( descartes_sturm_and_budan_01_alesina_rootcount_test )
{
  using namespace magnet::math;
  const Polynomial<1> x{0, 1};
  
  //The values 0.5, and 2.0 are strong tests of the algorithms, as
  //these are the division points for VCA and VAS algorithms.
  const double roots[] = {-1e5, -0.14159265, -0.0001,0.1, 0.3333, 0.5, 0.8, 1.001, 2.0, 3.14159265, 1e7};
  
  for (const double root1: roots)
    for (const double root2: roots)
      if (root1 != root2)
	for (const double root3: roots)
	  if (!std::set<double>{root1, root2}.count(root3))
	    for (const double root4: roots)
	      if (!std::set<double>{root1, root2, root3}.count(root4))
		for (const double root5: roots)
		  if (!std::set<double>{root1, root2, root3, root4}.count(root5))
		    for (int sign : {-1, +1}) {
		      //Test where all 5 roots of a 5th order
		      //Polynomial are real
		      auto f1 = sign * (x - root1) * (x - root2) * (x - root3) * (x - root4) * (x - root5);
		      
		      //Test with the same roots, but an additional 2
		      //imaginary roots
		      auto f2 = f1 * (x * x - 3 * x + 4);

		      auto roots_in_range = [&](double a, double b) {
			return size_t
			(((root1 > a) && (root1 < b))
			+((root2 > a) && (root2 < b))
			+((root3 > a) && (root3 < b))
			+((root4 > a) && (root4 < b))
			 +((root5 > a) && (root5 < b)))
			; };

		      size_t roots_in_01 = roots_in_range(0, 1);

		      auto chain1 = sturm_chain(f1);
		      auto chain2 = sturm_chain(f2);

		      //Test interval [0,1]
		      switch (roots_in_01) {
		      case 0: 
		      case 1: 
			BOOST_CHECK_EQUAL(budan_01_test(f1), roots_in_01);
			BOOST_CHECK_EQUAL(budan_01_test(f2), roots_in_01);
			BOOST_CHECK_EQUAL(alesina_galuzzi_test(f1, 0.0, 1.0), roots_in_01);
			BOOST_CHECK_EQUAL(alesina_galuzzi_test(f2, 0.0, 1.0), roots_in_01);
			break;
		      default: 
			BOOST_CHECK(budan_01_test(f1) >= roots_in_01); 
			BOOST_CHECK(budan_01_test(f2) >= roots_in_01); 
			BOOST_CHECK(alesina_galuzzi_test(f1, 0.0, 1.0) >= roots_in_01); 
			BOOST_CHECK(alesina_galuzzi_test(f2, 0.0, 1.0) >= roots_in_01); 
			break;
		      }
		      BOOST_CHECK_EQUAL(chain1.roots(0.0, 1.0), roots_in_01); 
		      BOOST_CHECK_EQUAL(chain2.roots(0.0, 1.0), roots_in_01); 
		      
		      //Test interval [0, \infty]
		      size_t positive_roots = roots_in_range(0, HUGE_VAL);
		      switch (positive_roots) {
		      case 0: 
		      case 1:
			BOOST_CHECK_EQUAL(descartes_rule_of_signs(f1), positive_roots); 
			BOOST_CHECK_EQUAL(descartes_rule_of_signs(f2), positive_roots);
			break;
		      default: 
			BOOST_CHECK(descartes_rule_of_signs(f1) >= positive_roots); 
			break;
		      }
		      BOOST_CHECK_EQUAL(chain1.roots(0.0, HUGE_VAL), positive_roots); 
		      BOOST_CHECK_EQUAL(chain2.roots(0.0, HUGE_VAL), positive_roots); 

		      //Try some others
		      BOOST_CHECK_EQUAL(chain1.roots(-HUGE_VAL, HUGE_VAL), 5);
		      BOOST_CHECK_EQUAL(chain2.roots(-HUGE_VAL, HUGE_VAL), 5); 
		      BOOST_CHECK(alesina_galuzzi_test(f1,-1.0, 30.0) >= roots_in_range(-1, 30));
		      BOOST_CHECK(alesina_galuzzi_test(f1,-0.01, 5.0) >= roots_in_range(-0.01, 5));
		    }
}

BOOST_AUTO_TEST_CASE( LMQ_upper_bound_test )
{
  using namespace magnet::math;
  const Polynomial<1> x{0, 1};

  const double roots[] = {-1e5, -0.14159265, 3.14159265, -0.0001,0.1, 0.3333, 0.6, 1.001, 2.0, 3.14159265, 1e7};

  //Test well-formed expressions
  for (const double root1: roots)
    for (const double root2: roots)
      for (const double root3: roots)
	for (const double root4: roots)
	  for (int sign : {-1, +1})
	    {
	      auto f = sign * (x - root1) * (x - root2) * (x - root3) * (x - root4);

	      double max_root = root1;
	      max_root = std::max(max_root, root2);
	      max_root = std::max(max_root, root3);
	      max_root = std::max(max_root, root4);
	    
	      double bound = LMQ_upper_bound(f);
	      if (max_root < 0)
		BOOST_CHECK_EQUAL(bound, 0);
	      else
		BOOST_CHECK(bound >= max_root);
	    }

  //Test expressions with zero coefficients
  for (const double root1: roots)
    for (const double root2: roots)
      for (int sign : {-1, +1})
	{
	  auto f = sign * (x - root1) * (x - root2) + 0 * x*x*x*x*x;
	  double max_root = std::max(root1, root2);
	  
	  double bound = LMQ_upper_bound(f);
	  if (max_root < 0)
	    BOOST_CHECK_EQUAL(bound, 0);
	  else
	    BOOST_CHECK(bound >= max_root);
	}

  //Test constant coefficients 
  BOOST_CHECK_EQUAL(LMQ_upper_bound(1 + 0 * x*x*x*x*x), 0);
}

BOOST_AUTO_TEST_CASE( LMQ_lower_bound_test )
{
  using namespace magnet::math;
  const Polynomial<1> x{0, 1};

  const double roots[] = {-1e5, -0.14159265, 3.14159265, -0.0001,0.1, 0.3333, 0.6, 1.001, 2.0, 3.14159265, 1e7};

  for (const double root1: roots)
    for (const double root2: roots)
      for (const double root3: roots)
	for (const double root4: roots)
	  for (int sign : {-1, +1})
	    {
	      auto f = sign * (x - root1) * (x - root2) * (x - root3) * (x - root4);

	      double min_pos_root = HUGE_VAL;
	      if (root1 >= 0)
		min_pos_root = std::min(min_pos_root, root1);
	      if (root2 >= 0)
		min_pos_root = std::min(min_pos_root, root2);
	      if (root3 >= 0)
		min_pos_root = std::min(min_pos_root, root3);
	      if (root4 >= 0)
		min_pos_root = std::min(min_pos_root, root4);
	    
	      double bound = LMQ_lower_bound(f);
	      if (min_pos_root == HUGE_VAL)
		BOOST_CHECK_EQUAL(bound, HUGE_VAL);
	      else
		BOOST_CHECK(bound <= min_pos_root);
	    }

  //Test expressions with zero coefficients
  for (const double root1: roots)
    for (const double root2: roots)
      for (int sign : {-1, +1})
	{
	  auto f = sign * (x - root1) * (x - root2) + 0 * x*x*x*x*x;

	  double min_pos_root = HUGE_VAL;
	  if (root1 >= 0)
	    min_pos_root = std::min(min_pos_root, root1);
	  if (root2 >= 0)
	    min_pos_root = std::min(min_pos_root, root2);
	  
	  double bound = LMQ_lower_bound(f);
	  if (min_pos_root == HUGE_VAL)
	    BOOST_CHECK_EQUAL(bound, HUGE_VAL);
	  else
	    BOOST_CHECK(bound <= min_pos_root);
	}

  //Test constant coefficients 
  BOOST_CHECK_EQUAL(LMQ_lower_bound(1 + 0 * x*x*x*x*x), HUGE_VAL);
}

template<class Real, size_t Order1, size_t Order2>
void check_roots(magnet::containers::StackVector<Real, Order1> sol, magnet::containers::StackVector<Real, Order2> standard)
{
  BOOST_CHECK_EQUAL(sol.size(), standard.size());
  if (sol.size() != standard.size())
    return;
  std::sort(standard.begin(), standard.end());
  
  for(size_t i(0); i < standard.size(); ++i)
    BOOST_CHECK_CLOSE(sol[i], standard[i], 1e-11);
}


BOOST_AUTO_TEST_CASE( generic_solve_real_roots )
{
  using namespace magnet::math;
  const Polynomial<1> x{0, 1};

  const double roots[] = {-1e5, -0.14159265, 3.14159265, -0.0001,0.1, 0.3333, 0.6, 1.001, 2.0, 3.14159265, 1e7};

  for (const double root1: roots)
    for (const double root2: roots)
      if (root1 != root2)
	for (const double root3: roots)
	  if (!std::set<double>{root1, root2}.count(root3))
	    for (const double root4: roots)
	      if (!std::set<double>{root1, root2, root3}.count(root4))
		for (const double root5: roots)
		  if (!std::set<double>{root1, root2, root3, root4}.count(root5))
		    for (int sign : {-1, +1})
		      {
			magnet::containers::StackVector<double, 5> test_roots{root1, root2, root3, root4, root5};

			//Test where all 5 roots of a 5th order Polynomial are real
			auto f1 = sign * (x - root1) * (x - root2) * (x - root3) * (x - root4) * (x - root5);
			//Test where 5 roots of a 7th order Polynomial are real
			auto f2 = f1 * (x * x - 3 * x + 4);

			//Test the default implementation (VAS + TOMS748)
			check_roots(solve_real_roots(f1), test_roots);
			check_roots(solve_real_roots(f2), test_roots);
			
			//Test VCA + TOMS748
			check_roots(solve_real_roots<PolyRootBounder::VCA, PolyRootBisector::TOMS748>(f1), test_roots);
			check_roots(solve_real_roots<PolyRootBounder::VCA, PolyRootBisector::TOMS748>(f2), test_roots);
			
			//Test VAS + BISECTION
			check_roots(solve_real_roots<PolyRootBounder::VAS, PolyRootBisector::BISECTION>(f1), test_roots);
			check_roots(solve_real_roots<PolyRootBounder::VAS, PolyRootBisector::BISECTION>(f2), test_roots);

			//Test the next_root implementation
			std::vector<double> pos_roots;
			for (double root : test_roots)
			  if (root >= 0)
			    pos_roots.push_back(root);
			
			std::sort(pos_roots.begin(), pos_roots.end());
			BOOST_CHECK_CLOSE(pos_roots[0], next_root(f1), 1e-11);
			BOOST_CHECK_CLOSE(pos_roots[0], next_root(f2), 1e-11);
		      }
}
