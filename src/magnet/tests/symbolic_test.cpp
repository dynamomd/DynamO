#define BOOST_TEST_MODULE Symbolic_math_test
#include <boost/test/included/unit_test.hpp>
#include <magnet/math/trigsymbols.hpp>


using namespace magnet::math;
const Polynomial<1> x{0, 1};

template<class T1, class T2>
bool compare_expression(const T1& f, const T2& g) {
  std::ostringstream os;
  os << f;
  std::string f_str = os.str();
  os.str(""); os.clear();
  os << g;
  std::string g_str = os.str();
  return f_str == g_str;
}

BOOST_AUTO_TEST_CASE( expand_null )
{
  //Test that expand does nothing when it has nothing to do
  auto poly1 = expand(2 * x * x);
  BOOST_CHECK_EQUAL(poly1[0], 0);
  BOOST_CHECK_EQUAL(poly1[1], 0);
  BOOST_CHECK_EQUAL(poly1[2], 2);    
}
  
BOOST_AUTO_TEST_CASE( expand_polynomials )
{
  //Test addition and simplification of Polynomials
  auto poly1 = expand(add(2 * x * x, x));
  //This should become a Polynomial class, with its coefficients
  //accessible by the array operator.
  BOOST_CHECK_EQUAL(poly1[0], 0);
  BOOST_CHECK_EQUAL(poly1[1], 1);
  BOOST_CHECK_EQUAL(poly1[2], 2);
}
  
BOOST_AUTO_TEST_CASE( polynomials_derivative_addition )
{
  //Test Polynomial derivatives on addition Operation types
  auto poly1 = derivative(add(2 * x * x,  x));
  //derivative will automatically combine polynomials
  BOOST_CHECK_EQUAL(poly1[0], 1);
  BOOST_CHECK_EQUAL(poly1[1], 4);
}

BOOST_AUTO_TEST_CASE( polynomials_derivative_subtraction )
{
  //Test Polynomial derivatives on subtraction Operation types
  auto poly1 = derivative(subtract(2 * x * x,  x));
  //derivative will automatically combine polynomials
  BOOST_CHECK_EQUAL(poly1[0], -1);
  BOOST_CHECK_EQUAL(poly1[1], 4);
}

BOOST_AUTO_TEST_CASE( polynomials_multiply_expansion )
{
  //Test Polynomial simplification on multiplication Operation types
  auto poly1 = expand(multiply(x + 1,  x + 3));

  //derivative will automatically combine polynomials
  BOOST_CHECK_EQUAL(poly1[0], 3);
  BOOST_CHECK_EQUAL(poly1[1], 4);
  BOOST_CHECK_EQUAL(poly1[2], 1);
}

BOOST_AUTO_TEST_CASE( function_basic )
{
  //Check basic Function operation
  BOOST_CHECK_CLOSE(magnet::math::sin(x)(0.5), std::sin(0.5), 1e-10);
  BOOST_CHECK_CLOSE(magnet::math::cos(x)(0.5), std::cos(0.5), 1e-10);
}

BOOST_AUTO_TEST_CASE( function_poly_multiplication )
{
  //Check function and Polynomial multiplication
  auto poly1 = sin(x + x) * x;
  BOOST_CHECK_CLOSE(poly1(0.5), std::sin(2 * 0.5) * 0.5, 1e-10);
  auto poly2 = x * sin(x + x);
  BOOST_CHECK_CLOSE(poly2(0.5), std::sin(2 * 0.5) * 0.5, 1e-10);
}

BOOST_AUTO_TEST_CASE( function_poly_derivatives )
{
  //Check function and Polynomial derivatives
  auto poly1 = derivative(x * sin(x));
  BOOST_CHECK_CLOSE(poly1(0.5), std::sin(0.5) + 0.5 * std::cos(0.5), 1e-10);
  auto poly2 = derivative(x * cos(x));
  BOOST_CHECK_CLOSE(poly2(0.5), -0.5 * std::sin(0.5) + std::cos(0.5), 1e-10);
}

BOOST_AUTO_TEST_CASE( function_poly_derivatives_special )
{ //Check special case derivatives of Functions with constant
  //arguments.
  auto poly1 = derivative(sin(Polynomial<0>{1}));
  BOOST_CHECK_EQUAL(poly1[0], 0);
  auto poly2 = derivative(cos(Polynomial<0>{1}));
  BOOST_CHECK_EQUAL(poly2[0], 0);
}

BOOST_AUTO_TEST_CASE( poly_specialised_multiply )
{
  //Check the specialised multiply operators are consistently
  //simplifying statements with Polynomial and Function types.

  //First check that negative matches are correctly determined by
  //compare_expression
  BOOST_CHECK(!compare_expression(x, sin(x)));

  //Now test
  BOOST_CHECK(compare_expression((sin(2*x) * x) * x, x * x * sin(2*x)));
  BOOST_CHECK(compare_expression((x * sin(2*x)) * x, x * x * sin(2*x)));
  BOOST_CHECK(compare_expression(x * (sin(2*x) * x), x * x * sin(2*x)));
  BOOST_CHECK(compare_expression(x * (x * sin(2*x)), x * x * sin(2*x)));
}
