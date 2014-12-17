#define BOOST_TEST_MODULE Symbolic_math_test
#include <boost/test/included/unit_test.hpp>
#include <magnet/math/polynomial.hpp>
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
  if (!(f_str == g_str))
    std::cerr << f << " != " << g << std::endl;
  return f_str == g_str;
}

BOOST_AUTO_TEST_CASE( Substitution_of_variables )
{
  Variable<'x'> x;
  Variable<'y'> y;
  
  BOOST_CHECK(compare_expression(substitution(x, x==y), "y"));
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
  auto poly1 = derivative(add(2 * x * x,  x), Variable<'x'>());
  //derivative will automatically combine polynomials
  BOOST_CHECK_EQUAL(poly1[0], 1);
  BOOST_CHECK_EQUAL(poly1[1], 4);
}

BOOST_AUTO_TEST_CASE( polynomials_derivative_subtraction )
{
  //Test Polynomial derivatives on subtraction Operation types
  auto poly1 = derivative(subtract(2 * x * x,  x), Variable<'x'>());
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
  BOOST_CHECK_CLOSE(eval(magnet::math::sin(x), 0.5), std::sin(0.5), 1e-10);
  BOOST_CHECK_CLOSE(eval(magnet::math::cos(x), 0.5), std::cos(0.5), 1e-10);

  //Test BinaryOP Addition and subtraction
  BOOST_CHECK_CLOSE(eval(x * sin(x) + x, 0.5), 0.5 * std::sin(0.5) + 0.5, 1e-10);
  BOOST_CHECK_CLOSE(eval(x * sin(x) - x, 0.5), 0.5 * std::sin(0.5) - 0.5, 1e-10);
}

BOOST_AUTO_TEST_CASE( function_poly_multiplication )
{
  //Check function and Polynomial multiplication
  auto poly1 = sin(x + x) * x;
  BOOST_CHECK_CLOSE(eval(poly1, 0.5), std::sin(2 * 0.5) * 0.5, 1e-10);
  auto poly2 = x * sin(x + x);
  BOOST_CHECK_CLOSE(eval(poly2, 0.5), std::sin(2 * 0.5) * 0.5, 1e-10);
}

BOOST_AUTO_TEST_CASE( function_poly_derivatives )
{
  //Check function and Polynomial derivatives
  auto f1 = derivative(x * sin(x), Variable<'x'>());
  BOOST_CHECK_CLOSE(eval(f1, 0.5), std::sin(0.5) + 0.5 * std::cos(0.5), 1e-10);
  auto f2 = derivative(x * cos(x), Variable<'x'>());
  BOOST_CHECK_CLOSE(eval(f2, 0.5), -0.5 * std::sin(0.5) + std::cos(0.5), 1e-10);
}

BOOST_AUTO_TEST_CASE( function_poly_derivatives_special )
{ //Check special case derivatives of Functions with constant
  //arguments.
  BOOST_CHECK(compare_expression(derivative(sin(Polynomial<0>{1}), Variable<'x'>()), NullSymbol()));
  BOOST_CHECK(compare_expression(derivative(cos(Polynomial<0>{1}), Variable<'x'>()), NullSymbol()));
}

BOOST_AUTO_TEST_CASE( power_basic )
{
  //Check evaluation of powers
  BOOST_CHECK_CLOSE(eval(pow<3>(x), 4.0), 4.0*4.0*4.0, 1e-10);
  BOOST_CHECK_CLOSE(eval(pow<3>(x), 0.75), std::pow(0.75, 3), 1e-10);

  //Test PowerOp algebraic operations
  BOOST_CHECK_CLOSE(eval(pow<3>(x) - x, 0.75), std::pow(0.75, 3) - 0.75, 1e-10);
  BOOST_CHECK_CLOSE(eval(pow<3>(x) + x, 0.75), std::pow(0.75, 3) + 0.75, 1e-10);
  BOOST_CHECK_CLOSE(eval(pow<3>(x) * x, 0.75), std::pow(0.75, 3) * 0.75, 1e-10);

  //Check special case derivatives
  BOOST_CHECK(compare_expression(derivative(pow<1>(x), Variable<'x'>()), 1));
  BOOST_CHECK(compare_expression(derivative(pow<2>(x), Variable<'x'>()), 2 * x));

  //Check expansion
  BOOST_CHECK(compare_expression(expand(pow<3>(x+2)), (x+2) * (x+2) * (x+2)));;
}

BOOST_AUTO_TEST_CASE( Null_tests )
{
  //Check that Null symbols have zero derivative and value
  BOOST_CHECK(compare_expression(NullSymbol(), "Null"));
  BOOST_CHECK(compare_expression(derivative(NullSymbol(), Variable<'x'>()), "Null"));
  BOOST_CHECK_EQUAL(eval(NullSymbol(), 100), 0);

  //Check derivatives of constants becomes Null
  BOOST_CHECK(compare_expression(derivative(2, Variable<'x'>()), "Null"));
  BOOST_CHECK(compare_expression(derivative(3.141, Variable<'x'>()), "Null"));
  BOOST_CHECK(compare_expression(derivative(Vector{1,2,3}, Variable<'x'>()), "Null"));
}

BOOST_AUTO_TEST_CASE( Unity_tests )
{
  //Check that Null symbols have zero derivative and value
  BOOST_CHECK(compare_expression(UnitySymbol(), "Unity"));
  BOOST_CHECK(compare_expression(derivative(UnitySymbol(),Variable<'x'>()), "Null"));
  BOOST_CHECK_EQUAL(eval(UnitySymbol(), 100), 1);

  //Check derivatives of Unity
  BOOST_CHECK(compare_expression(derivative(UnitySymbol(), Variable<'x'>()), "Null"));

  //Check simplification of multiplication with Unity
  BOOST_CHECK(compare_expression(UnitySymbol() * UnitySymbol(), "Unity"));
  BOOST_CHECK(compare_expression(UnitySymbol() * 2, 2));
  BOOST_CHECK(compare_expression(UnitySymbol() * x, x));
  BOOST_CHECK(compare_expression(UnitySymbol() * Vector{1,2,3}, Vector{1,2,3}));
  BOOST_CHECK(compare_expression(2 * UnitySymbol(), 2));
  BOOST_CHECK(compare_expression(x * UnitySymbol() * x, x * x));
  BOOST_CHECK(compare_expression(Vector{1,2,3} * UnitySymbol(), Vector{1,2,3}));
}

BOOST_AUTO_TEST_CASE( Var_tests )
{
  Variable<'x'> x;
  Variable<'y'> y;

  BOOST_CHECK(compare_expression(x, "x")); 
  BOOST_CHECK(compare_expression(y, "y"));
  BOOST_CHECK(compare_expression(derivative(x, Variable<'x'>()), "Unity"));
  BOOST_CHECK(compare_expression(derivative(y, Variable<'x'>()), "Null"));
  BOOST_CHECK(compare_expression(derivative(y, Variable<'y'>()), "Unity"));
  BOOST_CHECK(compare_expression(substitution(x, x == 3.14159265), 3.14159265));

  //Check that substitutions in the wrong variable do nothing
  BOOST_CHECK(compare_expression(substitution(y, x == 3.14159265), "y"));

  //Check default substitution is for x
  BOOST_CHECK(compare_expression(eval(y, 3.14159265), "y"));

  //Check that Var derivatives are correct
  BOOST_CHECK(compare_expression(derivative(sin(x), Variable<'x'>()), cos(x)));

  //Check derivatives of Unity
  BOOST_CHECK(compare_expression(derivative(UnitySymbol(), Variable<'x'>()), "Null"));
  BOOST_CHECK(compare_expression(derivative(x, Variable<'x'>()), "Unity"));
  BOOST_CHECK(compare_expression(derivative(x * sin(x), Variable<'x'>()), sin(x) + x * cos(x)));
}

BOOST_AUTO_TEST_CASE( reorder_operations )
{
  //Check the specialised multiply operators are consistently
  //simplifying statements.

  //Again, check that negative matches are correctly determined by
  //compare_expression
  BOOST_CHECK(!compare_expression(x, sin(x)));

  //Here we're looking for the two Polynomial terms to be reordered 
  BOOST_CHECK(compare_expression((sin(2*x) * x) * x, x * x * sin(2*x)));
  BOOST_CHECK(compare_expression((x * sin(2*x)) * x, x * x * sin(2*x)));
  BOOST_CHECK(compare_expression(x * (sin(2*x) * x), x * x * sin(2*x)));
  BOOST_CHECK(compare_expression(x * (x * sin(2*x)), x * x * sin(2*x)));

  //Here we check that constants (such as 2) will become NullSymbol
  //types when the derivative is taken, causing their terms to be
  //eliminated.
  BOOST_CHECK(compare_expression(derivative(2 * cos(x), Variable<'x'>()), -2 * sin(x)));
  BOOST_CHECK(compare_expression(derivative(2 * sin(x), Variable<'x'>()), 2 * cos(x)));
}
