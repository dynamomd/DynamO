#define BOOST_TEST_MODULE Numeric_math_test
#include <boost/test/included/unit_test.hpp>
#include <magnet/math/numeric.hpp>
#include <magnet/math/symbolic.hpp>

using namespace magnet::math;
const Variable<'x'> x = Variable<'x'>();

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

BOOST_AUTO_TEST_CASE( NewtonRaphson_root )
{
  //Simple check for positive roots, using manual derivatives
  auto f = simplify(x*x - 4.0);
  auto df = derivative(f, x);
  double xroot = 6.0;
  BOOST_CHECK(newton_raphson([&](double x){ return std::array<double, 2>{{eval(f, x), eval(df, x)}}; }, xroot));
  BOOST_CHECK_CLOSE(xroot, 2.0, 1e-10);
  
  xroot = 6.0;
  BOOST_CHECK(newton_raphson([&](double x){ return eval_derivatives<1>(f, x); }, xroot));
  BOOST_CHECK_CLOSE(xroot, 2.0, 1e-10);

  //Test "difficult" equations for NR.
  auto f2 = simplify(x*x*x - 2.0 * x + 2.0);
  xroot = 0.0;
  //This oscillates between 0 and 1, check that an early exit is found!
  BOOST_CHECK(!newton_raphson([&](double x){ return eval_derivatives<1>(f2, x); }, xroot));
}

BOOST_AUTO_TEST_CASE( Halley_root )
{
  //Simple check for positive roots, using manual derivatives
  auto f = simplify(x*x - 4.0);
  auto df = derivative(f, x);
  auto ddf = derivative(df, x);

  double xroot = 6.0;
  BOOST_CHECK(halleys_method([&](double x){ return std::array<double, 3>{{eval(f, x), eval(df, x), eval(ddf, x)}}; }, xroot));
  BOOST_CHECK_CLOSE(xroot, 2.0, 1e-10);
  
  xroot = 6.0;
  BOOST_CHECK(halleys_method([&](double x){ return eval_derivatives<2>(f, x); }, xroot));
  BOOST_CHECK_CLOSE(xroot, 2.0, 1e-10);

  //Test "difficult" equations for NR.
  auto f2 = simplify(1.0 - x*x);
  xroot = 0.0;
  BOOST_CHECK(!halleys_method([&](double x){ return eval_derivatives<2>(f2, x); }, xroot));

  xroot = 0.01;
  BOOST_CHECK(halleys_method([&](double x){ return eval_derivatives<2>(f2, x); }, xroot));
  BOOST_CHECK_CLOSE(xroot, 1.0, 1e-10);  
}
