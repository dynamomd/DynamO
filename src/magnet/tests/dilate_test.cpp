#define BOOST_TEST_MODULE Dilate_test
#include <algorithm>
#include <boost/test/included/unit_test.hpp>
#include <iomanip>
#include <magnet/math/dilated_int.hpp>
#include <random>
#include <sstream>

// std::string hexout(size_t val)
//{
//   std::ostringstream os;
//   os << "0x"<< std::hex << std::setw(sizeof(size_t)*2) << std::setfill ('0')
//   << val; return os.str();
// }
//
// template<size_t d, size_t i>
// std::string undilate_text()
//{
//   typedef magnet::math::dilatedinteger::c<i,d> c;
//   typedef magnet::math::dilatedinteger::z<i,d> z;
//
//   std::ostringstream os;
//   os << "\nUndilate<" << d <<"> round " << i << " ( val * " <<
//   hexout(c::result) << ") & " << hexout(z::result); return os.str();
// }
//
// template<size_t d, size_t i>
// std::string dilate_text()
//{
//   typedef magnet::math::dilatedinteger::b<i,d> b;
//   typedef magnet::math::dilatedinteger::y<i,d> y;
//
//   std::ostringstream os;
//   os << "Dilate<" << d << "> round " << i << " ( val * " << hexout(b::result)
//   << ") & " << hexout(y::result)
//      << "\n";
//   return os.str();
// }
//
// template<size_t i>
// std::string dilate_text_2()
//{
//   static size_t shiftval =
//   magnet::math::dilatedinteger::dilate<2>::shiftval<i>::result;
//
//   std::ostringstream os;
//   os << "Dilate<2> round " << i << " ( val | (val <<  " << shiftval << ")) &
//   "
//      << hexout(magnet::math::dilatedinteger::y<i,2>::result)
//      << "\n";
//   return os.str();
// }

template <size_t d> void test() {
  const size_t maxval =
      magnet::math::dilatedinteger::maxDilatableValue<d>::result;
  std::default_random_engine RNG;
  RNG.seed(std::random_device()());
  std::uniform_int_distribution<size_t> dist(0, maxval);

  for (size_t test(0); test < 1000000; ++test) {
    size_t i = dist(RNG);
    BOOST_CHECK_MESSAGE(
        i == magnet::math::undilate<d>(magnet::math::dilate<d>(i)),
        "Failed to dilate/undilate " << i << "\n");
  }

  BOOST_CHECK_MESSAGE(
      maxval == magnet::math::undilate<d>(magnet::math::dilate<d>(maxval)),
      "Undilate and dilate are not symmetric for the max value");
}

BOOST_AUTO_TEST_CASE(dilate2) { test<2>(); }
BOOST_AUTO_TEST_CASE(dilate3) { test<3>(); }
BOOST_AUTO_TEST_CASE(dilate4) { test<4>(); }
BOOST_AUTO_TEST_CASE(dilate5) { test<5>(); }
BOOST_AUTO_TEST_CASE(dilate6) { test<6>(); }
BOOST_AUTO_TEST_CASE(dilate7) { test<7>(); }
BOOST_AUTO_TEST_CASE(dilate8) { test<8>(); }
BOOST_AUTO_TEST_CASE(dilate9) { test<9>(); }
BOOST_AUTO_TEST_CASE(dilate10) { test<10>(); }
