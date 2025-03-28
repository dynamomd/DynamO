#define BOOST_TEST_MODULE StackVector_test
#include <boost/test/included/unit_test.hpp>
#include <magnet/containers/stack_vector.hpp>

using namespace magnet::containers;

BOOST_AUTO_TEST_CASE(StackVector_size) {
  StackVector<int, 3> vec;
  BOOST_CHECK(vec.size() == 0);
  BOOST_CHECK(vec.empty());

  vec.push_back(1);
  BOOST_CHECK(vec.size() == 1);
  BOOST_CHECK(!vec.empty());

  vec.push_back(2);
  BOOST_CHECK(vec.size() == 2);
  BOOST_CHECK(!vec.empty());
}

BOOST_AUTO_TEST_CASE(StackVector_initializer_list) {
  StackVector<double, 3> vec1{};
  BOOST_CHECK(vec1.size() == 0);
  BOOST_CHECK(vec1.empty());

  StackVector<double, 3> vec2{0.5, 0.25};
  BOOST_CHECK(vec2.size() == 2);
  BOOST_CHECK(!vec2.empty());
  BOOST_CHECK_EQUAL(vec2[0], 0.5);
  BOOST_CHECK_EQUAL(vec2[1], 0.25);

  StackVector<double, 3> vec3{0.5, 0.25, 0.125};
  BOOST_CHECK(vec3.size() == 3);
  BOOST_CHECK(!vec3.empty());
  BOOST_CHECK_EQUAL(vec3[0], 0.5);
  BOOST_CHECK_EQUAL(vec3[1], 0.25);
  BOOST_CHECK_EQUAL(vec3[2], 0.125);

  StackVector<double, 3> vec4{0.5, 0.25, 0.125, 0.1};
  BOOST_CHECK(vec4.size() == 3);
  BOOST_CHECK(!vec4.empty());
  BOOST_CHECK_EQUAL(vec4[0], 0.5);
  BOOST_CHECK_EQUAL(vec4[1], 0.25);
  BOOST_CHECK_EQUAL(vec4[2], 0.125);
}

BOOST_AUTO_TEST_CASE(StackVector_foreach) {
  StackVector<double, 3> vec4{0.5, 0.25, 0.125};

  double sum = 0.0;
  for (const auto &val : vec4)
    sum += val;
  BOOST_CHECK_CLOSE(sum, 0.875, 0.001);

  for (auto &val : vec4)
    val *= 2;
  BOOST_CHECK_CLOSE(vec4[0], 1.0, 0.001);
}
