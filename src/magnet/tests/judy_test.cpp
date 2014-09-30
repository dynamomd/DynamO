#define BOOST_TEST_MODULE JudyContainer_test
#include <boost/test/included/unit_test.hpp>
#include <iostream>
#include <magnet/containers/judy.hpp>
#include <cmath>
#include <set>
#include <map>

const size_t N = 10;
const size_t Nrange = 10000;

auto rng = std::mt19937(std::random_device()());
auto IDgen = std::bind(std::uniform_int_distribution<size_t>(0, Nrange), rng);

#define CompareContainers(test, reference)				\
  {BOOST_CHECK_EQUAL(test.size(), reference.size());			\
  typename std::remove_reference<decltype(reference)>::type copy(test.begin(), test.end()); \
  BOOST_CHECK(copy == reference);}

BOOST_AUTO_TEST_CASE( Judy_set )
{
  using namespace magnet::containers;
  JudySet<size_t> test;
  std::set<size_t> reference;
  
  //Test empty containers
  BOOST_CHECK(test.begin() == test.end());

  //Test inserts
  for (size_t i(0); i < N; ++i)
    {
      size_t value(IDgen());
      test.insert(value);
      reference.insert(value);
    }

  CompareContainers(test, reference);
  
  //Test erase of  1/10th of known elements
  for (size_t i(0); i < N / 10; ++i)
    {
      test.erase(*reference.begin());
      reference.erase(*reference.begin());
    }

  CompareContainers(test, reference);

  //Test erase of randomly generated elements
  for (size_t i(0); i < N; ++i)
    {
      size_t value(IDgen());
      test.erase(value);
      reference.erase(value);
    }

  CompareContainers(test, reference);

  //Test random access
  std::set<size_t> test2;
  for (size_t i(0); i < test.size(); ++i)
    test2.insert(*test.findNth(i));

  CompareContainers(test, test2);

  //Test clearing
  test.clear();
  BOOST_CHECK(test.size() == 0);
  BOOST_CHECK(test.begin() == test.end());  
}

BOOST_AUTO_TEST_CASE( Judy_map )
{
  using namespace magnet::containers;
  JudyMap<size_t, size_t> test;
  std::map<size_t, size_t> reference;

  //Test empty containers
  BOOST_CHECK(test.begin() == test.end());
  BOOST_CHECK(test.empty());
  BOOST_CHECK(test.size() == 0);

  //Test inserts
  for (size_t i(0); i < N; ++i)
    {
      JudyMap<size_t, size_t>::value_type value(IDgen(), IDgen());
      test.insert(value);
      reference.insert(value);
    }

  BOOST_CHECK(test.begin() != test.end());

  CompareContainers(test, reference);

  //Test erase of  1/10th of known elements
  for (size_t i(0); i < N / 10; ++i)
    {
      test.erase(reference.begin()->first);
      reference.erase(reference.begin()->first);
    }

  CompareContainers(test, reference);

  //Test erase of randomly generated elements
  for (size_t i(0); i < N; ++i)
    {
      size_t key(IDgen());
      test.erase(key);
      reference.erase(key);
    }

  CompareContainers(test, reference);

  //Test sequential findNth access
  decltype(reference) random_copy;
  for (size_t i(0); i < test.size(); ++i)
    random_copy.insert(*test.findNth(i));

  CompareContainers(test, random_copy);

  //Test array access operator read
  random_copy.clear();
  for (const auto& value : test)
    random_copy.insert(JudyMap<size_t, size_t>::value_type(value.first, test[value.first]));

  CompareContainers(test, random_copy);

  //Test clearing
  test.clear();
  BOOST_CHECK(test.size() == 0);
  BOOST_CHECK(test.empty());
  BOOST_CHECK(test.begin() == test.end());
}
