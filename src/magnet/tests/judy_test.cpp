#define BOOST_TEST_MODULE JudyContainer_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <magnet/containers/judy.hpp>
#include <iostream>
#include <cmath>
#include <set>
#include <map>

const size_t N = 10;
const size_t Nrange = 10000;

auto rng = std::mt19937(std::random_device()());
auto IDgen = std::bind(std::uniform_int_distribution<size_t>(0, Nrange), rng);

template<class C1, class C2>
void compareContainers(const C1& test, const C2& reference)
{
  BOOST_CHECK_EQUAL(test.size(), reference.size());
  typename std::remove_reference<decltype(reference)>::type copy(test.begin(), test.end());
  BOOST_CHECK(copy == reference);
}

template<class C1>
void compareContainers(const C1& test, const C1& reference)
{ BOOST_CHECK(test == reference); }

BOOST_AUTO_TEST_CASE( Judy_pair_set )
{
  magnet::containers::JudyPairSet test;
  std::set<std::pair<size_t, size_t> > reference;
  
  //Test empty containers
  BOOST_CHECK(test.begin() == test.end());

  //Test inserts
  for (size_t i(0); i < N; ++i)
    {
      std::pair<size_t, size_t> pair(IDgen(), IDgen());
      test.insert(pair);
      if (pair.first > pair.second) std::swap(pair.first, pair.second);
      reference.insert(pair);
    }

  compareContainers(test, reference);
  
  //Test erase of  1/10th of known elements
  for (size_t i(0); i < N / 10; ++i)
    {
      test.erase(*reference.begin());
      reference.erase(*reference.begin());
    }

  compareContainers(test, reference);

  //Test erase of randomly generated elements
  for (size_t i(0); i < N; ++i)
    {
      std::pair<size_t, size_t> pair(IDgen(), IDgen());
      test.erase(pair);
      reference.erase(pair);
    }

  compareContainers(test, reference);

  //Test random access
  std::set<std::pair<size_t, size_t> > test2;
  for (size_t i(0); i < test.size(); ++i)
    test2.insert(test[i]);

  compareContainers(test, test2);

  //Test clearing
  test.clear();
  BOOST_CHECK(test.size() == 0);
  BOOST_CHECK(test.begin() == test.end());  
}

BOOST_AUTO_TEST_CASE( Judy_pair_map )
{
  using namespace magnet::containers;
  JudyPairMap test;
  std::map<JudyPairMap::key_type, JudyPairMap::mapped_type> reference;

  //Test empty containers
  BOOST_CHECK(test.begin() == test.end());

  //Test inserts
  for (size_t i(0); i < N; ++i)
    {
      JudyPairMap::value_type value(JudyPairMap::key_type(IDgen(), IDgen()), IDgen() + 1);
      test.insert(value);
      if (value.first.first > value.first.second) std::swap(value.first.first, value.first.second);
      reference.insert(value);
    }

  compareContainers(test, reference);

  //Test erase of  1/10th of known elements
  for (size_t i(0); i < N / 10; ++i)
    {
      test.erase(reference.begin()->first);
      reference.erase(reference.begin()->first);
    }

  compareContainers(test, reference);

  //Test erase of randomly generated elements
  for (size_t i(0); i < N; ++i)
    {
      JudyPairMap::key_type pair(IDgen(), IDgen());
      test.erase(pair);
      reference.erase(pair);
    }

  compareContainers(test, reference);

  //Test random access
  std::map<std::pair<size_t, size_t>, size_t> test2;
  for (size_t i(0); i < test.size(); ++i)
    test2.insert(test[i]);

  compareContainers(test, test2);

  //Test clearing
  test.clear();
  BOOST_CHECK(test.size() == 0);
  BOOST_CHECK(test.begin() == test.end());  
}
