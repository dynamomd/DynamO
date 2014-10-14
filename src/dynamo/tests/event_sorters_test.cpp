#define BOOST_TEST_MODULE Sorters
#include <boost/test/included/unit_test.hpp>
#include <dynamo/eventtypes.hpp>
#include <random>
std::mt19937 RNG;

dynamo::Event genInteractionEvent(size_t N, const double meanFreeTime = 1.0, const size_t interactions = 2, size_t p1ID = std::numeric_limits<size_t>::max())
{
  const double dt = - meanFreeTime * std::log(1.0 - std::uniform_real_distribution<>()(RNG));
  std::uniform_int_distribution<size_t> pIDdist(0, N-1);
  std::uniform_int_distribution<size_t> iIDdist(0, interactions);

  if (p1ID == std::numeric_limits<size_t>::max())
    p1ID = pIDdist(RNG);

  size_t p2ID = p1ID;
  while (p1ID==p2ID)
    p2ID = pIDdist(RNG);
  
  return dynamo::Event(p1ID, dt, dynamo::INTERACTION, dynamo::CORE, iIDdist(RNG), p2ID);
}

template<class T>
std::vector<dynamo::Event> fillSorter(T& sorter, const size_t totalEvents, const size_t N) {
  std::vector<dynamo::Event> standard;
  for (size_t i(0); i < totalEvents; ++i)
    {
      const dynamo::Event e = genInteractionEvent(N);
      sorter.push(e);
      standard.push_back(e);
    }
  return standard;
}

#include <dynamo/schedulers/sorters/heapPEL.hpp>
#include <dynamo/schedulers/sorters/MinMaxPEL.hpp>
typedef boost::mpl::list<dynamo::HeapPEL,
			 dynamo::MinMaxPEL<2>,
			 dynamo::MinMaxPEL<3>,
			 dynamo::MinMaxPEL<4>
			 > PEL_types;

BOOST_AUTO_TEST_CASE_TEMPLATE(PEL_standard, T, PEL_types){
  RNG.seed(std::random_device()());
  T sorter;
  const size_t N=100;

  BOOST_REQUIRE(sorter.empty());
  BOOST_REQUIRE(sorter.size()==0);
  sorter.push(genInteractionEvent(N));
  BOOST_REQUIRE(!sorter.empty());
  BOOST_CHECK(sorter.size()==1);
  sorter.clear();
  BOOST_CHECK(sorter.empty());
  BOOST_CHECK(sorter.size()==0);

  const size_t totalEvents = 10;
  
  //Check that the sorter returns events in the correct order
  std::vector<dynamo::Event> standard = fillSorter(sorter, totalEvents, N);
  
  std::sort(standard.begin(), standard.end());  
  while (!standard.empty() && !sorter.empty()) {
    const dynamo::Event correct=standard.front();
    if (sorter.top()._type == dynamo::RECALCULATE) {
      BOOST_CHECK(correct._dt == sorter.top()._dt);
      sorter.pop();
      BOOST_CHECK(sorter.empty());
      break;
    }
    else
      BOOST_CHECK(correct == sorter.top());
    
    sorter.pop();
    standard.erase(standard.begin());
  }
}

#include <dynamo/schedulers/sorters/referenceFEL.hpp>
#include <dynamo/schedulers/sorters/CBTFEL.hpp>
#include <dynamo/schedulers/sorters/boundedPQFEL.hpp>
typedef boost::mpl::list<
  dynamo::ReferenceFEL
  ,dynamo::CBTFEL<dynamo::HeapPEL>
  ,dynamo::CBTFEL<dynamo::MinMaxPEL<2> >
  ,dynamo::CBTFEL<dynamo::MinMaxPEL<5> >
  ,dynamo::CBTFEL<dynamo::MinMaxPEL<30> >
  ,dynamo::BoundedPQFEL<dynamo::HeapPEL>
  ,dynamo::BoundedPQFEL<dynamo::MinMaxPEL<2> >
  ,dynamo::BoundedPQFEL<dynamo::MinMaxPEL<5> >
  ,dynamo::BoundedPQFEL<dynamo::MinMaxPEL<30> >
			 > FEL_types;

#define validateEvents(e1, e2) \
  BOOST_REQUIRE_CLOSE(e1._dt, e2._dt, 1e-7);\
  BOOST_REQUIRE_EQUAL(e1._particle1ID, e2._particle1ID);\
  BOOST_REQUIRE_EQUAL(e1._sourceID, e2._sourceID);\
  BOOST_REQUIRE_EQUAL(e1._type, e2._type);\
  BOOST_REQUIRE_EQUAL(e1._additionalData1, e2._additionalData1);\
  BOOST_REQUIRE_EQUAL(e1._source, e2._source);\
  if (e1._source != dynamo::INTERACTION) BOOST_REQUIRE_EQUAL(e1._additionalData2, e2._additionalData2);

BOOST_AUTO_TEST_CASE_TEMPLATE(FEL_standard, T, FEL_types){
  RNG.seed(std::random_device()());
  const size_t N = 100;
  const size_t eventsPerParticle = 10;
  T FEL;
  std::vector<dynamo::Event> reference;
  
  //Test empty(), clear().
  FEL.init(N);
  for (size_t i(0); i < 100; ++i)
    FEL.push(genInteractionEvent(N, 1.0, 1));
  BOOST_REQUIRE(!FEL.empty());
  FEL.clear();
  BOOST_REQUIRE(FEL.empty());

  //Load and fully drain the queue, testing that the FEL is fully sorted
  FEL.init(N);
  for (size_t i(0); i < N * eventsPerParticle; ++i) {
    const dynamo::Event e = genInteractionEvent(N, 1.0, 1);
    reference.push_back(e);
    FEL.push(e);
  }
  while (!reference.empty()) {
    BOOST_REQUIRE(!FEL.empty());
    const auto next_it = std::min_element(reference.begin(), reference.end());
    const dynamo::Event nextEvent = *next_it;
    const dynamo::Event testEvent = FEL.top();

    if (testEvent._type == dynamo::RECALCULATE) {
      FEL.pop();
      for (const dynamo::Event& e: reference)
	if (e._particle1ID == testEvent._particle1ID)
	  FEL.push(e);
      continue;
    }

    validateEvents(nextEvent, testEvent);
    reference.erase(next_it);
    FEL.pop();
  }
  BOOST_REQUIRE(FEL.empty());

  //Test scaling of event times
  FEL.clear();
  FEL.init(N);
  for (size_t i(0); i < N * eventsPerParticle; ++i) {
    const dynamo::Event e = genInteractionEvent(N, 1.0, 1);
    reference.push_back(e);
    FEL.push(e);
  }

  const double factor = 3.141;
  FEL.rescaleTimes(factor);
  for (dynamo::Event& e: reference)
    e._dt *= factor;
  while (!reference.empty()) {
    BOOST_REQUIRE(!FEL.empty());
    const auto next_it = std::min_element(reference.begin(), reference.end());
    const dynamo::Event nextEvent = *next_it;
    const dynamo::Event testEvent = FEL.top();

    if (testEvent._type == dynamo::RECALCULATE) {
      FEL.pop();
      for (const dynamo::Event& e: reference)
	if (e._particle1ID == testEvent._particle1ID)
	  FEL.push(e);
      continue;
    }

    validateEvents(nextEvent,testEvent);
    reference.erase(next_it);
    FEL.pop();
  }
  BOOST_REQUIRE(FEL.empty());

  //Test streaming of event times
  FEL.clear();
  FEL.init(N);
  for (size_t i(0); i < N * eventsPerParticle; ++i) {
    const dynamo::Event e = genInteractionEvent(N, 1.0, 1);
    reference.push_back(e);
    FEL.push(e);
  }
  FEL.stream(factor);
  for (dynamo::Event& e: reference)
    e._dt -= factor;
  while (!reference.empty()) {
    BOOST_REQUIRE(!FEL.empty());
    const auto next_it = std::min_element(reference.begin(), reference.end());
    const dynamo::Event nextEvent = *next_it;
    const dynamo::Event testEvent = FEL.top();

    if (testEvent._type == dynamo::RECALCULATE) {
      FEL.pop();
      for (const dynamo::Event& e: reference)
	if (e._particle1ID == testEvent._particle1ID)
	  FEL.push(e);
      continue;
    }

    validateEvents(nextEvent,testEvent);
    reference.erase(next_it);
    FEL.pop();
  }
  BOOST_REQUIRE(FEL.empty());

  //Test particle invalidation
  FEL.clear();
  FEL.init(N);
  for (size_t i(0); i < N * eventsPerParticle; ++i) {
    const dynamo::Event e = genInteractionEvent(N, 1.0, 1);
    reference.push_back(e);
    FEL.push(e);
  }
  while (!reference.empty()) {
    const auto next_it = std::min_element(reference.begin(), reference.end());
    const dynamo::Event nextEvent = *next_it;
    const dynamo::Event testEvent = FEL.top();

    if (testEvent._type == dynamo::RECALCULATE) {
      FEL.pop();
      for (const dynamo::Event& e: reference)
    	if (e._particle1ID == testEvent._particle1ID)
    	  FEL.push(e);
      continue;
    }
    
    validateEvents(nextEvent, testEvent);
    
    ////Delete the events from both the reference and the FEL
    auto test = [=](const dynamo::Event& e){ 
      return (e._particle1ID == testEvent._particle1ID) || ((e._source == dynamo::INTERACTION) && (e._particle2ID == testEvent._particle1ID));
    };
    reference.erase(std::remove_if(reference.begin(), reference.end(), test), reference.end());
    
    FEL.invalidate(testEvent._particle1ID);
  }

  //Perform a mock simulation
  FEL.clear();
  FEL.init(N);
  for (size_t i(0); i < N * eventsPerParticle; ++i) {
    const dynamo::Event e = genInteractionEvent(N, 1.0, 1);
    reference.push_back(e);
    FEL.push(e);
  }
  for (size_t i(0); (i < N * eventsPerParticle) && (!reference.empty()); ++i) {
    const auto next_it = std::min_element(reference.begin(), reference.end());
    const dynamo::Event nextEvent = *next_it;
    const dynamo::Event testEvent = FEL.top();
    
    if (testEvent._type == dynamo::RECALCULATE) {
      FEL.pop();
      for (const dynamo::Event& e: reference)
    	if (e._particle1ID == testEvent._particle1ID)
    	  FEL.push(e);
      continue;
    }
    
    validateEvents(nextEvent, testEvent);
    
    ////Delete the events from both the reference and the FEL
    auto test = [=](const dynamo::Event& e){
      return (e._particle1ID == testEvent._particle1ID) || (e._particle1ID == testEvent._particle2ID)
      || ((e._source == dynamo::INTERACTION) 
	  && ((e._particle2ID == testEvent._particle1ID) || (e._particle2ID == testEvent._particle2ID)));
    };
    reference.erase(std::remove_if(reference.begin(), reference.end(), test), reference.end());
    
    FEL.invalidate(testEvent._particle1ID);
    FEL.invalidate(testEvent._particle2ID);
  
    //Stream the lists forward.
    FEL.stream(testEvent._dt);
    for (dynamo::Event& e: reference)
      e._dt -= testEvent._dt;
    
    //Add the events to the FEL.
    for (size_t j(0); j < eventsPerParticle; j++) {
      const dynamo::Event newEvent = genInteractionEvent(N, 1.0, 1, testEvent._particle1ID);
      FEL.push(newEvent);
      reference.push_back(newEvent);
    }
    
    for (size_t j(0); j < eventsPerParticle; j++) {
      const dynamo::Event newEvent = genInteractionEvent(N, 1.0, 1, testEvent._particle2ID);
      FEL.push(newEvent);
      reference.push_back(newEvent);
    }
  }
}
