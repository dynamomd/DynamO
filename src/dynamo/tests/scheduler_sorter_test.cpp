#define BOOST_TEST_MODULE Scheduler_Sorter_test
#include <boost/test/included/unit_test.hpp>
#include <dynamo/BC/include.hpp>
#include <dynamo/dynamics/newtonian.hpp>
#include <dynamo/inputplugins/cells/include.hpp>
#include <dynamo/inputplugins/include.hpp>
#include <dynamo/interactions/hardsphere.hpp>
#include <dynamo/outputplugins/misc.hpp>
#include <dynamo/ranges/include.hpp>
#include <dynamo/schedulers/include.hpp>
#include <dynamo/schedulers/sorters/include.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/species/point.hpp>
#include <random>

std::mt19937 RNG;
typedef dynamo::FELBoundedPQ<dynamo::PELMinMax<3>> DefaultSorter;

template <class Scheduler, class Sorter> void runTest() {
  dynamo::Simulation Sim;

  RNG.seed(std::random_device()());
  Sim.ranGenerator.seed(std::random_device()());

  Sim.dynamics =
      dynamo::shared_ptr<dynamo::Dynamics>(new dynamo::DynNewtonian(&Sim));
  Sim.BCs = dynamo::shared_ptr<dynamo::BoundaryCondition>(
      new dynamo::BCPeriodic(&Sim));
  Sim.ptrScheduler =
      dynamo::shared_ptr<dynamo::Scheduler>(new Scheduler(&Sim, new Sorter()));
  Sim.primaryCellSize = dynamo::Vector(11, 11, 11);
  Sim.interactions.push_back(
      dynamo::shared_ptr<dynamo::Interaction>(new dynamo::IHardSphere(
          &Sim, 1.0, 1.0, new dynamo::IDPairRangeAll(), "Bulk")));
  Sim.addSpecies(dynamo::shared_ptr<dynamo::Species>(
      new dynamo::SpPoint(&Sim, new dynamo::IDRangeAll(&Sim), 1.0, "Bulk", 0)));

  Sim.particles.push_back(dynamo::Particle(dynamo::Vector{0.1, 0, 0},
                                           dynamo::Vector{0, 0, 0},
                                           Sim.particles.size()));
  Sim.particles.push_back(dynamo::Particle(dynamo::Vector{1.1, 0, 0},
                                           dynamo::Vector{1, 0, 0},
                                           Sim.particles.size()));
  Sim.particles.push_back(dynamo::Particle(dynamo::Vector{3.1, 0, 0},
                                           dynamo::Vector{0, 0, 0},
                                           Sim.particles.size()));
  Sim.particles.push_back(dynamo::Particle(dynamo::Vector{5.1, 0, 0},
                                           dynamo::Vector{0, 0, 0},
                                           Sim.particles.size()));
  Sim.particles.push_back(dynamo::Particle(dynamo::Vector{7.1, 0, 0},
                                           dynamo::Vector{0, 0, 0},
                                           Sim.particles.size()));
  Sim.particles.push_back(dynamo::Particle(dynamo::Vector{9.1, 0, 0},
                                           dynamo::Vector{0, 0, 0},
                                           Sim.particles.size()));

  Sim.ensemble = dynamo::Ensemble::loadEnsemble(Sim);
  Sim.endEventCount = 1000;
  Sim.addOutputPlugin("Misc");
  Sim.initialise();
  while (Sim.runSimulationStep()) {
  }

  // Grab the output plugins
  dynamo::OPMisc &opMisc = *Sim.getOutputPlugin<dynamo::OPMisc>();

  // Check the mean free time is roughly what is expected
  double MFT = opMisc.getMFT();
  BOOST_CHECK_CLOSE(MFT, 3, 0.000001);
  BOOST_CHECK_MESSAGE(
      Sim.checkSystem() <= 1,
      "There are more than one invalid states in the final configuration");
}

BOOST_AUTO_TEST_CASE(Dumb_Scheduler_CBT_Sorter) {
  runTest<dynamo::SDumb, dynamo::FELCBT>();
}

BOOST_AUTO_TEST_CASE(Dumb_Scheduler_BoundedPQ_Sorter) {
  runTest<dynamo::SDumb, dynamo::FELBoundedPQ<dynamo::PELMinMax<3>>>();
}

BOOST_AUTO_TEST_CASE(Neighbourlist_Scheduler_CBT_Sorter) {
  runTest<dynamo::SNeighbourList, dynamo::FELCBT>();
}

BOOST_AUTO_TEST_CASE(Neighbourlist_Scheduler_BoundedPQ_Sorter) {
  runTest<dynamo::SNeighbourList, dynamo::FELBoundedPQ<dynamo::PELMinMax<3>>>();
}
