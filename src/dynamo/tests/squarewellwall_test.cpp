#define BOOST_TEST_MODULE SquareWellWall_test
#include <boost/test/included/unit_test.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/BC/include.hpp>
#include <dynamo/ranges/include.hpp>
#include <dynamo/inputplugins/cells/include.hpp>
#include <dynamo/species/point.hpp>
#include <dynamo/dynamics/newtonian.hpp>
#include <dynamo/schedulers/include.hpp>
#include <dynamo/schedulers/sorters/include.hpp>
#include <dynamo/inputplugins/include.hpp>
#include <dynamo/inputplugins/compression.hpp>
#include <dynamo/interactions/squarewell.hpp>
#include <dynamo/systems/andersenThermostat.hpp>
#include <dynamo/outputplugins/misc.hpp>
#include <dynamo/locals/lwall.hpp>
#include <random>

std::mt19937 RNG;

dynamo::Vector getRandVelVec()
{
  //See http://mathworld.wolfram.com/SpherePointPicking.html
  std::normal_distribution<> normal_dist(0.0, (1.0 / sqrt(double(NDIM))));
  
  dynamo::Vector tmpVec;
  for (size_t iDim = 0; iDim < NDIM; iDim++)
    tmpVec[iDim] = normal_dist(RNG);
  
  return tmpVec;
}

void init(dynamo::Simulation& Sim)
{
  RNG.seed(std::random_device()());
  Sim.ranGenerator.seed(std::random_device()());

  Sim.primaryCellSize = dynamo::Vector{6.1, 10, 10};

  Sim.dynamics = dynamo::shared_ptr<dynamo::Dynamics>(new dynamo::DynNewtonian(&Sim));
  Sim.BCs = dynamo::shared_ptr<dynamo::BoundaryCondition>(new dynamo::BCPeriodicExceptX(&Sim));
  Sim.ptrScheduler = dynamo::shared_ptr<dynamo::SNeighbourList>(new dynamo::SNeighbourList(&Sim, new dynamo::FELCBT()));
  Sim.addSpecies(dynamo::shared_ptr<dynamo::Species>(new dynamo::SpPoint(&Sim, new dynamo::IDRangeAll(&Sim), 1.0, "Bulk", 0)));

  Sim.locals.push_back(dynamo::shared_ptr<dynamo::Local>(new dynamo::LWall(&Sim, 1.0, 1.0, dynamo::Vector{1, 0, 0}, dynamo::Vector{-3, 0, 0}, "LowWall", new dynamo::IDRangeAll(&Sim))));

  Sim.locals.push_back(dynamo::shared_ptr<dynamo::Local>(new dynamo::LWall(&Sim, 1.0, 1.0, dynamo::Vector{-1, 0, 0}, dynamo::Vector{3, 0, 0}, "HighWall", new dynamo::IDRangeAll(&Sim))));

  Sim.interactions.push_back(dynamo::shared_ptr<dynamo::Interaction>(new dynamo::ISquareWell(&Sim, 1, 1.5, 1, 1, new dynamo::IDPairRangeAll(), "Bulk")));

  Sim.particles.push_back(dynamo::Particle(dynamo::Vector{0.75, 0, 0}, dynamo::Vector{2,0,0}, Sim.particles.size()));
  Sim.particles.push_back(dynamo::Particle(dynamo::Vector{-0.75, 0, 0}, dynamo::Vector{-2,0,0}, Sim.particles.size()));

  Sim.ensemble = dynamo::Ensemble::loadEnsemble(Sim);

  BOOST_CHECK_EQUAL(Sim.N(), 2);
}

BOOST_AUTO_TEST_CASE( NVE_Simulation )
{
  {
    dynamo::Simulation Sim;
    init(Sim);
    Sim.writeXMLfile("SquareWellWall.xml");
  }

  dynamo::Simulation Sim;
  Sim.loadXMLfile("SquareWellWall.xml");

  //System is periodic over 5 collisions
  Sim.endEventCount = 9995;
  Sim.addOutputPlugin("Misc");
  Sim.initialise();
  const double totalEinit = Sim.getOutputPlugin<dynamo::OPMisc>()->getTotalEnergy();

  while (Sim.runSimulationStep()) {}

  const double totalEend = Sim.getOutputPlugin<dynamo::OPMisc>()->getTotalEnergy();
  dynamo::Vector momentum = Sim.getOutputPlugin<dynamo::OPMisc>()->getCurrentMomentum();

  BOOST_CHECK_CLOSE(totalEinit, totalEend, 0.000000001);
  //Check the final positions are close to the initial (system has a
  //period of 5 collisions)
  BOOST_CHECK_CLOSE(Sim.particles[0].getPosition()[0], 0.75, 0.000000001);
  BOOST_CHECK_CLOSE(Sim.particles[0].getPosition()[1], 0, 0.000000001);
  BOOST_CHECK_CLOSE(Sim.particles[0].getPosition()[2], 0, 0.000000001);
  BOOST_CHECK_CLOSE(Sim.particles[1].getPosition()[0], -0.75, 0.000000001);
  BOOST_CHECK_CLOSE(Sim.particles[1].getPosition()[1], 0, 0.000000001);
  BOOST_CHECK_CLOSE(Sim.particles[1].getPosition()[2], 0, 0.000000001);
  BOOST_CHECK_SMALL(momentum.nrm() / Sim.units.unitMomentum(), 0.0000000001);
  BOOST_CHECK_MESSAGE(Sim.checkSystem() <= 2, "There are more than two invalid states in the final configuration");
}

BOOST_AUTO_TEST_CASE( Null_Compression_Simulation )
{
  {
    dynamo::Simulation Sim;
    init(Sim);
    Sim.writeXMLfile("SquareWellWall.xml");
  }

  dynamo::Simulation Sim;
  Sim.loadXMLfile("SquareWellWall.xml");

  dynamo::shared_ptr<dynamo::IPCompression> compressPlug(new dynamo::IPCompression(&Sim, 0.0));
  compressPlug->MakeGrowth();
  Sim.endEventCount = 9995;
  Sim.addOutputPlugin("Misc");
  Sim.initialise();
  const double totalEinit = Sim.getOutputPlugin<dynamo::OPMisc>()->getTotalEnergy();

  while (Sim.runSimulationStep()) {}

  compressPlug->RestoreSystem();

  const double totalEend = Sim.getOutputPlugin<dynamo::OPMisc>()->getTotalEnergy();
  dynamo::Vector momentum = Sim.getOutputPlugin<dynamo::OPMisc>()->getCurrentMomentum();

  BOOST_CHECK_CLOSE(totalEinit, totalEend, 0.000000001);
  //Check the final positions are close to the initial (system has a
  //period of 5 collisions)
  BOOST_CHECK_CLOSE(Sim.particles[0].getPosition()[0], 0.75, 0.000000001);
  BOOST_CHECK_CLOSE(Sim.particles[0].getPosition()[1], 0, 0.000000001);
  BOOST_CHECK_CLOSE(Sim.particles[0].getPosition()[2], 0, 0.000000001);
  BOOST_CHECK_CLOSE(Sim.particles[1].getPosition()[0], -0.75, 0.000000001);
  BOOST_CHECK_CLOSE(Sim.particles[1].getPosition()[1], 0, 0.000000001);
  BOOST_CHECK_CLOSE(Sim.particles[1].getPosition()[2], 0, 0.000000001);
  BOOST_CHECK_SMALL(momentum.nrm() / Sim.units.unitMomentum(), 0.0000000001);
  BOOST_CHECK_MESSAGE(Sim.checkSystem() <= 2, "There are more than two invalid states in the final configuration");
}
