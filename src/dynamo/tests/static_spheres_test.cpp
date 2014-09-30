#define BOOST_TEST_MODULE Static_Spheres_test
#include <boost/test/included/unit_test.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/BC/include.hpp>
#include <dynamo/ranges/include.hpp>
#include <dynamo/inputplugins/cells/include.hpp>
#include <dynamo/species/point.hpp>
#include <dynamo/species/fixedCollider.hpp>
#include <dynamo/dynamics/gravity.hpp>
#include <dynamo/schedulers/include.hpp>
#include <dynamo/schedulers/sorters/include.hpp>
#include <dynamo/locals/lwall.hpp>
#include <dynamo/inputplugins/include.hpp>
#include <dynamo/inputplugins/compression.hpp>
#include <dynamo/interactions/hardsphere.hpp>
#include <dynamo/interactions/nullInteraction.hpp>
#include <dynamo/outputplugins/misc.hpp>
#include <dynamo/outputplugins/msd.hpp>
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

void init(dynamo::Simulation& Sim, const double density)
{
  RNG.seed(std::random_device()());
  Sim.ranGenerator.seed(std::random_device()());
  Sim.dynamics = dynamo::shared_ptr<dynamo::Dynamics>(new dynamo::DynGravity(&Sim, dynamo::Vector(0,-1,0)));
  Sim.BCs = dynamo::shared_ptr<dynamo::BoundaryCondition>(new dynamo::BCNone(&Sim));
  Sim.ptrScheduler = dynamo::shared_ptr<dynamo::SNeighbourList>(new dynamo::SNeighbourList(&Sim, new dynamo::FELCBT()));
  Sim.primaryCellSize = dynamo::Vector(52,52,52);

  Sim.addSpecies(dynamo::shared_ptr<dynamo::Species>(new dynamo::SpPoint(&Sim, new dynamo::IDRangeRange(0, 0), 1.0, "Bulk", 0)));
  Sim.addSpecies(dynamo::shared_ptr<dynamo::Species>(new dynamo::SpFixedCollider(&Sim, new dynamo::IDRangeRange(1, 8), "FixedColliders", 1)));

  Sim.interactions.push_back(dynamo::shared_ptr<dynamo::Interaction>(new dynamo::IHardSphere(&Sim, 1.0, 1, new dynamo::IDPairRangePair(new dynamo::IDRangeAll(&Sim), new dynamo::IDRangeRange(0,0)), "Bulk")));

  Sim.interactions.push_back(dynamo::shared_ptr<dynamo::Interaction>(new dynamo::INull(&Sim, new dynamo::IDPairRangeAll(), "NoInteraction")));
  
  Sim.locals.push_back(dynamo::shared_ptr<dynamo::Local>(new dynamo::LWall(&Sim, 1.0, 1.0, dynamo::Vector(0,1,0), dynamo::Vector(0,-2.67753263802375e+01,0), "GroundPlate", new dynamo::IDRangeAll(&Sim))));
  
  Sim.particles.push_back(dynamo::Particle(dynamo::Vector(0, 4, 0), dynamo::Vector(0,0,0), Sim.particles.size()));
  Sim.particles.push_back(dynamo::Particle(dynamo::Vector(0.6, 1, 0), dynamo::Vector(0,0,0), Sim.particles.size()));
  Sim.particles.push_back(dynamo::Particle(dynamo::Vector(-1.51, 1, 0), dynamo::Vector(0,0,0), Sim.particles.size()));
  Sim.particles.push_back(dynamo::Particle(dynamo::Vector(-2.51, 1.5, 0), dynamo::Vector(0,0,0), Sim.particles.size()));
  Sim.particles.push_back(dynamo::Particle(dynamo::Vector(-3.51, 2, 0), dynamo::Vector(0,0,0), Sim.particles.size()));
  Sim.particles.push_back(dynamo::Particle(dynamo::Vector(-3.51, 3.5, 0), dynamo::Vector(0,0,0), Sim.particles.size()));
  Sim.particles.push_back(dynamo::Particle(dynamo::Vector(1.6, 2, 0), dynamo::Vector(0,0,0), Sim.particles.size()));
  Sim.particles.push_back(dynamo::Particle(dynamo::Vector(2, 3.5, 0), dynamo::Vector(0,0,0), Sim.particles.size()));
  Sim.particles.push_back(dynamo::Particle(dynamo::Vector(-0.75, 0.5, 0), dynamo::Vector(0,0,0), Sim.particles.size()));

  Sim.ensemble = dynamo::Ensemble::loadEnsemble(Sim);
}

BOOST_AUTO_TEST_CASE( Test_Simulation )
{
  {
    dynamo::Simulation Sim;
    init(Sim, 0.5);
    Sim.initialise();
    Sim.writeXMLfile("staticsphere.xml");
  }

  dynamo::Simulation Sim;
  Sim.loadXMLfile("staticsphere.xml");

  Sim.endEventCount = 500000;
  Sim.addOutputPlugin("Misc");
  Sim.initialise();
  while (Sim.runSimulationStep()) {}

  const double expectedMFT = 7.81945252098576;
  dynamo::OPMisc& opMisc = *Sim.getOutputPlugin<dynamo::OPMisc>();
  //Check the mean free time is roughly what is expected
  double MFT = opMisc.getMFT();
  BOOST_CHECK_CLOSE(MFT, expectedMFT, 0.1);
  BOOST_CHECK_MESSAGE(Sim.checkSystem() <= 1, "There are more than two invalid states in the final configuration");
}
