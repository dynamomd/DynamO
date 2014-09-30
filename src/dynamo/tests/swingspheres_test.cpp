#define BOOST_TEST_MODULE SwingSpheres_test
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
#include <dynamo/inputplugins/include.hpp>
#include <dynamo/inputplugins/compression.hpp>
#include <dynamo/interactions/hardsphere.hpp>
#include <dynamo/interactions/squarebond.hpp>
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
  size_t N = 11;
  double bond_elasticity = 0.9;

  RNG.seed(std::random_device()());
  Sim.ranGenerator.seed(std::random_device()());

  Sim.primaryCellSize = dynamo::Vector(60, 60, 60);

  Sim.dynamics = dynamo::shared_ptr<dynamo::Dynamics>(new dynamo::DynGravity(&Sim, dynamo::Vector(0, -1, 0), 0, 0.01));
  Sim.BCs = dynamo::shared_ptr<dynamo::BoundaryCondition>(new dynamo::BCNone(&Sim));
  Sim.ptrScheduler = dynamo::shared_ptr<dynamo::SNeighbourList>(new dynamo::SNeighbourList(&Sim, new dynamo::FELCBT()));

  dynamo::shared_ptr<dynamo::ParticleProperty> D(new dynamo::ParticleProperty(N, dynamo::Property::Units::Length(), "D", 1.0));
  dynamo::shared_ptr<dynamo::ParticleProperty> M(new dynamo::ParticleProperty(N, dynamo::Property::Units::Mass(), "M", 1.0));
  Sim._properties.push(D);
  Sim._properties.push(M);

  Sim.addSpecies(dynamo::shared_ptr<dynamo::Species>(new dynamo::SpFixedCollider(&Sim, new dynamo::IDRangeRange(0,0), "FixedColliders", 0)));
  Sim.addSpecies(dynamo::shared_ptr<dynamo::Species>(new dynamo::SpPoint(&Sim, new dynamo::IDRangeRange(1, N - 1), "M", "Bulk", 1)));

  Sim.interactions.push_back(dynamo::shared_ptr<dynamo::Interaction>(new dynamo::ISquareBond(&Sim, "D", 1.06, bond_elasticity, new dynamo::IDPairRangeChains(0, N - 1, N), "Bonds")));
  Sim.interactions.push_back(dynamo::shared_ptr<dynamo::Interaction>(new dynamo::IHardSphere(&Sim, "D", 1.0, new dynamo::IDPairRangeAll(), "Bulk")));
  
  for (size_t i(0); i < N; ++i)
    {
      Sim.particles.push_back(dynamo::Particle(dynamo::Vector(i * 1.05, 0, 0), dynamo::Vector(0, 0, 0), Sim.particles.size()));
      D->getProperty(i) = 1;
      M->getProperty(i) = 1;
    }

  D->getProperty(N-1) = 2;
  M->getProperty(N-1) = 100;
  

  Sim.ensemble = dynamo::Ensemble::loadEnsemble(Sim);
}

BOOST_AUTO_TEST_CASE( Simulation )
{
  //This tests properties, the tc model gravity bonds, and static
  //spheres/bonds in gravity.
  {
    dynamo::Simulation Sim;
    init(Sim);
    Sim.writeXMLfile("SwingSpheres.xml");
  }

  dynamo::Simulation Sim;
  Sim.loadXMLfile("SwingSpheres.xml");

  Sim.endEventCount = 500000;
  Sim.initialise();
  while (Sim.runSimulationStep()) {}

  BOOST_CHECK_MESSAGE(Sim.checkSystem() <= 2, "There are more than two invalid states in the final configuration");
}

