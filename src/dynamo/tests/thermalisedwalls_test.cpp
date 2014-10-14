#define BOOST_TEST_MODULE ThermalisedWalls_test
#include <boost/test/included/unit_test.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/BC/include.hpp>
#include <dynamo/ranges/include.hpp>
#include <dynamo/inputplugins/cells/include.hpp>
#include <dynamo/species/point.hpp>
#include <dynamo/dynamics/gravity.hpp>
#include <dynamo/schedulers/include.hpp>
#include <dynamo/schedulers/sorters/boundedPQFEL.hpp>
#include <dynamo/schedulers/sorters/MinMaxPEL.hpp>
#include <dynamo/inputplugins/include.hpp>
#include <dynamo/inputplugins/compression.hpp>
#include <dynamo/interactions/hardsphere.hpp>
#include <dynamo/outputplugins/misc.hpp>
#include <dynamo/locals/lwall.hpp>
#include <random>

std::mt19937 RNG;
typedef dynamo::BoundedPQFEL<dynamo::MinMaxPEL<3> > DefaultSorter;

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

  const double elasticity = 1.0;
  const size_t cells = 7;
  const double wallkT = 1.0;

  std::unique_ptr<dynamo::UCell> packptr(new dynamo::CUFCC(std::array<long, 3>{{cells, cells, cells}}, dynamo::Vector{1,1,1}, new dynamo::UParticle()));
  packptr->initialise();
  std::vector<dynamo::Vector> latticeSites(packptr->placeObjects(dynamo::Vector{0,0,0}));
  const size_t N = latticeSites.size();
  const double boxL = std::cbrt(N / density);
  Sim.primaryCellSize = dynamo::Vector{boxL, boxL, boxL};

  dynamo::shared_ptr<dynamo::ParticleProperty> D(new dynamo::ParticleProperty(N, dynamo::Property::Units::Length(), "D", 1.0));
  dynamo::shared_ptr<dynamo::ParticleProperty> M(new dynamo::ParticleProperty(N, dynamo::Property::Units::Mass(), "M", 1.0));
  Sim._properties.push(D);
  Sim._properties.push(M);

  Sim.dynamics = dynamo::shared_ptr<dynamo::Dynamics>(new dynamo::DynGravity(&Sim, dynamo::Vector{0, -1, 0}));
  Sim.BCs = dynamo::shared_ptr<dynamo::BoundaryCondition>(new dynamo::BCNone(&Sim));
  Sim.ptrScheduler = dynamo::shared_ptr<dynamo::SNeighbourList>(new dynamo::SNeighbourList(&Sim, new DefaultSorter()));

  Sim.interactions.push_back(dynamo::shared_ptr<dynamo::Interaction>(new dynamo::IHardSphere(&Sim, "D", elasticity, new dynamo::IDPairRangeAll(), "Bulk")));
  Sim.addSpecies(dynamo::shared_ptr<dynamo::Species>(new dynamo::SpPoint(&Sim, new dynamo::IDRangeAll(&Sim), "M", "Bulk", 0)));

  Sim.locals.push_back(dynamo::shared_ptr<dynamo::Local>(new dynamo::LWall(&Sim, elasticity, "D", dynamo::Vector{1, 0, 0}, dynamo::Vector{-0.5 * Sim.primaryCellSize[0] - 1, 0, 0}, "XwallLow", new dynamo::IDRangeAll(&Sim))));
  Sim.locals.push_back(dynamo::shared_ptr<dynamo::Local>(new dynamo::LWall(&Sim, elasticity, "D", dynamo::Vector{-1, 0, 0}, dynamo::Vector{0.5 * Sim.primaryCellSize[0] + 1, 0, 0}, "XwallHigh", new dynamo::IDRangeAll(&Sim))));

  Sim.locals.push_back(dynamo::shared_ptr<dynamo::Local>(new dynamo::LWall(&Sim, elasticity, "D", dynamo::Vector{0, 0, 1}, dynamo::Vector{0, 0, -0.5 * Sim.primaryCellSize[2] - 1}, "ZwallLow", new dynamo::IDRangeAll(&Sim))));
  Sim.locals.push_back(dynamo::shared_ptr<dynamo::Local>(new dynamo::LWall(&Sim, elasticity, "D", dynamo::Vector{0, 0, -1}, dynamo::Vector{0, 0, 0.5 * Sim.primaryCellSize[2] + 1}, "ZwallHigh", new dynamo::IDRangeAll(&Sim))));

  Sim.locals.push_back(dynamo::shared_ptr<dynamo::Local>(new dynamo::LWall(&Sim, elasticity, "D", dynamo::Vector{0, 1, 0}, dynamo::Vector{0, - 0.5 * Sim.primaryCellSize[1] - 1, 0}, "GroundPlate", new dynamo::IDRangeAll(&Sim), wallkT)));
  
  for (size_t i(0); i < latticeSites.size(); ++i)
    {
      Sim.particles.push_back(dynamo::Particle(boxL * latticeSites[i], getRandVelVec(), Sim.particles.size()));
      D->getProperty(i) = (i < 100) ? 1 : 0.5;
      M->getProperty(i) = (i < 100) ? 1 : 0.5 * 0.5 * 0.5;
    }

  Sim.ensemble = dynamo::Ensemble::loadEnsemble(Sim);

  dynamo::InputPlugin(&Sim, "Rescaler").zeroMomentum();
  dynamo::InputPlugin(&Sim, "Rescaler").rescaleVels(1.0);

  BOOST_CHECK_EQUAL(Sim.N(), 1372);
  BOOST_CHECK_CLOSE(Sim.getNumberDensity() * Sim.units.unitVolume(), density, 0.000000001);
}

BOOST_AUTO_TEST_CASE( Simulation )
{
  {
    dynamo::Simulation Sim;
    init(Sim, 0.1);
    Sim.writeXMLfile("ThermalisedPlate.xml");
  }

  dynamo::Simulation Sim;
  Sim.loadXMLfile("ThermalisedPlate.xml");

  Sim.endEventCount = 100000;
  Sim.addOutputPlugin("Misc");
  Sim.initialise();
  while (Sim.runSimulationStep()) {}

  Sim.reset();
  Sim.endEventCount = 400000;
  Sim.addOutputPlugin("Misc"); 
  Sim.initialise();
  while (Sim.runSimulationStep()) {}

  const double expectedMFT = 0.89464166687653845;

  //Grab the output plugins
  dynamo::OPMisc& opMisc = *Sim.getOutputPlugin<dynamo::OPMisc>();

  //Check the mean free time is roughly what is expected
  double MFT = opMisc.getMFT();
  BOOST_CHECK_CLOSE(MFT, expectedMFT, 6);

  BOOST_CHECK_MESSAGE(Sim.checkSystem() <= 1, "There are more than two invalid states in the final configuration");
}

