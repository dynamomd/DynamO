#define BOOST_TEST_MODULE Infinite_Mass_Spheres_test
#include <boost/test/included/unit_test.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/BC/include.hpp>
#include <dynamo/ranges/include.hpp>
#include <dynamo/inputplugins/cells/include.hpp>
#include <dynamo/species/point.hpp>
#include <dynamo/species/fixedCollider.hpp>
#include <dynamo/dynamics/newtonian.hpp>
#include <dynamo/schedulers/include.hpp>
#include <dynamo/schedulers/sorters/include.hpp>
#include <dynamo/inputplugins/include.hpp>
#include <dynamo/inputplugins/compression.hpp>
#include <dynamo/interactions/hardsphere.hpp>
#include <dynamo/outputplugins/misc.hpp>
#include <dynamo/outputplugins/msd.hpp>
#include <random>

std::mt19937 RNG;
typedef dynamo::FELBoundedPQ<dynamo::PELMinMax<3> > DefaultSorter;

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

  double L = std::cbrt(1372 / density);

  Sim.dynamics = dynamo::shared_ptr<dynamo::Dynamics>(new dynamo::DynNewtonian(&Sim));
  Sim.BCs = dynamo::shared_ptr<dynamo::BoundaryCondition>(new dynamo::BCPeriodic(&Sim));
  Sim.ptrScheduler = dynamo::shared_ptr<dynamo::SNeighbourList>(new dynamo::SNeighbourList(&Sim, new DefaultSorter()));

  std::unique_ptr<dynamo::UCell> packptr(new dynamo::CUFCC(std::array<long, 3>{{7,7,7}}, dynamo::Vector{L,L,L}, new dynamo::UParticle()));

  packptr->initialise();
  std::vector<dynamo::Vector> latticeSites(packptr->placeObjects(dynamo::Vector{0,0,0}));
  Sim.primaryCellSize = dynamo::Vector{L,L,L};

  Sim.addSpecies(dynamo::shared_ptr<dynamo::Species>(new dynamo::SpFixedCollider(&Sim, new dynamo::IDRangeRange(0,1), "HeavySpheres", 0)));

  Sim.addSpecies(dynamo::shared_ptr<dynamo::Species>(new dynamo::SpPoint(&Sim, new dynamo::IDRangeRange(2,1371), 1.0, "Bulk", 1)));

  Sim.interactions.push_back(dynamo::shared_ptr<dynamo::Interaction>(new dynamo::IHardSphere(&Sim, 1.0, 1, new dynamo::IDPairRangeSingle(new dynamo::IDRangeRange(0,1)), "HeavyHeavy")));

  Sim.interactions.push_back(dynamo::shared_ptr<dynamo::Interaction>(new dynamo::IHardSphere(&Sim, 0.6, 1, new dynamo::IDPairRangePair(new dynamo::IDRangeRange(0,1), new dynamo::IDRangeRange(2,1371)), "HeavyBulk")));

  Sim.interactions.push_back(dynamo::shared_ptr<dynamo::Interaction>(new dynamo::IHardSphere(&Sim, 0.2, 1, new dynamo::IDPairRangeAll(), "BulkBulk")));

  unsigned long nParticles = 0;
  Sim.particles.reserve(latticeSites.size());
  for (const dynamo::Vector & position : latticeSites)
    Sim.particles.push_back(dynamo::Particle(position, getRandVelVec() * Sim.units.unitVelocity(), nParticles++));

  Sim.particles[0].getVelocity() = dynamo::Vector{0,0,0};
  Sim.particles[1].getVelocity() = dynamo::Vector{-1,-1,0};
  Sim.ensemble = dynamo::Ensemble::loadEnsemble(Sim);

  dynamo::InputPlugin(&Sim, "Rescaler").zeroMomentum();
  dynamo::InputPlugin(&Sim, "Rescaler").rescaleVels(1.0);

  BOOST_CHECK_EQUAL(Sim.N(), 1372);
  BOOST_CHECK_CLOSE(Sim.getNumberDensity() * Sim.units.unitVolume(), density, 0.000000001);
}

BOOST_AUTO_TEST_CASE( Equilibrium_Simulation )
{
  {
    dynamo::Simulation Sim;
    init(Sim, 0.5);
    Sim.initialise();
    Sim.writeXMLfile("infmass.xml");
  }

  dynamo::Simulation Sim;
  Sim.loadXMLfile("infmass.xml");

  Sim.endEventCount = 100000;
  Sim.addOutputPlugin("Misc");
  Sim.initialise();
  while (Sim.runSimulationStep()) {}

  const double expectedMFT = 5.74807417926229;
  dynamo::OPMisc& opMisc = *Sim.getOutputPlugin<dynamo::OPMisc>();
  //Check the mean free time is roughly what is expected
  double MFT = opMisc.getMFT();
  BOOST_CHECK_CLOSE(MFT, expectedMFT, 2);
  //6.2715616026913548
  BOOST_CHECK_MESSAGE(Sim.checkSystem() <= 1, "There are more than two invalid states in the final configuration");
}

BOOST_AUTO_TEST_CASE( Compression_Simulation )
{
  dynamo::Simulation Sim;
  init(Sim, 0.1);

  const double growthRate = 0.1;
  const double targetDensity = 0.9;
  Sim.endEventCount = 1000000;
  Sim.addOutputPlugin("Misc");

  dynamo::shared_ptr<dynamo::IPCompression> compressPlug(new dynamo::IPCompression(&Sim, growthRate));
  compressPlug->MakeGrowth();
  
  compressPlug->limitDensity(targetDensity);
  
  //Not needed in this system
  //compressPlug->CellSchedulerHack();

  Sim.initialise();
  while (Sim.runSimulationStep()) {}
  compressPlug->RestoreSystem();
  
  BOOST_CHECK_CLOSE(Sim.getNumberDensity() * Sim.units.unitVolume(), targetDensity, 0.000000001);
  //BOOST_CHECK_CLOSE(Sim.getPackingFraction(), Sim.getNumberDensity() * Sim.units.unitVolume() * M_PI / 6.0, 0.000000001);

  BOOST_CHECK_MESSAGE(Sim.checkSystem() <= 1, "After compression, there are more than one invalid states in the final configuration");
}
