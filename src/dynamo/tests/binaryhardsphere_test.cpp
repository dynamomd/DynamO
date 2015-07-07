#define BOOST_TEST_MODULE BinaryHardSphere_test
#include <boost/test/included/unit_test.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/BC/include.hpp>
#include <dynamo/ranges/include.hpp>
#include <dynamo/ranges/IDRangeRange.hpp>
#include <dynamo/inputplugins/cells/include.hpp>
#include <dynamo/species/point.hpp>
#include <dynamo/dynamics/newtonian.hpp>
#include <dynamo/schedulers/include.hpp>
#include <dynamo/schedulers/sorters/boundedPQFEL.hpp>
#include <dynamo/schedulers/sorters/MinMaxPEL.hpp>
#include <dynamo/inputplugins/include.hpp>
#include <dynamo/inputplugins/compression.hpp>
#include <dynamo/interactions/hardsphere.hpp>
#include <dynamo/outputplugins/misc.hpp>
#include <dynamo/outputplugins/msd.hpp>
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

  double massFrac = 0.001, sizeRatio = 0.5;
  size_t Na=100;
  
  Sim.dynamics = dynamo::shared_ptr<dynamo::Dynamics>(new dynamo::DynNewtonian(&Sim));
  Sim.BCs = dynamo::shared_ptr<dynamo::BoundaryCondition>(new dynamo::BCPeriodic(&Sim));
  Sim.ptrScheduler = dynamo::shared_ptr<dynamo::SNeighbourList>(new dynamo::SNeighbourList(&Sim, new DefaultSorter()));

  std::unique_ptr<dynamo::UCell> packptr(new dynamo::CUFCC(std::array<long, 3>{{10, 10, 10}}, dynamo::Vector{1, 1, 1}, new dynamo::UParticle()));
  packptr->initialise();
  std::vector<dynamo::Vector> latticeSites(packptr->placeObjects(dynamo::Vector{0,0,0}));
  Sim.primaryCellSize = dynamo::Vector{1,1,1};

  double simVol = 1.0;
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    simVol *= Sim.primaryCellSize[iDim];

  double particleDiam = std::cbrt(simVol * density / latticeSites.size());
  
  Sim.interactions.push_back(dynamo::shared_ptr<dynamo::Interaction>(new dynamo::IHardSphere(&Sim, particleDiam, new dynamo::IDPairRangeSingle(new dynamo::IDRangeRange(0, Na - 1)), "AAInt")));
  
  Sim.interactions.push_back(dynamo::shared_ptr<dynamo::Interaction>(new dynamo::IHardSphere(&Sim, ((1.0 + sizeRatio) / 2.0) * particleDiam, new dynamo::IDPairRangePair(new dynamo::IDRangeRange(0, Na - 1), new dynamo::IDRangeRange(Na, latticeSites.size() - 1)), "ABInt")));

  Sim.interactions.push_back(dynamo::shared_ptr<dynamo::Interaction>(new dynamo::IHardSphere(&Sim, sizeRatio * particleDiam, new dynamo::IDPairRangeAll(), "BBInt")));

  Sim.addSpecies(dynamo::shared_ptr<dynamo::Species>(new dynamo::SpPoint(&Sim, new dynamo::IDRangeRange(0, Na - 1), 1.0, "A", 0)));

  Sim.addSpecies(dynamo::shared_ptr<dynamo::Species>(new dynamo::SpPoint(&Sim, new dynamo::IDRangeRange(Na, latticeSites.size() - 1), massFrac, "B", 0)));

  Sim.units.setUnitLength(particleDiam);

  unsigned long nParticles = 0;
  Sim.particles.reserve(latticeSites.size());

  for (const dynamo::Vector & position : latticeSites)
    Sim.particles.push_back(dynamo::Particle(position, getRandVelVec() * Sim.units.unitVelocity(), nParticles++));

  Sim.ensemble = dynamo::Ensemble::loadEnsemble(Sim);

  dynamo::InputPlugin(&Sim, "Rescaler").zeroMomentum();
  dynamo::InputPlugin(&Sim, "Rescaler").rescaleVels(1.0);

  BOOST_CHECK_EQUAL(Sim.N(), 4000);
  //BOOST_CHECK_CLOSE(Sim.getNumberDensity() * Sim.units.unitVolume(), density, 0.000000001);
  //BOOST_CHECK_CLOSE(Sim.getPackingFraction(), Sim.getNumberDensity() * Sim.units.unitVolume() * M_PI / 6.0, 0.000000001);
}

BOOST_AUTO_TEST_CASE( Equilibrium_Simulation )
{
  {
    dynamo::Simulation Sim;
    init(Sim, 1.4);
    Sim.writeXMLfile("BHSequil.xml");
  }

  dynamo::Simulation Sim;
  Sim.loadXMLfile("BHSequil.xml");

  Sim.endEventCount = 1000000;
  Sim.addOutputPlugin("Misc");
  Sim.initialise();
  while (Sim.runSimulationStep()) {}

  Sim.reset();
  Sim.endEventCount = 1000000;
  Sim.addOutputPlugin("Misc"); 
  //Sim.addOutputPlugin("MSD");
  Sim.initialise();
  while (Sim.runSimulationStep()) {}
  

  //Taken from Lue 2005 DOI:10.1063/1.1834498
  const double expectedMFT = 0.0098213311089127;
  //const double expectedD = 0.247;

  //Grab the output plugins
  dynamo::OPMisc& opMisc = *Sim.getOutputPlugin<dynamo::OPMisc>();
  //dynamo::OPMSD& opMSD = *Sim.getOutputPlugin<dynamo::OPMSD>();

  //Check the mean free time is roughly what is expected
  double MFT = opMisc.getMFT();
  BOOST_CHECK_CLOSE(MFT, expectedMFT, 1);

  //Check the diffusion coefficient as well
  //double D = opMSD.calcD(*Sim.species[0]->getRange()) / Sim.units.unitDiffusion();
  //BOOST_CHECK_CLOSE(D, expectedD, 6);

  //Check the temperature is constant at 1
  std::cerr << "Temperature is " << opMisc.getCurrentkT() << std::endl;
  std::cerr << "Unit Temperature is " << Sim.units.unitEnergy() << std::endl;
  const double Temperature = opMisc.getCurrentkT() / Sim.units.unitEnergy();
  BOOST_CHECK_CLOSE(Temperature, 1.0, 0.000000001);

  //Check that the momentum is around 0
  dynamo::Vector momentum = Sim.getOutputPlugin<dynamo::OPMisc>()->getCurrentMomentum();
  BOOST_CHECK_SMALL(momentum.nrm() / Sim.units.unitMomentum(), 0.0000000001);
  
  BOOST_CHECK_MESSAGE(Sim.checkSystem() <= 1, "There are more than two invalid states in the final configuration");
}

//BOOST_AUTO_TEST_CASE( Compression_Simulation )
//{
//  dynamo::Simulation Sim;
//  init(Sim, 0.1);
//
//  const double growthRate = 1;
//  const double targetDensity = 0.9;
//  Sim.endEventCount = 1000000;
//  Sim.addOutputPlugin("Misc");
//
//  dynamo::shared_ptr<dynamo::IPCompression> compressPlug(new dynamo::IPCompression(&Sim, growthRate));
//  compressPlug->MakeGrowth();
//  
//  compressPlug->limitDensity(targetDensity);
//  
//  //Not needed in this system
//  //compressPlug->CellSchedulerHack();
//
//  Sim.initialise();
//  while (Sim.runSimulationStep()) {}
//  compressPlug->RestoreSystem();
//
//  BOOST_CHECK_CLOSE(Sim.getNumberDensity() * Sim.units.unitVolume(), targetDensity, 0.000000001);
//  BOOST_CHECK_CLOSE(Sim.getPackingFraction(), Sim.getNumberDensity() * Sim.units.unitVolume() * M_PI / 6.0, 0.000000001);
//
//  BOOST_CHECK_MESSAGE(Sim.checkSystem() <= 1, "After compression, there are more than one invalid states in the final configuration");
//}
