#define BOOST_TEST_MODULE Hardsphere_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <dynamo/simulation.hpp>
#include <dynamo/BC/include.hpp>
#include <dynamo/ranges/include.hpp>
#include <dynamo/inputplugins/cells/include.hpp>
#include <dynamo/species/point.hpp>
#include <dynamo/dynamics/newtonian.hpp>
#include <dynamo/globals/cells.hpp>
#include <dynamo/globals/cellsShearing.hpp>
#include <dynamo/schedulers/include.hpp>
#include <dynamo/schedulers/sorters/include.hpp>
#include <dynamo/inputplugins/include.hpp>
#include <dynamo/interactions/squarewell.hpp>
#include <dynamo/systems/andersenThermostat.hpp>
#include <dynamo/outputplugins/misc.hpp>
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

void init(dynamo::Simulation& Sim)
{
  RNG.seed(std::random_device()());
  Sim.ranGenerator.seed(std::random_device()());

  const double density = 0.5;
  const double elasticity = 1.0;
  const double lambda = 1.5;
  const double welldepth = 1.0;

  Sim.dynamics = dynamo::shared_ptr<dynamo::Dynamics>(new dynamo::DynNewtonian(&Sim));
  Sim.BCs = dynamo::shared_ptr<dynamo::BoundaryCondition>(new dynamo::BCPeriodic(&Sim));
  Sim.ptrScheduler = dynamo::shared_ptr<dynamo::SNeighbourList>(new dynamo::SNeighbourList(&Sim, new DefaultSorter()));

  std::unique_ptr<dynamo::UCell> packptr(new dynamo::CUFCC(std::array<long, 3>{{7,7,7}}, dynamo::Vector(1,1,1), new dynamo::UParticle()));
  packptr->initialise();
  std::vector<dynamo::Vector> latticeSites(packptr->placeObjects(dynamo::Vector(0,0,0)));
  Sim.primaryCellSize = dynamo::Vector(1,1,1);

  double simVol = 1.0;
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    simVol *= Sim.primaryCellSize[iDim];

  double particleDiam = std::cbrt(simVol * density / latticeSites.size());

  Sim.interactions.push_back(dynamo::shared_ptr<dynamo::Interaction>(new dynamo::ISquareWell(&Sim, particleDiam, lambda, welldepth, elasticity, new dynamo::IDPairRangeAll(), "Bulk")));
  Sim.addSpecies(dynamo::shared_ptr<dynamo::Species>(new dynamo::SpPoint(&Sim, new dynamo::IDRangeAll(&Sim), 1.0, "Bulk", 0, "Bulk")));
  Sim.units.setUnitLength(particleDiam);
  Sim.units.setUnitTime(particleDiam); 

  unsigned long nParticles = 0;
  Sim.particles.reserve(latticeSites.size());
  for (const dynamo::Vector & position : latticeSites)
    Sim.particles.push_back(dynamo::Particle(position, getRandVelVec() * Sim.units.unitVelocity(), nParticles++));

  Sim.globals.push_back(dynamo::shared_ptr<dynamo::Global>(new dynamo::GCells(&Sim,"SchedulerNBList")));

  Sim.ensemble = dynamo::Ensemble::loadEnsemble(Sim);

  dynamo::InputPlugin(&Sim, "Rescaler").zeroMomentum();
  dynamo::InputPlugin(&Sim, "Rescaler").rescaleVels(1.0);

  BOOST_CHECK_EQUAL(Sim.N(), 1372);
  BOOST_CHECK_CLOSE(Sim.getNumberDensity() * Sim.units.unitVolume(), 0.5, 0.000000001);
  BOOST_CHECK_CLOSE(Sim.getPackingFraction(), Sim.getNumberDensity() * Sim.units.unitVolume() * M_PI / 6.0, 0.000000001);
}

BOOST_AUTO_TEST_CASE( NVE_Simulation )
{
  dynamo::Simulation Sim;
  init(Sim);

  Sim.status = dynamo::CONFIG_LOADED;
  Sim.endEventCount = 100000;
  Sim.addOutputPlugin("Misc");
  Sim.initialise();
  const double totalEinit = Sim.getOutputPlugin<dynamo::OPMisc>()->getTotalEnergy();
  while (Sim.runSimulationStep()) {}
  const double totalEequil = Sim.getOutputPlugin<dynamo::OPMisc>()->getTotalEnergy();
  BOOST_CHECK_CLOSE(totalEinit, totalEequil, 0.000000001);

  Sim.endEventCount += 100000;
  Sim.reset();
  Sim.addOutputPlugin("Misc");
  Sim.initialise();
  const double totalEprerun = Sim.getOutputPlugin<dynamo::OPMisc>()->getTotalEnergy();
  BOOST_CHECK_CLOSE(totalEequil, totalEprerun, 0.000000001);
  while (Sim.runSimulationStep()) {}
  const double totalEfinal = Sim.getOutputPlugin<dynamo::OPMisc>()->getTotalEnergy();
  BOOST_CHECK_CLOSE(totalEprerun, totalEfinal, 0.000000001);
  
  //Check that the momentum is around 0
  dynamo::Vector momentum = Sim.getOutputPlugin<dynamo::OPMisc>()->getCurrentMomentum();
  BOOST_CHECK_SMALL(momentum.nrm() / Sim.units.unitMomentum(), 0.0000000001);
}

BOOST_AUTO_TEST_CASE( NVT_Simulation )
{
  dynamo::Simulation Sim;
  init(Sim);

  Sim.systems.push_back(dynamo::shared_ptr<dynamo::System>(new dynamo::SysAndersen(&Sim, 0.036 / Sim.N(), 1.0 * Sim.units.unitEnergy(), "Thermostat")));
  Sim.ensemble = dynamo::Ensemble::loadEnsemble(Sim);

  Sim.status = dynamo::CONFIG_LOADED;
  Sim.eventPrintInterval = 50000;
  Sim.endEventCount = 300000;
  Sim.addOutputPlugin("Misc");
  Sim.initialise();
  while (Sim.runSimulationStep()) {}

  Sim.reset();
  Sim.endEventCount = 100000;
  Sim.addOutputPlugin("Misc");
  Sim.initialise();

  while (Sim.runSimulationStep()) {}

  //Check the mean free time is roughly what is expected
  double MFT = Sim.getOutputPlugin<dynamo::OPMisc>()->getMFT();
  BOOST_CHECK_CLOSE(MFT, 0.0368185, 5);

  //Check the temperature is constant at 1
  double Temperature = Sim.getOutputPlugin<dynamo::OPMisc>()->getCurrentkT() / Sim.units.unitEnergy();
  BOOST_CHECK_CLOSE(Temperature, 1.0, 5);
}
