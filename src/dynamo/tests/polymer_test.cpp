#define BOOST_TEST_MODULE Polymer_test
#include <boost/test/included/unit_test.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/BC/include.hpp>
#include <dynamo/ranges/include.hpp>
#include <dynamo/inputplugins/cells/include.hpp>
#include <dynamo/species/point.hpp>
#include <dynamo/dynamics/newtonian.hpp>
#include <dynamo/schedulers/include.hpp>
#include <dynamo/schedulers/sorters/heapPEL.hpp>
#include <dynamo/schedulers/sorters/CBTFEL.hpp>
#include <dynamo/inputplugins/include.hpp>
#include <dynamo/inputplugins/compression.hpp>
#include <dynamo/interactions/squarebond.hpp>
#include <dynamo/interactions/squarewell.hpp>
#include <dynamo/outputplugins/misc.hpp>
#include <dynamo/outputplugins/msd.hpp>
#include <dynamo/systems/andersenThermostat.hpp>
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
  const double diameter = 1.6;
  const double lambda = 1.0;
  const double welldepth = 1.5;
  const double elasticity = 1.0;
  const double bondinner = 0.9;
  const double bondouter = 1.1;
  const size_t N = 50;
  const size_t kT = 1;

  RNG.seed(std::random_device()());
  Sim.ranGenerator.seed(std::random_device()());

  Sim.dynamics = dynamo::shared_ptr<dynamo::Dynamics>(new dynamo::DynNewtonian(&Sim));
  Sim.BCs = dynamo::shared_ptr<dynamo::BoundaryCondition>(new dynamo::BCNone(&Sim));
  Sim.ptrScheduler = dynamo::shared_ptr<dynamo::SNeighbourList>(new dynamo::SNeighbourList(&Sim, new dynamo::CBTFEL<dynamo::HeapPEL>()));
  Sim.primaryCellSize = dynamo::Vector{50, 50, 50};
  Sim.interactions.push_back(dynamo::shared_ptr<dynamo::Interaction>(new dynamo::ISquareBond(&Sim, bondinner, bondouter / bondinner, elasticity, new dynamo::IDPairRangeChains(0, N - 1, N), "Bonds")));
  Sim.interactions.push_back(dynamo::shared_ptr<dynamo::Interaction>(new dynamo::ISquareWell(&Sim, diameter, lambda, welldepth, elasticity, new dynamo::IDPairRangeAll(), "Bulk")));
  Sim.addSpecies(dynamo::shared_ptr<dynamo::Species>(new dynamo::SpPoint(&Sim, new dynamo::IDRangeAll(&Sim), 1.0, "Bulk", 0)));

  for (size_t i = 0; i < N; ++i)
    Sim.particles.push_back(dynamo::Particle(dynamo::Vector{0 + (bondinner + 0.95 * (bondouter - bondinner)) * i, 0, 0}, getRandVelVec() * Sim.units.unitVelocity(), Sim.particles.size()));

  Sim.systems.push_back(dynamo::shared_ptr<dynamo::System>(new dynamo::SysAndersen(&Sim, 0.001 / Sim.N(), kT, "Thermostat")));

  Sim.ensemble = dynamo::Ensemble::loadEnsemble(Sim);
  dynamo::InputPlugin(&Sim, "Rescaler").zeroMomentum();
  dynamo::InputPlugin(&Sim, "Rescaler").rescaleVels(1.0);
  BOOST_CHECK_EQUAL(Sim.N(), N);
}

BOOST_AUTO_TEST_CASE( Equilibrium_Simulation )
{
  {
    dynamo::Simulation Sim;
    init(Sim);
    Sim.writeXMLfile("Polymer.xml");
  }

  dynamo::Simulation Sim;
  Sim.loadXMLfile("Polymer.xml");

  Sim.endEventCount = 1000000;
  Sim.initialise();
  while (Sim.runSimulationStep()) {}

  Sim.reset();
  Sim.endEventCount = 1000000;
  Sim.addOutputPlugin("Misc"); 
  Sim.initialise();
  while (Sim.runSimulationStep()) {}
  
  const double expectedMFT = 0.054058793117007897;
  //Grab the output plugins
  dynamo::OPMisc& opMisc = *Sim.getOutputPlugin<dynamo::OPMisc>();
  //Check the mean free time is roughly what is expected
  double MFT = opMisc.getMFT();
  BOOST_CHECK_CLOSE(MFT, expectedMFT, 1);

  BOOST_CHECK_MESSAGE(Sim.checkSystem() <= 2, "There are more than three invalid states in the final configuration");
}
