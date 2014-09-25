#define BOOST_TEST_MODULE 2DSteppedPotential_test
#include <boost/test/unit_test.hpp>
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
#include <dynamo/interactions/stepped.hpp>
#include <dynamo/interactions/potentials/potential.hpp>
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

  Sim.dynamics = dynamo::shared_ptr<dynamo::Dynamics>(new dynamo::DynNewtonian(&Sim));
  Sim.BCs = dynamo::shared_ptr<dynamo::BoundaryCondition>(new dynamo::BCPeriodic(&Sim));
  Sim.ptrScheduler = dynamo::shared_ptr<dynamo::SNeighbourList>(new dynamo::SNeighbourList(&Sim, new DefaultSorter()));

  std::unique_ptr<dynamo::UCell> packptr(new dynamo::CUSC(std::array<long, 3>{{128, 128, 1}}, dynamo::Vector(1,1,1), new dynamo::UParticle()));
  packptr->initialise();
  std::vector<dynamo::Vector> latticeSites(packptr->placeObjects(dynamo::Vector(0,0,0)));
  
  //As we're in 2D we need to take the square root
  double particleDiam = std::sqrt(density / latticeSites.size());
  Sim.units.setUnitLength(particleDiam);
  Sim.units.setUnitTime(particleDiam); 

  Sim.primaryCellSize = dynamo::Vector(1, 1, 4 * particleDiam);
  
  typedef std::pair<double,double> Step;
  std::vector<Step> steps;
  steps.push_back(Step{1.0,0.1});
  steps.push_back(Step{0.9,0.2});
  steps.push_back(Step{0.8,0.3});
  steps.push_back(Step{0.7,0.4});
  steps.push_back(Step{0.6,0.5});
  steps.push_back(Step{0.5,0.6});
  steps.push_back(Step{0.4,0.7});
  steps.push_back(Step{0.3,0.8});
  steps.push_back(Step{0.2,0.9});
  steps.push_back(Step{0.1,1.0});
  Sim.interactions.push_back(dynamo::shared_ptr<dynamo::Interaction>(new dynamo::IStepped(&Sim, dynamo::shared_ptr<dynamo::Potential>(new dynamo::PotentialStepped(steps, false)), new dynamo::IDPairRangeAll(), "Bulk", particleDiam, 1.0)));
    
  Sim.addSpecies(dynamo::shared_ptr<dynamo::Species>(new dynamo::SpPoint(&Sim, new dynamo::IDRangeAll(&Sim), 1.0, "Bulk", 0)));

  unsigned long nParticles = 0;
  Sim.particles.reserve(latticeSites.size());
  for (const dynamo::Vector & position : latticeSites)
    Sim.particles.push_back(dynamo::Particle(position, getRandVelVec() * Sim.units.unitVelocity(), nParticles++));

  for (auto& particle : Sim.particles)
    particle.getVelocity()[2] = 0;

  dynamo::InputPlugin(&Sim, "Rescaler").zeroMomentum();
  dynamo::InputPlugin(&Sim, "Rescaler").rescaleVels(2.0 / 3.0);
    
  Sim.ensemble = dynamo::Ensemble::loadEnsemble(Sim);

  BOOST_CHECK_EQUAL(Sim.N(), 128 * 128);
}

BOOST_AUTO_TEST_CASE( Equilibrium_Simulation )
{
  {
    dynamo::Simulation Sim;
    init(Sim, 1.0);
    Sim.writeXMLfile("2Dstepped.xml");
  }
  std::cout.flush();
  dynamo::Simulation Sim;
  Sim.loadXMLfile("2Dstepped.xml");

  Sim.endEventCount = 1000000;
  Sim.addOutputPlugin("Misc");
  Sim.initialise();
  while (Sim.runSimulationStep()) {}

  Sim.reset();
  Sim.endEventCount = 1000000;
  Sim.addOutputPlugin("Misc"); 
  Sim.initialise();
  while (Sim.runSimulationStep()) {}
  
  const double expectedMFT = 0.0419518;

  //Grab the output plugins
  dynamo::OPMisc& opMisc = *Sim.getOutputPlugin<dynamo::OPMisc>();

  //Check the mean free time is roughly what is expected
  double MFT = opMisc.getMFT();
  BOOST_CHECK_CLOSE(MFT, expectedMFT, 1);
  
  BOOST_CHECK_MESSAGE(Sim.checkSystem() <= 2, "There are more than two invalid states in the final configuration");
}
