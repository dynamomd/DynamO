#define BOOST_TEST_MODULE Hardsphere_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <dynamo/simulation.hpp>
#include <dynamo/BC/include.hpp>
#include <dynamo/species/point.hpp>
#include <dynamo/systems/andersenThermostat.hpp>
#include <dynamo/schedulers/include.hpp>
#include <dynamo/inputplugins/include.hpp>
#include <magnet/exception.hpp>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

BOOST_AUTO_TEST_CASE( mean_free_time )
{
  dynamo::Simulation Sim;
  Sim.ranGenerator.seed(42);

  const double density = 0.5;
  const double elasticity = 1.0;

  Sim.dynamics = dynamo::shared_ptr<dynamo::Dynamics>(new dynamo::DynNewtonian(&Sim));
  Sim.BCs = dynamo::shared_ptr<dynamo::BoundaryCondition>(new dynamo::BCPeriodic(&Sim));
  Sim.ptrScheduler = dynamo::shared_ptr<dynamo::SNeighbourList>(new dynamo::SNeighbourList(&Sim, new dynamo::DefaultSorter()));
  std::unique_ptr<UCell> packptr(standardPackingHelper(new dynamo::UParticle()));
  packptr->initialise();
  std::vector<dynamo::Vector> latticeSites(packptr->placeObjects(Vector(0,0,0)));
  Sim.primaryCellSize = dynamo::Vector(1,1,1);

  double simVol = 1.0;
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    simVol *= Sim.primaryCellSize[iDim];

  double particleDiam = std::cbrt(simVol * density / latticeSites.size());
  Sim.interactions.push_back(dynamo::shared_ptr<dynamo::Interaction>(new dynamo::IHardSphere(&Sim, particleDiam, elasticity, new dynamo::IDPairRangeAll(), "Bulk")));
  Sim.addSpecies(dynamo::shared_ptr<dynamo::Species>(new dynamo::SpPoint(&Sim, new dynamo::IDRangeAll(&Sim), 1.0, "Bulk", 0, "Bulk")));
  Sim.units.setUnitLength(particleDiam);

  unsigned long nParticles = 0;
  Sim.particles.reserve(latticeSites.size());
  for (const dynamo::Vector & position : latticeSites)
    Sim.particles.push_back(dynamo::Particle(position, getRandVelVec() * Sim.units.unitVelocity(), nParticles++));
  
  dynamo::InputPlugin(&sim, "Rescaler").zeroMomentum();
  dynamo::InputPlugin(&sim, "Rescaler").rescaleVels(1.0);
}
