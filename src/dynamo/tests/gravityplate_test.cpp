#define BOOST_TEST_MODULE GravityPlate_test
#include <boost/test/included/unit_test.hpp>
#include <dynamo/BC/include.hpp>
#include <dynamo/dynamics/gravity.hpp>
#include <dynamo/inputplugins/cells/include.hpp>
#include <dynamo/inputplugins/compression.hpp>
#include <dynamo/inputplugins/include.hpp>
#include <dynamo/interactions/hardsphere.hpp>
#include <dynamo/locals/lwall.hpp>
#include <dynamo/outputplugins/misc.hpp>
#include <dynamo/ranges/IDPairRangeAll.hpp>
#include <dynamo/ranges/IDRangeAll.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/species/point.hpp>
#include <random>

std::mt19937 RNG;

dynamo::Vector getRandVelVec() {
  // See http://mathworld.wolfram.com/SpherePointPicking.html
  std::normal_distribution<> normal_dist(0.0, (1.0 / sqrt(double(NDIM))));

  dynamo::Vector tmpVec;
  for (size_t iDim = 0; iDim < NDIM; iDim++)
    tmpVec[iDim] = normal_dist(RNG);

  return tmpVec;
}

void init(dynamo::Simulation &Sim, const double density) {
  RNG.seed(std::random_device()());
  Sim.ranGenerator.seed(std::random_device()());

  const double elasticity = 1.0;

  std::unique_ptr<dynamo::UCell> packptr(
      new dynamo::CUFCC(std::array<long, 3>{{7, 7, 7}}, dynamo::Vector{1, 1, 1},
                        new dynamo::UParticle()));
  packptr->initialise();
  std::vector<dynamo::Vector> latticeSites(
      packptr->placeObjects(dynamo::Vector{0, 0, 0}));

  Sim.primaryCellSize = dynamo::Vector{1, 1, 1};

  double particleDiam = std::cbrt(density / latticeSites.size());

  Sim.dynamics = dynamo::shared_ptr<dynamo::Dynamics>(
      new dynamo::DynGravity(&Sim, dynamo::Vector{0, -particleDiam, 0}));
  Sim.BCs =
      dynamo::shared_ptr<dynamo::BoundaryCondition>(new dynamo::BCNone(&Sim));

  Sim.interactions.push_back(dynamo::shared_ptr<dynamo::Interaction>(
      new dynamo::IHardSphere(&Sim, particleDiam, elasticity,
                              new dynamo::IDPairRangeAll(), "Bulk")));
  Sim.addSpecies(dynamo::shared_ptr<dynamo::Species>(
      new dynamo::SpPoint(&Sim, new dynamo::IDRangeAll(&Sim), 1.0, "Bulk", 0)));
  Sim.units.setUnitLength(particleDiam);

  Sim.locals.push_back(dynamo::shared_ptr<dynamo::Local>(new dynamo::LWall(
      &Sim, 1.0, particleDiam, dynamo::Vector{0, 1, 0},
      dynamo::Vector{0, -0.5 * Sim.primaryCellSize[1] - 0.5 * particleDiam, 0},
      "GroundPlate", new dynamo::IDRangeAll(&Sim))));

  unsigned long nParticles = 0;
  Sim.particles.reserve(latticeSites.size());
  for (const dynamo::Vector &position : latticeSites)
    Sim.particles.push_back(dynamo::Particle(
        0.999 * position, getRandVelVec() * Sim.units.unitVelocity(),
        nParticles++));

  Sim.ensemble = dynamo::Ensemble::loadEnsemble(Sim);

  dynamo::InputPlugin(&Sim, "Rescaler").zeroMomentum();
  dynamo::InputPlugin(&Sim, "Rescaler").rescaleVels(1.0);

  BOOST_CHECK_EQUAL(Sim.N(), 1372);
  BOOST_CHECK_CLOSE(Sim.getNumberDensity() * Sim.units.unitVolume(), density,
                    0.000000001);
  BOOST_CHECK_CLOSE(Sim.getPackingFraction(),
                    Sim.getNumberDensity() * Sim.units.unitVolume() * M_PI /
                        6.0,
                    0.000000001);
}

BOOST_AUTO_TEST_CASE(Simulation) {
  {
    dynamo::Simulation Sim;
    init(Sim, 0.1);
    Sim.writeXMLfile("HSgravityplate.xml");
  }

  dynamo::Simulation Sim;
  Sim.loadXMLfile("HSgravityplate.xml");

  Sim.endEventCount = 100000;
  Sim.addOutputPlugin("Misc");
  Sim.initialise();
  while (Sim.runSimulationStep()) {
  }
  const double expectedMFT = 3.55501052762802;

  // Grab the output plugins
  dynamo::OPMisc &opMisc = *Sim.getOutputPlugin<dynamo::OPMisc>();

  // Check the mean free time is roughly what is expected
  double MFT = opMisc.getMFT();
  BOOST_CHECK_CLOSE(MFT, expectedMFT, 10);

  BOOST_CHECK_MESSAGE(
      Sim.checkSystem() <= 1,
      "There are more than two invalid states in the final configuration");
}
