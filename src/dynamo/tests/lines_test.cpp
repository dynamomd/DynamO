#define BOOST_TEST_MODULE Lines_test
#include <boost/test/included/unit_test.hpp>
#include <dynamo/BC/LEBC.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/inputplugins/cells/include.hpp>
#include <dynamo/inputplugins/include.hpp>
#include <dynamo/interactions/lines.hpp>
#include <dynamo/outputplugins/misc.hpp>
#include <dynamo/ranges/IDPairRangeAll.hpp>
#include <dynamo/ranges/IDRangeAll.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/species/sphericalTop.hpp>
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
  const size_t N = 1000;

  dynamo::CURandom packroutine(N, dynamo::Vector{1, 1, 1},
                               new dynamo::UParticle());
  packroutine.initialise();
  std::vector<dynamo::Vector> latticeSites(
      packroutine.placeObjects(dynamo::Vector{0, 0, 0}));
  Sim.BCs = dynamo::shared_ptr<dynamo::BoundaryCondition>(
      new dynamo::BCPeriodic(&Sim));
  double particleDiam = std::cbrt(density / N);
  Sim.interactions.push_back(dynamo::shared_ptr<dynamo::Interaction>(
      new dynamo::ILines(&Sim, particleDiam, elasticity,
                         new dynamo::IDPairRangeAll(), "Bulk")));
  Sim.addSpecies(dynamo::shared_ptr<dynamo::Species>(new dynamo::SpSphericalTop(
      &Sim, new dynamo::IDRangeAll(&Sim), 1.0, "Bulk", 0,
      particleDiam * particleDiam / 12.0)));
  Sim.units.setUnitLength(particleDiam);

  unsigned long nParticles = 0;
  Sim.particles.reserve(latticeSites.size());
  for (const dynamo::Vector &position : latticeSites)
    Sim.particles.push_back(dynamo::Particle(
        position, getRandVelVec() * Sim.units.unitVelocity(), nParticles++));

  Sim.dynamics->initOrientations();

  Sim.ensemble = dynamo::Ensemble::loadEnsemble(Sim);

  dynamo::InputPlugin(&Sim, "Rescaler").zeroMomentum();
  dynamo::InputPlugin(&Sim, "Rescaler").rescaleVels(1.0);

  BOOST_CHECK_EQUAL(Sim.N(), N);
  BOOST_CHECK_CLOSE(Sim.getNumberDensity() * Sim.units.unitVolume(), density,
                    0.000000001);
  BOOST_CHECK_CLOSE(Sim.getPackingFraction(), 0, 0.000000001);
}

BOOST_AUTO_TEST_CASE(Equilibrium_Simulation) {
  try {
    const double density = 0.1;
    // Create the sim and test save/loading
    {
      dynamo::Simulation Sim;
      init(Sim, density);
      Sim.writeXMLfile("lines.xml");
    }

    dynamo::Simulation Sim;
    Sim.loadXMLfile("lines.xml");

    Sim.eventPrintInterval = 10000;
    Sim.endEventCount = 100000;
    Sim.addOutputPlugin("Misc");
    Sim.initialise();
    while (Sim.runSimulationStep()) {
    }

    const double expectedMFT = 1.0 / (1.237662399 * density);

    // Grab the output plugins
    dynamo::OPMisc &opMisc = *Sim.getOutputPlugin<dynamo::OPMisc>();
    // Check the mean free time is roughly what is expected
    double MFT = opMisc.getMFT() / Sim.units.unitTime();
    BOOST_CHECK_CLOSE(MFT, expectedMFT, 4);

    // Check that the momentum is still around 0
    dynamo::Vector momentum =
        Sim.getOutputPlugin<dynamo::OPMisc>()->getCurrentMomentum();
    BOOST_CHECK_SMALL(momentum.nrm() / Sim.units.unitMomentum(), 0.0000000001);
  } catch (std::exception &cep) {
    std::cout << cep.what() << std::endl;
    throw;
  }
}
