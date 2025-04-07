#define BOOST_TEST_MODULE SquareWell_test
#include <boost/test/included/unit_test.hpp>
#include <dynamo/inputplugins/cells/include.hpp>
#include <dynamo/inputplugins/compression.hpp>
#include <dynamo/inputplugins/include.hpp>
#include <dynamo/interactions/squarewell.hpp>
#include <dynamo/outputplugins/misc.hpp>
#include <dynamo/ranges/include.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/species/point.hpp>
#include <dynamo/systems/andersenThermostat.hpp>
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

void init(dynamo::Simulation &Sim, double density = 0.5) {
  RNG.seed(std::random_device()());
  Sim.ranGenerator.seed(std::random_device()());

  const double elasticity = 1.0;
  const double lambda = 1.5;
  const double welldepth = 1.0;

  std::unique_ptr<dynamo::UCell> packptr(
      new dynamo::CUFCC(std::array<long, 3>{{7, 7, 7}}, dynamo::Vector{1, 1, 1},
                        new dynamo::UParticle()));
  packptr->initialise();
  std::vector<dynamo::Vector> latticeSites(
      packptr->placeObjects(dynamo::Vector{0, 0, 0}));
  Sim.primaryCellSize = dynamo::Vector{1, 1, 1};

  double simVol = 1.0;
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    simVol *= Sim.primaryCellSize[iDim];

  double particleDiam = std::cbrt(simVol * density / latticeSites.size());

  Sim.interactions.push_back(dynamo::shared_ptr<dynamo::Interaction>(
      new dynamo::ISquareWell(&Sim, particleDiam, lambda, welldepth, elasticity,
                              new dynamo::IDPairRangeAll(), "Bulk")));
  Sim.addSpecies(dynamo::shared_ptr<dynamo::Species>(
      new dynamo::SpPoint(&Sim, new dynamo::IDRangeAll(&Sim), 1.0, "Bulk", 0)));
  Sim.units.setUnitLength(particleDiam);
  Sim.units.setUnitTime(particleDiam);

  unsigned long nParticles = 0;
  Sim.particles.reserve(latticeSites.size());
  for (const dynamo::Vector &position : latticeSites)
    Sim.particles.push_back(dynamo::Particle(
        position, getRandVelVec() * Sim.units.unitVelocity(), nParticles++));

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

BOOST_AUTO_TEST_CASE(NVE_Simulation) {
  {
    dynamo::Simulation Sim;
    init(Sim);
    Sim.writeXMLfile("SWNVE.xml");
  }

  dynamo::Simulation Sim;
  Sim.loadXMLfile("SWNVE.xml");

  Sim.endEventCount = 100000;
  Sim.addOutputPlugin("Misc");
  Sim.initialise();
  const double totalEinit =
      Sim.getOutputPlugin<dynamo::OPMisc>()->getTotalEnergy();
  while (Sim.runSimulationStep()) {
  }
  const double totalEequil =
      Sim.getOutputPlugin<dynamo::OPMisc>()->getTotalEnergy();
  BOOST_CHECK_CLOSE(totalEinit, totalEequil, 0.000000001);

  Sim.endEventCount += 100000;
  Sim.reset();
  Sim.addOutputPlugin("Misc");
  Sim.initialise();
  const double totalEprerun =
      Sim.getOutputPlugin<dynamo::OPMisc>()->getTotalEnergy();
  BOOST_CHECK_CLOSE(totalEequil, totalEprerun, 0.000000001);
  while (Sim.runSimulationStep()) {
  }
  const double totalEfinal =
      Sim.getOutputPlugin<dynamo::OPMisc>()->getTotalEnergy();
  BOOST_CHECK_CLOSE(totalEprerun, totalEfinal, 0.000000001);

  // Check that the momentum is around 0
  dynamo::Vector momentum =
      Sim.getOutputPlugin<dynamo::OPMisc>()->getCurrentMomentum();
  BOOST_CHECK_SMALL(momentum.nrm() / Sim.units.unitMomentum(), 0.0000000001);

  BOOST_CHECK_MESSAGE(
      Sim.checkSystem() <= 2,
      "There are more than two invalid states in the final configuration");
}

BOOST_AUTO_TEST_CASE(NVT_Simulation) {
  dynamo::Simulation Sim;
  init(Sim);

  Sim.systems.push_back(
      dynamo::shared_ptr<dynamo::System>(new dynamo::SysAndersen(
          &Sim, 0.036 / Sim.N(), 1.0 * Sim.units.unitEnergy(), "Thermostat")));
  Sim.ensemble = dynamo::Ensemble::loadEnsemble(Sim);

  Sim.eventPrintInterval = 50000;
  Sim.endEventCount = 300000;
  Sim.addOutputPlugin("Misc");
  Sim.initialise();
  while (Sim.runSimulationStep()) {
  }

  Sim.reset();
  Sim.endEventCount = 100000;
  Sim.addOutputPlugin("Misc");
  Sim.initialise();

  while (Sim.runSimulationStep()) {
  }

  // Check the mean free time is roughly what is expected
  double MFT = Sim.getOutputPlugin<dynamo::OPMisc>()->getMFT();
  BOOST_CHECK_CLOSE(MFT, 0.0368185, 5);

  // Check the temperature is constant at 1
  double Temperature = Sim.getOutputPlugin<dynamo::OPMisc>()->getCurrentkT() /
                       Sim.units.unitEnergy();
  BOOST_CHECK_CLOSE(Temperature, 1.0, 8);

  BOOST_CHECK_MESSAGE(
      Sim.checkSystem() <= 2,
      "There are more than two invalid states in the final configuration");
}

BOOST_AUTO_TEST_CASE(Compression_Simulation) {
  dynamo::Simulation Sim;
  init(Sim, 0.1);
  Sim.ensemble = dynamo::Ensemble::loadEnsemble(Sim);

  const double growthRate = 1;
  const double targetDensity = 0.9;
  Sim.endEventCount = 1000000;
  Sim.addOutputPlugin("Misc");

  dynamo::shared_ptr<dynamo::IPCompression> compressPlug(
      new dynamo::IPCompression(&Sim, growthRate));
  compressPlug->MakeGrowth();

  compressPlug->limitDensity(targetDensity);

  Sim.initialise();
  while (Sim.runSimulationStep()) {
  }

  while (Sim.runSimulationStep()) {
  }
  compressPlug->RestoreSystem();

  BOOST_CHECK_CLOSE(Sim.getNumberDensity() * Sim.units.unitVolume(),
                    targetDensity, 0.000000001);
  BOOST_CHECK_CLOSE(Sim.getPackingFraction(),
                    Sim.getNumberDensity() * Sim.units.unitVolume() * M_PI /
                        6.0,
                    0.000000001);
  BOOST_CHECK_MESSAGE(Sim.checkSystem() <= 1,
                      "After compression, there are more than one invalid "
                      "states in the final configuration");
}
