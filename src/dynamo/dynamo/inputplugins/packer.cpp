/*  dynamo:- Event driven molecular dynamics simulator
    http://www.dynamomd.org
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

    This program is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    version 3 as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <boost/tokenizer.hpp>
#include <cmath>
#include <dynamo/BC/include.hpp>
#include <dynamo/dynamics/include.hpp>
#include <dynamo/ensemble.hpp>
#include <dynamo/globals/include.hpp>
#include <dynamo/inputplugins/cells/include.hpp>
#include <dynamo/inputplugins/packer.hpp>
#include <dynamo/interactions/include.hpp>
#include <dynamo/interactions/potentials/lennard_jones.hpp>
#include <dynamo/locals/lcylinder.hpp>
#include <dynamo/locals/lwall.hpp>
#include <dynamo/locals/oscillatingplate.hpp>
#include <dynamo/particle.hpp>
#include <dynamo/ranges/include.hpp>
#include <dynamo/schedulers/include.hpp>
#include <dynamo/schedulers/sorters/CBTFEL.hpp>
#include <dynamo/schedulers/sorters/heapPEL.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/species/include.hpp>
#include <dynamo/systems/DSMCspheres.hpp>
#include <dynamo/systems/andersenThermostat.hpp>
#include <dynamo/systems/rescale.hpp>
#include <dynamo/systems/rotateGravity.hpp>
#include <dynamo/systems/sleep.hpp>
#include <dynamo/topology/include.hpp>
#include <magnet/exception.hpp>
#include <magnet/math/matrix.hpp>
#include <memory>

namespace dynamo {

namespace {
struct speciesData {
  double diameter;
  double lambda;
  double mass;
  double wellDepth;
  double molfraction;
  size_t idStart;
  size_t idEnd;
};
} // namespace

IPPacker::IPPacker(po::variables_map &vm2, dynamo::Simulation *tmp)
    : SimBase(tmp, "SysPacker"), vm(vm2) {}

bool mySortPredictate(const Vector &v1, const Vector &v2) {
  return v1[0] > v2[0];
}

po::options_description IPPacker::getOptions() {
  po::options_description retval("Packer options");

  retval.add_options()(
      "pack-mode,m", po::value<size_t>(),
      "Chooses the system to pack (construct)"
      "\nPacker Modes:"
      "\n0:  Monocomponent hard spheres"
      "\n1:  Mono/Multi-component square wells"
      "\n2:  Random walk of an isolated attractive polymer"
      "\n3:  Load a config and pack it, you will need to reset the "
      "interactions etc."
      "\n4:  Monocomponent (in)elastic hard spheres in LEBC (shearing)"
      "\n5:  Walk an isolated spiral/helix"
      "\n6:  Monocomponent hard spheres confined by two walls, aspect ratio is "
      "set by the number of cells"
      "\n7:  Ring/Linear polymer, dropped as a straight rod"
      "\n8:  Binary Hard Spheres"
      "\n9:  Hard needle system"
      "\n10: Monocomponent hard spheres using DSMC interactions"
      "\n11: (DEPRECATED) Monocomponent hard spheres sheared using DSMC "
      "interactions"
      "\n12: Binary hard spheres using DSMC interactions"
      "\n13: Crystal pack of sheared lines"
      "\n14: Packing of spheres and linear rods made from stiff polymers"
      "\n15: Monocomponent hard-parallel cubes"
      "\n16: Stepped Potential"
      "\n17: (DEPRECATED) Monocomponent hard spheres using Ring DSMC "
      "interactions"
      "\n18: (DEPRECATED) Monocomponent sheared hard spheres using Ring DSMC "
      "interactions"
      "\n19: Oscillating plates bounding a system"
      "\n20: Load a set of triangles and plate it with spheres"
      "\n21: Pack a cylinder with spheres"
      "\n22: Infinite system with spheres falling onto a plate with gravity"
      "\n23: Funnel test for static spheres in gravity"
      "\n24: Random walk of an isolated MJ model polymer"
      "\n25: Funnel and cup simulation (with sleepy particles)"
      "\n26: Polydisperse (Gaussian) hard spheres in LEBC (shearing)"
      "\n27: Crystal pack of snowmen molecules"
      "\n28: Rotating drum made out of particles.");

  return retval;
}

void IPPacker::initialise() {
  std::string defaultOptionText =
      " Options\n"
      "  -C [ --NCells ] arg (=7)    Set the default number of lattice "
      "unit-cells in each direction.\n"
      "  -x [ --xcell ] arg          Number of unit-cells in the x dimension.\n"
      "  -y [ --ycell ] arg          Number of unit-cells in the y dimension.\n"
      "  -z [ --zcell ] arg          Number of unit-cells in the z dimension.\n"
      "  --rectangular-box           Set the simulation box to be rectangular "
      "so that the x,y,z cells also specify the simulation aspect ratio.\n"
      "  -d [ --density ] arg (=0.5) System density.\n"
      "  --i1 arg (=FCC)             Lattice type (0=FCC, 1=BCC, 2=SC, "
      "3=HCP)\n";

  switch (vm["pack-mode"].as<size_t>()) {
  case 0: {
    if (vm.count("help")) {
      std::cout << "\nMode 0: Monocomponent hard spheres\n"
                << defaultOptionText
                << "  --i2 arg (disabled)    Adds a temperature rescale event "
                   "every x events\n"
                   "  --f1 arg (=1)             Sets the elasticity of the "
                   "hard spheres\n"
                   "  --f2 arg (=1)             Sets the tangential elasticity "
                   "of the hard spheres (=1 disables rotation)\n";
      exit(1);
    }
    // Pack of hard spheres
    // Pack the system, determine the number of particles
    std::unique_ptr<UCell> packptr(standardPackingHelper(new UParticle()));
    packptr->initialise();

    std::vector<Vector> latticeSites(packptr->placeObjects(Vector{0, 0, 0}));

    Sim->primaryCellSize = packptr->systemDims();

    double simVol = 1.0;

    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      simVol *= Sim->primaryCellSize[iDim];

    double particleDiam =
        pow(simVol * vm["density"].as<double>() / latticeSites.size(),
            double(1.0 / 3.0));

    bool twoD = false;

    if (vm.count("rectangular-box") &&
        (vm.count("i1") && vm["i1"].as<size_t>() == 2)) {
      std::array<long, 3> cells = getCells();
      if ((cells[0] == 1) || (cells[1] == 1) || (cells[2] == 1)) {
        twoD = true;
        derr << "Warning! Now assuming that you're trying to set up a 2D "
                "simulation!\n"
                "I'm going to temporarily calculate the density by the 2D "
                "definition!"
             << std::endl;

        size_t dimension;
        if (cells[0] == 1)
          dimension = 0;
        if (cells[1] == 1)
          dimension = 1;
        if (cells[2] == 1)
          dimension = 2;

        particleDiam =
            std::sqrt(simVol * vm["density"].as<double>() /
                      (Sim->primaryCellSize[dimension] * latticeSites.size()));

        dout << "I'm changing what looks like the unused box dimension ("
             << dimension
             << ") to the smallest value allowed by the neighbourlist "
                "implementation (slightly more than 4 particle diameters)"
             << std::endl;

        Sim->primaryCellSize[dimension] = 4.0000001 * particleDiam;
      }
    }

    double elasticity = 1.0;
    if (vm.count("f1"))
      elasticity = vm["f1"].as<double>();

    if (vm.count("f2") && (vm["f2"].as<double>() != 1.0)) {
      // Simulation with rotation
      Sim->interactions.push_back(shared_ptr<Interaction>(
          new IHardSphere(Sim, particleDiam, elasticity, vm["f2"].as<double>(),
                          new IDPairRangeAll(), "Bulk")));
      Sim->addSpecies(shared_ptr<Species>(
          new SpSphericalTop(Sim, new IDRangeAll(Sim), 1.0, "Bulk", 0,
                             2.0 * particleDiam * particleDiam / (5.0 * 4.0))));
      Sim->units.setUnitLength(particleDiam);

      size_t nParticles = 0;
      Sim->particles.reserve(latticeSites.size());
      for (const Vector &position : latticeSites) {
        Sim->particles.push_back(
            Particle(position, getRandVelVec() * Sim->units.unitVelocity(),
                     nParticles++));
        if (twoD)
          Sim->particles.back().getVelocity()[2] = 0;
      }

      Sim->dynamics->initOrientations();
    } else {
      // Simulation without rotation
      Sim->interactions.push_back(shared_ptr<Interaction>(new IHardSphere(
          Sim, particleDiam, elasticity, new IDPairRangeAll(), "Bulk")));
      Sim->addSpecies(shared_ptr<Species>(
          new SpPoint(Sim, new IDRangeAll(Sim), 1.0, "Bulk", 0)));
      Sim->units.setUnitLength(particleDiam);

      size_t nParticles = 0;
      Sim->particles.reserve(latticeSites.size());
      for (const Vector &position : latticeSites) {
        Sim->particles.push_back(
            Particle(position, getRandVelVec() * Sim->units.unitVelocity(),
                     nParticles++));
        if (twoD)
          Sim->particles.back().getVelocity()[2] = 0;
      }
    }

    const double kT = 1.0 * Sim->units.unitEnergy();
    if (vm.count("i2"))
      Sim->systems.push_back(shared_ptr<System>(
          new SysRescale(Sim, vm["i2"].as<size_t>(), "RescalerEvent", kT)));

    break;
  }
  case 1: {
    if (vm.count("help")) {
      std::cout
          << "\nMode 1: Mono/Multi-component square wells\n"
          << defaultOptionText
          << "  --f1 arg (=1.5)             Well width factor (also known as "
             "lambda)\n"
             "  --f2 arg (=1)               Well Depth (negative values create "
             "square shoulders)\n"
             "  --s1 arg (monocomponent)    Instead of f1 and f2, you can "
             "specify a multicomponent system using this option. You need to "
             "pass the the parameters for each species as follows --s1 "
             "\"diameter(d),lambda(l),mass(m),welldepth(e),molefrac(x):d,l,m,e,"
             "x[:...]\"\n";
      exit(1);
    }
    // Pack of square well molecules
    // Pack the system, determine the number of particles
    std::unique_ptr<UCell> packptr(
        new CURandomise(standardPackingHelper(new UParticle())));
    packptr->initialise();

    std::vector<Vector> latticeSites(packptr->placeObjects(Vector{0, 0, 0}));

    Sim->primaryCellSize = packptr->systemDims();

    double simVol = 1.0;

    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      simVol *= Sim->primaryCellSize[iDim];

    double particleDiam =
        pow(simVol * vm["density"].as<double>() / latticeSites.size(),
            double(1.0 / 3.0));

    // Set the unit energy to 1 (assuming the unit of mass is 1);
    Sim->units.setUnitLength(particleDiam);
    Sim->units.setUnitTime(particleDiam);

    if (!vm.count("s1")) {
      // Only one species
      double lambda = 1.5, wellDepth = 1.0;
      if (vm.count("f1"))
        lambda = vm["f1"].as<double>();

      if (vm.count("f2"))
        wellDepth = vm["f2"].as<double>();

      Sim->interactions.push_back(shared_ptr<Interaction>(
          new ISquareWell(Sim, particleDiam, lambda, wellDepth, 1.0,
                          new IDPairRangeAll(), "Bulk")));

      Sim->addSpecies(shared_ptr<Species>(
          new SpPoint(Sim, new IDRangeAll(Sim), 1.0, "Bulk", 0)));
    } else {
      // Multiple species specified by a string

      std::vector<speciesData> speciesList;

      typedef boost::tokenizer<boost::char_separator<char>> tokenizer;

      // Split up over species delimited by ":"
      tokenizer speciesStrings(vm["s1"].as<std::string>(),
                               boost::char_separator<char>(":"));

      double totMoleFrac = 0;
      for (const std::string &species : speciesStrings) {
        tokenizer speciesStringData(species, boost::char_separator<char>(","));
        tokenizer::iterator value_iter = speciesStringData.begin();

        speciesData dat;

        try {
          if (value_iter == speciesStringData.end())
            throw std::runtime_error("Stray : in species definition");
          dat.diameter = boost::lexical_cast<double>(*value_iter);

          if (++value_iter == speciesStringData.end())
            throw std::runtime_error("No lambda specified for a species");
          dat.lambda = boost::lexical_cast<double>(*(value_iter));

          if (++value_iter == speciesStringData.end())
            throw std::runtime_error("No mass specified for a species");
          dat.mass = boost::lexical_cast<double>(*(value_iter));

          if (++value_iter == speciesStringData.end())
            throw std::runtime_error("No well depth specified for a species");
          dat.wellDepth = boost::lexical_cast<double>(*(value_iter));

          if (++value_iter == speciesStringData.end())
            throw std::runtime_error(
                "No mole fraction specified for a species");
          dat.molfraction = boost::lexical_cast<double>(*(value_iter));
          totMoleFrac += dat.molfraction;

          if (++value_iter != speciesStringData.end())
            throw std::runtime_error("Too many comma's");

          speciesList.push_back(dat);
        } catch (std::exception &except) {
          M_throw() << "Malformed square well species data, \"" << species
                    << "\"\n"
                    << except.what();
        }
      }

      // Normalize the mole fraction and calculate the range
      const size_t N = latticeSites.size();
      size_t idStart = 0;
      for (speciesData &dat : speciesList) {
        dat.molfraction /= totMoleFrac;
        dat.idStart = idStart;
        // The minus 0.5 is to make the double to size_t into a
        //"round to nearest minus 1"
        dat.idEnd = (N * dat.molfraction - 0.5) + idStart;
        idStart = dat.idEnd + 1;
      }

      // Chuck the rounding error amount of spheres into the last species
      speciesList.back().idEnd = N - 1;

      // Now we have the particle ranges we should build the species and
      // interactions
      for (size_t spID1(0); spID1 < speciesList.size(); ++spID1) {
        const speciesData &spdat1 = speciesList[spID1];
        std::string sp1Name = "A";
        sp1Name[0] += spID1;
        for (size_t spID2(spID1); spID2 < speciesList.size(); ++spID2) {
          const speciesData &spdat2 = speciesList[spID2];
          std::string sp2Name = "A";
          sp2Name[0] += spID2;
          Sim->interactions.push_back(shared_ptr<Interaction>(new ISquareWell(
              Sim, 0.5 * particleDiam * (spdat1.diameter + spdat2.diameter),
              0.5 * (spdat1.lambda + spdat2.lambda),
              0.5 * (spdat1.wellDepth + spdat2.wellDepth), 1.0,
              new IDPairRangePair(
                  new IDRangeRange(spdat1.idStart, spdat1.idEnd),
                  new IDRangeRange(spdat2.idStart, spdat2.idEnd)),
              sp1Name + sp2Name)));
        }
      }

      for (size_t spID1(0); spID1 < speciesList.size(); ++spID1) {
        const speciesData &spdat1 = speciesList[spID1];
        std::string sp1Name = "A";
        sp1Name[0] += spID1;
        Sim->addSpecies(shared_ptr<Species>(
            new SpPoint(Sim, new IDRangeRange(spdat1.idStart, spdat1.idEnd),
                        spdat1.mass, sp1Name, spID1)));
      }
    }

    unsigned long nParticles = 0;

    Sim->particles.reserve(latticeSites.size());
    for (const Vector &position : latticeSites)
      Sim->particles.push_back(Particle(
          position, getRandVelVec() * Sim->units.unitVelocity(), nParticles++));
    break;
  }
  case 2: {
    if (vm.count("help")) {
      std::cout
          << "\nMode 2: Create an isolated, homo or HP polymer using a random "
             "self-avoiding walk\n"
             "  --i1 arg (=20)              Chain length (No. of monomers)\n"
             "  --f1 arg (=1.6)             Monomer diameter\n"
             "  --f2 arg (=1.5)             Monomer well width factor (also "
             "called lambda)\n"
             "  --f3 arg (=0.9)             Bond inner core\n"
             "  --f4 arg (=1.1)             Bond outer well\n"
             "  --s1 arg (homopolymer)      HP sequence to use (eg 0001010), "
             "defaults to homopolymer if unset\n";
      exit(1);
    }
    // Random walk an isolated attractive homopolymer
    size_t chainlength = 20;

    double sigmin(0.9), sigmax(1.1), sigma(1.6), lambda(1.5);

    if (vm.count("f1"))
      sigma = vm["f1"].as<double>();

    if (vm.count("f2"))
      lambda = vm["f2"].as<double>();

    if (vm.count("f3"))
      sigmin = vm["f3"].as<double>();

    if (vm.count("f4"))
      sigmax = vm["f4"].as<double>();

    if (vm.count("i1"))
      chainlength = vm["i1"].as<size_t>();

    std::string stringseq;
    if (vm.count("s1")) {
      stringseq = vm["s1"].as<std::string>();
      if (!vm.count("i1"))
        chainlength = stringseq.size();
      else if (chainlength != stringseq.size())
        M_throw() << "Error, mismatch between chain length and sequence "
                     "length. You can remove --i1 and let the chain length be "
                     "determined from the sequence length if needed?";
    }

    // Sit the particles 95% away of max distance from each other
    // to help with seriously overlapping wells
    double diamScale = 1.0 / (4.0 * std::max(1.0, std::sqrt(chainlength)) *
                              std::max(lambda * sigma, sigmax));

    CURandWalk sysPack(chainlength,
                       (sigmin + 0.95 * (sigmax - sigmin)) * diamScale,
                       sigma * diamScale, new UParticle());

    sysPack.initialise();

    // Drop them in the middle of the sim
    std::vector<Vector> latticeSites(sysPack.placeObjects(Vector{0, 0, 0}));

    Sim->interactions.push_back(shared_ptr<Interaction>(new ISquareBond(
        Sim, sigmin * diamScale, sigmax / sigmin, 1.0,
        new IDPairRangeChains(0, latticeSites.size() - 1, latticeSites.size()),
        "Bonds")));

    if (vm.count("s1")) {
      // A sequence has been supplied
      std::vector<size_t> seq;
      seq.resize(chainlength, 0);

      // Transcribe the sequence
      bool has0(false), has1(false);
      for (size_t i = 0; i < chainlength; ++i) {
        seq[i] = boost::lexical_cast<size_t>(stringseq[i % stringseq.size()]);
        if (seq[i])
          has1 = true;
        else
          has0 = true;
        if (seq[i] > 1)
          M_throw()
              << "Dynamod only supports 2 types of monomers, make a sample "
                 "chain and edit the configuration file by hand to use more";
      }

      if (has1 && has0) {
        Sim->interactions.push_back(shared_ptr<Interaction>(
            new ISWSequence(Sim, sigma * diamScale, lambda, 1.0, seq,
                            new IDPairRangeAll(), "Bulk")));

        ISWSequence &interaction(
            dynamic_cast<ISWSequence &>(*Sim->interactions["Bulk"]));
        interaction.getAlphabet().at(0).at(0) = 1.0;

        interaction.getAlphabet().at(1).at(0) = 0.5;

        interaction.getAlphabet().at(0).at(1) = 0.5;
      } else if (has0 && !has1)
        Sim->interactions.push_back(shared_ptr<Interaction>(
            new ISquareWell(Sim, sigma * diamScale, lambda, 1.0, 1.0,
                            new IDPairRangeAll(), "Bulk")));
      else if (has1 && !has0)
        Sim->interactions.push_back(shared_ptr<Interaction>(new IHardSphere(
            Sim, sigma * diamScale, new IDPairRangeAll(), "Bulk")));
    } else
      Sim->interactions.push_back(shared_ptr<Interaction>(
          new ISquareWell(Sim, sigma * diamScale, lambda, 1.0, 1.0,
                          new IDPairRangeAll(), "Bulk")));

    Sim->addSpecies(shared_ptr<Species>(
        new SpPoint(Sim, new IDRangeAll(Sim), 1.0, "Bulk", 0)));

    Sim->units.setUnitLength(diamScale);
    // Set the unit energy to 1 (assuming the unit of mass is 1);
    Sim->units.setUnitTime(diamScale);

    Sim->topology.push_back(
        shared_ptr<Topology>(new TChain(Sim, 1, "HelixPolymer")));

    Sim->topology.back()->addMolecule(new IDRangeAll(Sim));

    Sim->BCs = shared_ptr<BoundaryCondition>(new BCNone(Sim));

    unsigned long nParticles = 0;

    Sim->particles.reserve(latticeSites.size());
    for (const Vector &position : latticeSites)
      Sim->particles.push_back(Particle(
          position, getRandVelVec() * Sim->units.unitVelocity(), nParticles++));
    break;
  }
  case 3: {
    if (vm.count("help")) {
      std::cout << "\nMode 3: Takes an existing configuration and packs images "
                   "of it on a lattice (useful for packing polymers and other "
                   "complex molecules)\n"
                   "        BE WARNED: It is very easy to produce overlaps "
                   "with this mode if the density is too high!"
                << defaultOptionText
                << "  --s1 arg                    Filename of the "
                   "configuration to use as the image\n"
                   "  --f1 arg (=0)               Fraction of images that are "
                   "mirrored before placement (from 0.0 to 1.0)\n";
      exit(1);
    }
    // This packs a system using a file for the unit cell, the
    // density should just be adjusted by hand
    std::string fileName;

    if (!vm.count("s1"))
      M_throw()
          << "You must specify the config file to pack using the --s1 option!";

    fileName = vm["s1"].as<std::string>();

    // Just load the existing configuration
    Sim->loadXMLfile(fileName);

    // Now read the positions that were entered in
    size_t NUnit = Sim->particles.size();

    // Figure out how many units there are
    UCell *tmpPtr = standardPackingHelper(new UParticle());
    tmpPtr->initialise();
    size_t NUnitSites = tmpPtr->placeObjects(Vector{0, 0, 0}).size();
    delete tmpPtr;

    double diamScale = pow(vm["density"].as<double>() / (NUnitSites * NUnit),
                           double(1.0 / 3.0));

    // Now set the size of the system
    Sim->primaryCellSize = Vector{1, 1, 1} / diamScale;
    if (vm.count("rectangular-box"))
      Sim->primaryCellSize = getNormalisedCellDimensions() / diamScale;

    {
      std::vector<Vector> positions;
      for (const Particle &p : Sim->particles)
        positions.push_back(p.getPosition());
      tmpPtr = new UList(positions, diamScale, new UParticle());
    }

    // Delete any loaded capture maps
    for (const shared_ptr<Interaction> &ptr : Sim->interactions)
      if (std::dynamic_pointer_cast<ICapture>(ptr))
        std::dynamic_pointer_cast<ICapture>(ptr)->forgetMap();

    // Use the mirror unit cell if needed
    if (vm.count("f1") && (vm["f1"].as<double>() != 0)) {
      if ((vm["f1"].as<double>() < 0) || (vm["f1"].as<double>() > 1))
        M_throw() << "You must specify a chiral fraction between 0.0 and 1.0";

      tmpPtr = new CUMirror(vm["f1"].as<double>(), tmpPtr);
    }

    Sim->BCs = shared_ptr<BoundaryCondition>(new BCPeriodic(Sim));
    std::unique_ptr<UCell> packptr(standardPackingHelper(tmpPtr));

    packptr->initialise();

    std::vector<Vector> latticeSites(packptr->placeObjects(Vector{0, 0, 0}));

    unsigned long nParticles = 0;
    Sim->particles.clear();
    Sim->particles.reserve(latticeSites.size());
    for (const Vector &position : latticeSites)
      Sim->particles.push_back(
          Particle(position / diamScale,
                   getRandVelVec() * Sim->units.unitVelocity(), nParticles++));
    break;
  }
  case 4: {
    if (vm.count("help")) {
      std::cout << "\nMode 4: Monocomponent (in)elastic hard spheres in LEBC "
                   "(shearing)\n"
                << defaultOptionText
                << "  --f1 arg (=1.0)             Sets the elasticity of the "
                   "hard spheres\n";
      exit(1);
    }
    // FCC simple cubic pack of hard spheres with inelasticity and shearing
    // Pack the system, determine the number of particles
    std::unique_ptr<UCell> packptr(standardPackingHelper(new UParticle()));
    packptr->initialise();

    std::vector<Vector> latticeSites(packptr->placeObjects(Vector{0, 0, 0}));

    Sim->primaryCellSize = packptr->systemDims();

    double simVol = 1.0;

    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      simVol *= Sim->primaryCellSize[iDim];

    double particleDiam =
        pow(simVol * vm["density"].as<double>() / latticeSites.size(),
            double(1.0 / 3.0));

    double alpha = 1.0;

    if (vm.count("f1"))
      alpha = vm["f1"].as<double>();

    Sim->BCs = shared_ptr<BoundaryCondition>(new BCLeesEdwards(Sim));
    const double shearRate = 1;

    Sim->interactions.push_back(shared_ptr<Interaction>(new IHardSphere(
        Sim, particleDiam, alpha, new IDPairRangeAll(), "Bulk")));
    Sim->addSpecies(shared_ptr<Species>(
        new SpPoint(Sim, new IDRangeAll(Sim), 1.0, "Bulk", 0)));
    Sim->units.setUnitLength(particleDiam);

    unsigned long nParticles = 0;
    Sim->particles.reserve(latticeSites.size());
    for (const Vector &position : latticeSites)
      Sim->particles.push_back(Particle(
          position, getRandVelVec() * Sim->units.unitVelocity(), nParticles++));

    // Insert a linear profile, zero momentum then add a vel gradient
    Sim->setCOMVelocity();
    for (Particle &part : Sim->particles)
      part.getVelocity()[0] += part.getPosition()[1] * shearRate;
    break;
  }
  case 5: {
    if (vm.count("help")) {
      std::cout
          << "\nMode 5: Create an isolated, homopolymer using a spiraling "
             "walk\n"
             "  --i1 arg (=20)              Chain length (No. of monomers)\n"
             "  --i2 arg (=9)               Ring length (monomers in one turn "
             "of the spiral)\n"
             "  --f1 arg (=1.6)             Monomer diameter\n"
             "  --f2 arg (=1.5)             Monomer well width factor (also "
             "called lambda)\n"
             "  --f3 arg (=0.9)             Bond inner core\n"
             "  --f4 arg (=1.1)             Bond outer well\n"
             "  --f5 arg (=0.05)            Relative tightness of the helix (0 "
             "is as close as possible, 1 is as far apart as possible)\n";
      exit(1);
    }
    // Helix/spiral layout of molecules
    size_t chainlength = 20;

    if (vm.count("i1"))
      chainlength = vm["i1"].as<size_t>();

    size_t ringlength = 9;

    if (vm.count("i2"))
      ringlength = vm["i2"].as<size_t>();

    double sigmin(0.9), sigmax(1.1), sigma(1.6), lambda(1.5), tightness(0.05);

    if (vm.count("f1"))
      sigma = vm["f1"].as<double>();

    if (vm.count("f2"))
      lambda = vm["f2"].as<double>();

    if (vm.count("f3"))
      sigmin = vm["f3"].as<double>();

    if (vm.count("f4"))
      sigmax = vm["f4"].as<double>();

    if (vm.count("f5"))
      tightness = vm["f5"].as<double>();

    // Sit the particles 95% away of max distance from each other
    // to help with seriously overlapping wells
    double diamScale = 1.0 / chainlength;

    // Space the hard spheres 2% further apart than minimum and set
    // the bonds to 2% max length to coil this as much as possible
    CUHelix sysPack(chainlength, ringlength,
                    (sigmin + tightness * (sigmax - sigmin)) * diamScale,
                    (1.0 + tightness) * sigma * diamScale, new UParticle());

    sysPack.initialise();

    // Drop them in the middle of the sim
    std::vector<Vector> latticeSites(sysPack.placeObjects(Vector{0, 0, 0}));

    Sim->interactions.push_back(shared_ptr<Interaction>(new ISquareBond(
        Sim, sigmin * diamScale, sigmax / sigmin, 1.0,
        new IDPairRangeChains(0, latticeSites.size() - 1, latticeSites.size()),
        "Bonds")));

    Sim->interactions.push_back(shared_ptr<Interaction>(
        new ISquareWell(Sim, sigma * diamScale, lambda, 1.0, 1.0,
                        new IDPairRangeAll(), "Bulk")));

    Sim->addSpecies(shared_ptr<Species>(
        new SpPoint(Sim, new IDRangeAll(Sim), 1.0, "Bulk", 0)));

    Sim->units.setUnitLength(diamScale);
    // Set the unit energy to 1 (assuming the unit of mass is 1);
    Sim->units.setUnitTime(diamScale);

    Sim->topology.push_back(
        shared_ptr<Topology>(new TChain(Sim, 1, "HelixPolymer")));

    Sim->topology.back()->addMolecule(new IDRangeAll(Sim));

    unsigned long nParticles = 0;

    Sim->particles.reserve(latticeSites.size());
    for (const Vector &position : latticeSites)
      Sim->particles.push_back(Particle(
          position, getRandVelVec() * Sim->units.unitVelocity(), nParticles++));
    break;
  }
  case 6: {
    if (vm.count("help")) {
      std::cout << "Mode specific options:\n"
                   "  6:  Monocomponent hard spheres confined by two walls, "
                   "aspect ratio is set by the number of cells\n"
                   "       --f1 : Elasticity of the particle and wall "
                   "collisions [1]\n";
      exit(1);
    }
    std::unique_ptr<UCell> packptr(
        standardPackingHelper(new UParticle(), true));
    packptr->initialise();

    std::vector<Vector> latticeSites(packptr->placeObjects(Vector{0, 0, 0}));

    Sim->primaryCellSize = packptr->systemDims();

    // Cut off the x periodic boundaries
    Sim->BCs = shared_ptr<BoundaryCondition>(new BCPeriodicExceptX(Sim));

    double simVol = 1.0;

    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      simVol *= Sim->primaryCellSize[iDim];

    double particleDiam =
        pow(simVol * vm["density"].as<double>() / latticeSites.size(),
            double(1.0 / 3.0));

    Sim->units.setUnitLength(particleDiam);

    double elasticity = 1;
    if (vm.count("f1"))
      elasticity = vm["f1"].as<double>();

    Sim->locals.push_back(shared_ptr<Local>(new LWall(
        Sim, elasticity, particleDiam, Vector{1, 0, 0},
        Vector{-Sim->primaryCellSize[0] / 2 - 0.5 * particleDiam, 0, 0},
        "LowWall", new IDRangeAll(Sim))));
    Sim->locals.push_back(shared_ptr<Local>(new LWall(
        Sim, elasticity, particleDiam, Vector{-1, 0, 0},
        Vector{Sim->primaryCellSize[0] / 2 + 0.5 * particleDiam, 0, 0},
        "HighWall", new IDRangeAll(Sim))));

    Sim->interactions.push_back(shared_ptr<Interaction>(new IHardSphere(
        Sim, particleDiam, elasticity, new IDPairRangeAll(), "Bulk")));
    Sim->addSpecies(shared_ptr<Species>(
        new SpPoint(Sim, new IDRangeAll(Sim), 1.0, "Bulk", 0)));

    unsigned long nParticles = 0;
    Sim->particles.reserve(latticeSites.size());
    for (const Vector &position : latticeSites)
      Sim->particles.push_back(Particle(
          position, getRandVelVec() * Sim->units.unitVelocity(), nParticles++));
    break;
  }
  case 7: {
    if (vm.count("help")) {
      std::cout
          << "Mode specific options:\n"
             "  7:  Ring/Linear polymer, dropped as a straight rod\n"
             "       --i1 : Chain length (number supplied is multiplied by 2, "
             "e.g. default of 10 gives a 20mer) [10]\n"
             "       --f1 : Bond inner core (>0) [1.0]\n"
             "       --f2 : Bond outer well (>0) [1.05]\n"
             "       --f3 : Well width factor, values <= 1 use a hard sphere "
             "[1.5]\n"
             "       --b1 : If set it drops a linear chain instead of a ring\n";
      exit(1);
    }
    // Just drops a ring polymer, you should crystalize it then
    // pack it for bulk

    size_t chainlength = 10;

    if (vm.count("i1"))
      chainlength = vm["i1"].as<size_t>();

    double sigma(1.0), sigmin(1.0), sigmax(1.05), lambda(1.5);

    if (vm.count("f1"))
      sigmin = vm["f1"].as<double>();

    if (vm.count("f2"))
      sigmax = vm["f2"].as<double>();

    if (vm.count("f3"))
      lambda = vm["f3"].as<double>();

    // 10 % more than double whats needed
    double diamScale = 0.5 / (sigmax * chainlength + 2 * sigma);

    // CUringRod sysPack(chainlength, ((sigmax - sigmin) * 0.95 + sigmin)
    //                                * diamScale, new UParticle());

    CUringSnake sysPack(chainlength,
                        ((sigmax - sigmin) * 0.95 + sigmin) * diamScale,
                        new UParticle());

    sysPack.initialise();

    // Drop them in the middle of the sim
    std::vector<Vector> latticeSites(sysPack.placeObjects(Vector{0, 0, 0}));

    Sim->interactions.push_back(shared_ptr<Interaction>(new ISquareBond(
        Sim, sigmin * diamScale, sigmax / sigmin, 1.0,
        (vm.count("b1"))
            ? static_cast<IDPairRange *>(new IDPairRangeChains(
                  0, latticeSites.size() - 1, latticeSites.size()))
            : static_cast<IDPairRange *>(new IDPairRangeRings(
                  0, latticeSites.size() - 1, latticeSites.size())),
        "Bonds")));

    if (lambda >= 1.0) {
      Sim->units.setUnitLength(diamScale);
      // Set the unit energy to 1 (assuming the unit of mass is 1);
      Sim->units.setUnitTime(diamScale);
      Sim->interactions.push_back(shared_ptr<Interaction>(
          new ISquareWell(Sim, sigma * diamScale, lambda, 1.0, 1.0,
                          new IDPairRangeAll(), "Bulk")));
    } else {
      Sim->units.setUnitLength(diamScale);
      Sim->interactions.push_back(shared_ptr<Interaction>(
          new IHardSphere(Sim, diamScale, new IDPairRangeAll(), "Bulk")));
    }

    Sim->addSpecies(shared_ptr<Species>(
        new SpPoint(Sim, new IDRangeAll(Sim), 1.0, "Bulk", 0)));

    Sim->topology.push_back(shared_ptr<Topology>(new TChain(Sim, 1, "Ring")));
    Sim->topology.back()->addMolecule(new IDRangeAll(Sim));

    Sim->BCs = shared_ptr<BoundaryCondition>(new BCNone(Sim));

    unsigned long nParticles = 0;

    Sim->particles.reserve(latticeSites.size());
    for (const Vector &position : latticeSites)
      Sim->particles.push_back(Particle(
          position, getRandVelVec() * Sim->units.unitVelocity(), nParticles++));
    break;
  }
  case 8: {
    // Pack of binary hard spheres
    if (vm.count("help")) {
      std::cout << "Mode specific options:\n"
                   "  8:  Binary Hard Spheres\n"
                   "       --i1 : Picks the packing routine to use [0] "
                   "(0:FCC,1:BCC,2:SC)\n"
                   "       --f1 : Size Ratio (B/A), must be (0,1] [0.1]\n"
                   "       --f2 : Mass Ratio (B/A) [0.001]\n"
                   "       --i2 : Number of large particles [100]\n";
      exit(1);
    }
    // Pack the system, determine the number of particles
    std::unique_ptr<UCell> packptr(
        new CURandomise(standardPackingHelper(new UParticle())));

    packptr->initialise();

    std::vector<Vector> latticeSites(packptr->placeObjects(Vector{0, 0, 0}));

    double massFrac = 0.001, sizeRatio = 0.1;
    size_t Na = 100;

    if (vm.count("f1"))
      sizeRatio = vm["f1"].as<double>();

    if (vm.count("f2"))
      massFrac = vm["f2"].as<double>();

    if (vm.count("i2"))
      Na = vm["i2"].as<size_t>();

    Sim->primaryCellSize = packptr->systemDims();

    double simVol = 1.0;

    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      simVol *= Sim->primaryCellSize[iDim];

    double particleDiam =
        pow(simVol * vm["density"].as<double>() / latticeSites.size(),
            double(1.0 / 3.0));

    if (Na > latticeSites.size())
      M_throw() << "Too many large particles for the selected packing";

    Sim->interactions.push_back(shared_ptr<Interaction>(new IHardSphere(
        Sim, particleDiam, new IDPairRangeSingle(new IDRangeRange(0, Na - 1)),
        "AAInt")));

    Sim->interactions.push_back(shared_ptr<Interaction>(new IHardSphere(
        Sim, ((1.0 + sizeRatio) / 2.0) * particleDiam,
        new IDPairRangePair(new IDRangeRange(0, Na - 1),
                            new IDRangeRange(Na, latticeSites.size() - 1)),
        "ABInt")));

    Sim->interactions.push_back(shared_ptr<Interaction>(new IHardSphere(
        Sim, sizeRatio * particleDiam, new IDPairRangeAll(), "BBInt")));

    Sim->addSpecies(shared_ptr<Species>(
        new SpPoint(Sim, new IDRangeRange(0, Na - 1), 1.0, "A", 0)));
    Sim->addSpecies(shared_ptr<Species>(new SpPoint(
        Sim, new IDRangeRange(Na, latticeSites.size() - 1), massFrac, "B", 0)));

    Sim->units.setUnitLength(particleDiam);

    unsigned long nParticles = 0;
    Sim->particles.reserve(latticeSites.size());
    for (const Vector &position : latticeSites)
      Sim->particles.push_back(Particle(
          position, getRandVelVec() * Sim->units.unitVelocity(), nParticles++));
    break;
  }
  case 9: {
    if (vm.count("help")) {
      std::cout << "Mode specific options:\n"
                   "  9:  Hard needle system\n"
                   "       --f1 : Inelasticity [1.0]\n"
                   "       --f2 : Inertia multiplicative factor [1.0]\n";
      exit(1);
    }
    // Pack of lines
    // Pack the system, determine the number of particles
    CURandom packroutine(vm["NCells"].as<unsigned long>(), Vector{1, 1, 1},
                         new UParticle());

    packroutine.initialise();

    std::vector<Vector> latticeSites(packroutine.placeObjects(Vector{0, 0, 0}));

    double particleDiam = pow(vm["density"].as<double>() / latticeSites.size(),
                              double(1.0 / 3.0));

    // Set up a standard simulation
    // We pick a scheduler algorithm based on the density of the system
    if (vm["density"].as<double>() * 8.0 >= vm["NCells"].as<unsigned long>())
      M_throw() << "Unable to simulate systems where box volume is <= (2L)^3";

    double elasticity = (vm.count("f1")) ? vm["f1"].as<double>() : 1.0;

    Sim->interactions.push_back(shared_ptr<Interaction>(new ILines(
        Sim, particleDiam, elasticity, new IDPairRangeAll(), "Bulk")));

    double inertiaMultiplicativeFactor =
        (vm.count("f2")) ? vm["f2"].as<double>() : 1.0;

    Sim->addSpecies(shared_ptr<Species>(new SpSphericalTop(
        Sim, new IDRangeAll(Sim), 1.0, "Bulk", 0,
        inertiaMultiplicativeFactor * particleDiam * particleDiam / 12.0)));

    Sim->units.setUnitLength(particleDiam);

    unsigned long nParticles = 0;
    Sim->particles.reserve(latticeSites.size());
    for (const Vector &position : latticeSites)
      Sim->particles.push_back(Particle(
          position, getRandVelVec() * Sim->units.unitVelocity(), nParticles++));
    Sim->dynamics->initOrientations();
    break;
  }
  case 10: {
    if (vm.count("help")) {
      std::cout << "Mode specific options:\n"
                   "  10: Monocomponent hard spheres using DSMC interactions\n"
                   "       --i1 : Picks the packing routine to use [0] "
                   "(0:FCC,1:BCC,2:SC)\n";
      exit(1);
    }
    // Pack of DSMC hard spheres
    // Pack the system, determine the number of particles
    std::unique_ptr<UCell> packptr(standardPackingHelper(new UParticle()));
    packptr->initialise();

    std::vector<Vector> latticeSites(packptr->placeObjects(Vector{0, 0, 0}));

    Sim->primaryCellSize = packptr->systemDims();

    double simVol = 1.0;

    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      simVol *= Sim->primaryCellSize[iDim];

    double particleDiam =
        pow(simVol * vm["density"].as<double>() / latticeSites.size(),
            double(1.0 / 3.0));

    Sim->units.setUnitLength(particleDiam);

    // This is to stop interactions being used for these particles
    Sim->interactions.push_back(shared_ptr<Interaction>(
        new INull(Sim, new IDPairRangeAll(), "Catchall")));

    // This is to provide data on the particles
    Sim->interactions.push_back(shared_ptr<Interaction>(
        new IHardSphere(Sim, particleDiam, new IDPairRangeAll(), "Bulk")));

    // double chi = 1.0 /
    //(4.0 * tij * vm["density"].as<double>() * std::sqrt(M_PI));

    double packfrac = vm["density"].as<double>() * M_PI / 6.0;
    double chi = (1.0 - 0.5 * packfrac) / std::pow(1.0 - packfrac, 3);
    double tij =
        1.0 / (4.0 * std::sqrt(M_PI) * vm["density"].as<double>() * chi);

    // No thermostat added yet
    Sim->systems.push_back(shared_ptr<System>(new SysDSMCSpheres(
        Sim, particleDiam, 2.0 * tij / latticeSites.size(), chi, 1.0,
        "Thermostat", new IDRangeAll(Sim), new IDRangeAll(Sim))));

    Sim->addSpecies(shared_ptr<Species>(
        new SpPoint(Sim, new IDRangeAll(Sim), 1.0, "Bulk", 0)));

    unsigned long nParticles = 0;
    Sim->particles.reserve(latticeSites.size());
    for (const Vector &position : latticeSites)
      Sim->particles.push_back(Particle(
          position, getRandVelVec() * Sim->units.unitVelocity(), nParticles++));
    break;
  }
  case 11: {
    M_throw() << "Option no longer supported";
  }
  case 12: {
    // Pack of DSMC hard sphere mixture
    if (vm.count("help")) {
      std::cout
          << "Mode specific options:\n"
             "  12: Binary hard spheres using DSMC interactions\n"
             "       --i1 : Picks the packing routine to use [0] "
             "(0:FCC,1:BCC,2:SC)\n"
             "       --i2 : Picks the g(r) to use (0:BMCSL, 1:VS, 2:HC2)\n"
             "       --f1 : Size Ratio (B/A), must be (0,1] [0.1]\n"
             "       --f2 : Mass Ratio (B/A) [0.001]\n"
             "       --f3 : Mol Fraction of large system (A) [0.95]\n";
      exit(1);
    }
    // Pack the system, determine the number of particles
    std::unique_ptr<UCell> packptr(
        new CURandomise(standardPackingHelper(new UParticle())));
    packptr->initialise();

    std::vector<Vector> latticeSites(packptr->placeObjects(Vector{0, 0, 0}));

    Sim->primaryCellSize = packptr->systemDims();

    double molFrac = 0.01, massFrac = 0.001, sizeRatio = 0.1;

    if (vm.count("f1"))
      sizeRatio = vm["f1"].as<double>();

    if (vm.count("f2"))
      massFrac = vm["f2"].as<double>();

    if (vm.count("f3"))
      molFrac = vm["f3"].as<double>();

    double simVol = 1.0;

    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      simVol *= Sim->primaryCellSize[iDim];

    double particleDiam =
        pow(simVol * vm["density"].as<double>() / latticeSites.size(),
            double(1.0 / 3.0));

    Sim->units.setUnitLength(particleDiam);

    // This is to stop interactions being used for these particles
    Sim->interactions.push_back(shared_ptr<Interaction>(
        new INull(Sim, new IDPairRangeAll(), "Catchall")));

    size_t nA = static_cast<size_t>(molFrac * latticeSites.size());

    double chiAA(1.0), chiAB(1.0), chiBB(1.0);

    size_t chimode(0);
    if (vm.count("i2"))
      chimode = vm["i2"].as<size_t>();

    double xi1 = (1.0 / 6.0) * M_PI * vm["density"].as<double>() *
                 (molFrac + (1.0 - molFrac) * sizeRatio);

    double xi2 = (1.0 / 6.0) * M_PI * vm["density"].as<double>() *
                 (molFrac + (1.0 - molFrac) * sizeRatio * sizeRatio);

    double xi3 =
        (1.0 / 6.0) * M_PI * vm["density"].as<double>() *
        (molFrac + (1.0 - molFrac) * sizeRatio * sizeRatio * sizeRatio);

    switch (chimode) {
    case 0:
      // BMCSL
      chiAA =
          (1.0 / (1.0 - xi3)) * (1.0 + 3.0 * xi2 / (2.0 * (1.0 - xi3)) +
                                 xi2 * xi2 / (2.0 * (1.0 - xi3) * (1.0 - xi3)));

      chiAB = (1.0 / (1.0 - xi3)) *
              (1.0 +
               3.0 * xi2 / (2.0 * (1.0 - xi3)) * sizeRatio /
                   (0.5 + 0.5 * sizeRatio) +
               xi2 * xi2 * std::pow(sizeRatio / (0.5 + 0.5 * sizeRatio), 2) /
                   (2.0 * (1.0 - xi3) * (1.0 - xi3)));

      chiBB = (1.0 / (1.0 - xi3)) *
              (1.0 + 3.0 * xi2 / (2.0 * (1.0 - xi3)) * sizeRatio +
               xi2 * xi2 * sizeRatio * sizeRatio /
                   (2.0 * (1.0 - xi3) * (1.0 - xi3)));
      break;
    case 1:
      // VS
      chiAA = (1.0 / (1.0 - xi3)) +
              (3.0 - xi3 + xi3 * xi3 * 0.5) * xi2 /
                  (2.0 * (1.0 - xi3) * (1.0 - xi3)) +
              (2.0 - xi3 - xi3 * xi3 * 0.5) * (2.0 * xi2 * xi2 + xi1 * xi3) /
                  (6.0 * std::pow(1.0 - xi3, 3.0));

      chiAB = (1.0 / (1.0 - xi3)) +
              (3.0 - xi3 + xi3 * xi3 * 0.5) * xi2 *
                  (sizeRatio / (0.5 + 0.5 * sizeRatio)) /
                  (2.0 * (1.0 - xi3) * (1.0 - xi3)) +
              (2.0 - xi3 - xi3 * xi3 * 0.5) * (2.0 * xi2 * xi2 + xi1 * xi3) *
                  (sizeRatio / (0.5 + 0.5 * sizeRatio)) *
                  (sizeRatio / (0.5 + 0.5 * sizeRatio)) /
                  (6.0 * std::pow(1.0 - xi3, 3.0));

      chiBB = (1.0 / (1.0 - xi3)) +
              (3.0 - xi3 + xi3 * xi3 * 0.5) * xi2 * sizeRatio /
                  (2.0 * (1.0 - xi3) * (1.0 - xi3)) +
              (2.0 - xi3 - xi3 * xi3 * 0.5) * (2.0 * xi2 * xi2 + xi1 * xi3) *
                  sizeRatio * sizeRatio / (6.0 * std::pow(1.0 - xi3, 3.0));
      break;
    case 2: {
      // HC2
      double x(3.0 * (xi2 - xi3) * 0.5);

      double R(1.0 / sizeRatio);

      chiAA = (1.0 / (1.0 - xi3)) +
              (3.0 - xi3 + xi3 * xi3 * 0.5) * xi2 /
                  (2.0 * (1.0 - xi3) * (1.0 - xi3)) +
              (2.0 - xi3 - xi3 * xi3 * 0.5) * (2.0 * xi2 * xi2 + xi1 * xi3) /
                  (6.0 * std::pow(1.0 - xi3, 3.0)) +
              std::exp(x) - 1.0 - x - x * x * 0.5;

      chiAB =
          (1.0 / (1.0 - xi3)) *
              (1.0 +
               3.0 * xi2 / (2.0 * (1.0 - xi3)) * sizeRatio /
                   (0.5 + 0.5 * sizeRatio) +
               xi2 * xi2 * std::pow(sizeRatio / (0.5 + 0.5 * sizeRatio), 2) /
                   (2.0 * (1.0 - xi3) * (1.0 - xi3))) +
          xi2 * xi2 * sizeRatio * sizeRatio * (R * R - 1.0) /
              (std::pow(1.0 - xi3, 3.0) * (R + 1.0) * (R + 1.0)) -
          xi2 * xi2 * xi2 * sizeRatio * sizeRatio * sizeRatio *
              (R * R * R - 1.0) /
              ((1.0 - xi3) * (1.0 - xi3) * (1.0 - xi3) * (R + 1.0) * (R + 1.0) *
               (R + 1.0));

      chiBB = (1.0 / (1.0 - xi3)) *
              (1.0 + 3.0 * xi2 / (2.0 * (1.0 - xi3)) * sizeRatio +
               xi2 * xi2 * sizeRatio * sizeRatio /
                   (2.0 * (1.0 - xi3) * (1.0 - xi3)));
    } break;
    default:
      M_throw() << "Unknown mode to set the chi's";
    }

    chiAB *= 2.0;

    double tAA = std::sqrt(M_PI) /
                 (chiAA * 4.0 * M_PI * molFrac * vm["density"].as<double>());

    double tAB =
        std::sqrt(2.0 * M_PI * massFrac / (1.0 + massFrac)) /
        (chiAB * 4.0 * M_PI * (1.0 - molFrac) * vm["density"].as<double>() *
         (0.5 + 0.5 * sizeRatio) * (0.5 + 0.5 * sizeRatio));

    double tBB = std::sqrt(M_PI * massFrac) /
                 (chiBB * 4.0 * M_PI * (1.0 - molFrac) *
                  vm["density"].as<double>() * sizeRatio * sizeRatio);

    // This is to provide data on the particles
    Sim->interactions.push_back(shared_ptr<Interaction>(new IHardSphere(
        Sim, particleDiam, new IDPairRangeSingle(new IDRangeRange(0, nA - 1)),
        "AAInt")));

    Sim->interactions.push_back(shared_ptr<Interaction>(new IHardSphere(
        Sim, sizeRatio * particleDiam,
        new IDPairRangeSingle(new IDRangeRange(nA, latticeSites.size() - 1)),
        "BBInt")));

    Sim->systems.push_back(shared_ptr<System>(new SysDSMCSpheres(
        Sim, particleDiam, tAA / (2.0 * nA), chiAA, 1.0, "AADSMC",
        new IDRangeRange(0, nA - 1), new IDRangeRange(0, nA - 1))));

    Sim->systems.push_back(shared_ptr<System>(new SysDSMCSpheres(
        Sim, ((1.0 + sizeRatio) / 2.0) * particleDiam, tAB / (2.0 * nA), chiAB,
        1.0, "ABDSMC", new IDRangeRange(0, nA - 1),
        new IDRangeRange(nA, latticeSites.size() - 1))));

    Sim->systems.push_back(shared_ptr<System>(new SysDSMCSpheres(
        Sim, sizeRatio * particleDiam, tBB / (2.0 * (latticeSites.size() - nA)),
        chiBB, 1.0, "BBDSMC", new IDRangeRange(nA, latticeSites.size() - 1),
        new IDRangeRange(nA, latticeSites.size() - 1))));

    Sim->addSpecies(shared_ptr<Species>(
        new SpPoint(Sim, new IDRangeRange(0, nA - 1), 1.0, "A", 0)));

    Sim->addSpecies(shared_ptr<Species>(new SpPoint(
        Sim, new IDRangeRange(nA, latticeSites.size() - 1), massFrac, "B", 0)));

    unsigned long nParticles = 0;
    Sim->particles.reserve(latticeSites.size());
    for (const Vector &position : latticeSites)
      Sim->particles.push_back(Particle(
          position, getRandVelVec() * Sim->units.unitVelocity(), nParticles++));
    break;
  }
  case 13: {
    // Pack of lines
    if (vm.count("help")) {
      std::cout << "  13: Crystal pack of sheared lines\n"
                   "       --i1 : Picks the packing routine to use [0] "
                   "(0:FCC,1:BCC,2:SC)\n"
                   "       --f1 : Inelasticity [1.0]\n";
      exit(1);
    }
    // Pack the system, determine the number of particles
    CURandom packroutine(vm["NCells"].as<unsigned long>(), Vector{1, 1, 1},
                         new UParticle());
    packroutine.initialise();
    std::vector<Vector> latticeSites(packroutine.placeObjects(Vector{0, 0, 0}));
    Sim->BCs = shared_ptr<BoundaryCondition>(new BCLeesEdwards(Sim));
    double particleDiam = pow(vm["density"].as<double>() / latticeSites.size(),
                              double(1.0 / 3.0));
    double elasticity = (vm.count("f1")) ? vm["f1"].as<double>() : 1.0;
    Sim->interactions.push_back(shared_ptr<Interaction>(new ILines(
        Sim, particleDiam, elasticity, new IDPairRangeAll(), "Bulk")));
    Sim->addSpecies(shared_ptr<Species>(
        new SpSphericalTop(Sim, new IDRangeAll(Sim), 1.0, "Bulk", 0,
                           particleDiam * particleDiam / 12.0)));
    Sim->units.setUnitLength(particleDiam);

    unsigned long nParticles = 0;
    Sim->particles.reserve(latticeSites.size());
    for (const Vector &position : latticeSites)
      Sim->particles.push_back(Particle(
          position, getRandVelVec() * Sim->units.unitVelocity(), nParticles++));

    Sim->dynamics->initOrientations();
    break;
  }
  case 14: {
    // Pack of Mings system
    if (vm.count("help")) {
      std::cout << "Mode specific options:\n"
                   "  14: Packing of spheres and linear rods made from stiff "
                   "polymers\n"
                   "       --i1 : Picks the packing routine to use [0] "
                   "(0:FCC,1:BCC,2:SC)\n"
                   "       --i2 : Number of spheres in chain\n"
                   "       --f1 : Mol fraction of spheres [0.5]\n"
                   "       --f2 : Rod Length [1.0]\n";
      exit(1);
    }
    // Pack the system, determine the number of particles
    double molfrac(0.5), massFrac(1.0), rodlength(1.0);
    size_t chainlength(10);
    size_t nPart;

    if (vm.count("f1"))
      molfrac = vm["f1"].as<double>();

    if (vm.count("f2"))
      rodlength = vm["f2"].as<double>();

    if (vm.count("i2"))
      chainlength = vm["i2"].as<size_t>();

    {
      std::unique_ptr<UCell> packptr(
          new CURandomise(standardPackingHelper(new UParticle())));

      packptr->initialise();

      std::vector<Vector> latticeSites(packptr->placeObjects(Vector{0, 0, 0}));

      nPart = latticeSites.size();

      Sim->primaryCellSize = packptr->systemDims();
    }

    size_t nPartA = size_t(nPart * molfrac);

    double simVol = 1.0;

    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      simVol *= Sim->primaryCellSize[iDim];

    double particleDiam =
        pow(simVol * vm["density"].as<double>() / nPart, 1.0 / 3.0);

    double particleDiamB = rodlength * particleDiam / chainlength;

    std::unique_ptr<UCell> packptr(standardPackingHelper(new CUBinary(
        nPartA, new UParticle(),
        new CUlinearRod(chainlength, 1.05 * particleDiamB, new UParticle()))));

    packptr->initialise();

    std::vector<Vector> latticeSites(packptr->placeObjects(Vector{0, 0, 0}));

    Sim->interactions.push_back(shared_ptr<Interaction>(new IHardSphere(
        Sim, particleDiam,
        new IDPairRangeSingle(new IDRangeRange(0, nPartA - 1)), "AAInt")));

    Sim->interactions.push_back(shared_ptr<Interaction>(new IHardSphere(
        Sim, (particleDiam + particleDiamB) / 2.0,
        new IDPairRangePair(new IDRangeRange(0, nPartA - 1),
                            new IDRangeRange(nPartA, latticeSites.size() - 1)),
        "ABInt")));

    Sim->interactions.push_back(shared_ptr<Interaction>(new ISquareBond(
        Sim, 0.9 * particleDiamB, 1.1 / 0.9, 1.0,
        new IDPairRangeChains(nPartA, latticeSites.size() - 1, chainlength),
        "Bonds")));

    Sim->interactions.push_back(shared_ptr<Interaction>(new IHardSphere(
        Sim, (chainlength - 1) * particleDiamB,
        new IDPairRangeChainEnds(nPartA, latticeSites.size() - 1, chainlength),
        "RodEnds")));

    Sim->interactions.push_back(shared_ptr<Interaction>(
        new IHardSphere(Sim, particleDiamB, new IDPairRangeAll(), "BBInt")));
    Sim->addSpecies(shared_ptr<Species>(
        new SpPoint(Sim, new IDRangeRange(0, nPartA - 1), 1.0, "A", 0)));

    Sim->addSpecies(shared_ptr<Species>(
        new SpPoint(Sim, new IDRangeRange(nPartA, latticeSites.size() - 1),
                    massFrac / chainlength, "B", 0)));

    Sim->units.setUnitLength(particleDiam);
    unsigned long nParticles = 0;
    Sim->particles.reserve(latticeSites.size());
    for (const Vector &position : latticeSites)
      Sim->particles.push_back(Particle(
          position, getRandVelVec() * Sim->units.unitVelocity(), nParticles++));
    break;
  }
  case 15: {
    // Pack of hard spheres
    if (vm.count("help")) {
      std::cout
          << "Mode specific options:\n"
             "  15: Monocomponent hard-parallel cubes\n"
             "       --i1 : Picks the packing routine to use [0] "
             "(0:FCC,1:BCC,2:SC)\n"
             "       --b1 : If set it enables the single occupancy model\n";
      exit(1);
    }
    // Pack the system, determine the number of particles

    if (!vm.count("i1") || vm["i1"].as<size_t>() != 2)
      M_throw()
          << "You should initialise cubes with simple cubic packing \"--i1 2\"";

    std::unique_ptr<UCell> packptr(standardPackingHelper(new UParticle()));
    packptr->initialise();

    std::vector<Vector> latticeSites(packptr->placeObjects(Vector{0, 0, 0}));

    if (latticeSites.size() % 2)
      M_throw() << "To make sure the system has zero momentum and +-1 "
                   "velocities, you must"
                   " use an even number of particles";

    Sim->primaryCellSize = packptr->systemDims();

    double simVol = 1.0;

    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      simVol *= Sim->primaryCellSize[iDim];

    double particleDiam =
        pow(simVol * vm["density"].as<double>() / latticeSites.size(),
            double(1.0 / 3.0));

    // Set up a standard simulation
    if (vm.count("b1"))
      Sim->globals.push_back(shared_ptr<Global>(new GSOCells(Sim, "SOCells")));

    Sim->interactions.push_back(shared_ptr<Interaction>(new IParallelCubes(
        Sim, particleDiam, 1.0, new IDPairRangeAll(), "Bulk")));

    Sim->addSpecies(shared_ptr<Species>(
        new SpPoint(Sim, new IDRangeAll(Sim), 1.0, "Bulk", 0)));
    Sim->units.setUnitLength(particleDiam);

    size_t nParticles = 0;
    Sim->particles.reserve(latticeSites.size());
    for (const Vector &position : latticeSites)
      Sim->particles.push_back(
          Particle(position,
                   Vector{Sim->units.unitVelocity(), Sim->units.unitVelocity(),
                          Sim->units.unitVelocity()},
                   nParticles++));

    {
      std::uniform_real_distribution<> uniform_dist(-0.5, 0.5);

      std::array<long, 3> tmp = getCells();

      Vector wobblespacing;

      for (size_t iDim(0); iDim < NDIM; ++iDim)
        wobblespacing[iDim] =
            (Sim->primaryCellSize[iDim] - particleDiam * tmp[iDim]) / tmp[iDim];

      for (Particle &part : Sim->particles)
        for (size_t iDim(0); iDim < NDIM; ++iDim)
          part.getPosition()[iDim] +=
              uniform_dist(Sim->ranGenerator) * wobblespacing[iDim];
    }

    {
      std::uniform_int_distribution<size_t> uniform_dist(0, nParticles - 1);

      size_t ID(uniform_dist(Sim->ranGenerator));

      for (size_t iDim(0); iDim < NDIM; ++iDim)
        for (size_t i(0); i < nParticles / 2; ++i) {
          while (Sim->particles[ID].getVelocity()[iDim] < 0)
            ID = uniform_dist(Sim->ranGenerator);

          Sim->particles[ID].getVelocity()[iDim] = -Sim->units.unitVelocity();
        }
    }
    break;
  }
  case 16: {
    // Pack of Lennard Jones stepped molecules
    if (vm.count("help")) {
      std::cout
          << "Mode specific options:\n"
             "  16: Stepped Potential, default is a Lennard-Jones potential by "
             "Chapela et al.\n"
             "       --i1 : Picks the packing routine to use [0] "
             "(0:FCC,1:BCC,2:SC)\n"
             "       --i2 : Set the potential type [Manual] (0: Manual, "
             "1:Lennard-Jones)\n"
             "       --i3 : Step placement algorithm [Energetic] (0:Radial, "
             "1:Energetic, 2: Volumetric)\n"
             "       --i4 : Step energy algorithm [Volumetric] (0:Mid-point, "
             "1:Left, 2:Right, 3:Volume, 4:Virial, 5:MidVolume)\n"
             "       --f1 : Order of approximation (Nsteps) [5.8]\n"
             "       --f2 : Cut-off [3]\n"
             "       --s1 : Manual entry of the potential (e.g., r1,E1:r2,E2) "
             "[default, Chapela et al potential 6]\n";
      exit(1);
    }

    // Pack the system, determine the number of particles
    std::unique_ptr<UCell> packptr(standardPackingHelper(new UParticle()));
    packptr->initialise();
    std::vector<Vector> latticeSites(packptr->placeObjects(Vector{0, 0, 0}));

    Sim->primaryCellSize = packptr->systemDims();

    double simVol = 1.0;
    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      simVol *= Sim->primaryCellSize[iDim];

    double particleDiam =
        pow(simVol * vm["density"].as<double>() / latticeSites.size(),
            double(1.0 / 3.0));

    bool twoD = false;
    if (vm.count("rectangular-box") &&
        (vm.count("i1") && vm["i1"].as<size_t>() == 2)) {
      std::array<long, 3> cells = getCells();
      if ((cells[0] == 1) || (cells[1] == 1) || (cells[2] == 1)) {
        twoD = true;
        derr << "Warning! Now assuming that you're trying to set up a 2D "
                "simulation!\n"
                "I'm going to temporarily calculate the density by the 2D "
                "definition!"
             << std::endl;

        size_t dimension;
        if (cells[0] == 1)
          dimension = 0;
        if (cells[1] == 1)
          dimension = 1;
        if (cells[2] == 1)
          dimension = 2;

        particleDiam =
            std::sqrt(simVol * vm["density"].as<double>() /
                      (Sim->primaryCellSize[dimension] * latticeSites.size()));

        dout << "I'm changing what looks like the unused box dimension ("
             << dimension
             << ") to the smallest value allowed by the neighbourlist "
             << "implementation (slightly more than 4 particle diameters)"
             << std::endl;

        Sim->primaryCellSize[dimension] = 4.0000001 * particleDiam;
      }
    }

    Sim->units.setUnitLength(particleDiam);
    Sim->units.setUnitTime(particleDiam);

    int potential_mode = 0;
    if (vm.count("i2"))
      potential_mode = vm["i2"].as<size_t>();

    int placement_mode = PotentialLennardJones::DELTAU;
    if (vm.count("i3"))
      placement_mode = vm["i3"].as<size_t>();

    int energy_mode = PotentialLennardJones::VOLUME;
    if (vm.count("i4"))
      energy_mode = vm["i4"].as<size_t>();

    double Nsteps = 5.8;
    if (vm.count("f1"))
      Nsteps = vm["f1"].as<double>();

    double cutoff = 3;
    if (vm.count("f2"))
      cutoff = vm["f2"].as<double>();

    shared_ptr<Potential> potential;
    switch (potential_mode) {
    case 0: // Manual entry of the potential
    {
      typedef std::pair<double, double> locpair;
      std::vector<locpair> diamvec;

      if (vm.count("s1")) {
        typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
        tokenizer steps(vm["s1"].as<std::string>(),
                        boost::char_separator<char>(":"));
        for (const std::string &step : steps) {
          tokenizer stepData(step, boost::char_separator<char>(","));
          tokenizer::iterator value_iter = stepData.begin();
          locpair data;
          try {
            data.first = boost::lexical_cast<double>(*value_iter);
            if (++value_iter == stepData.end())
              throw std::runtime_error("No comma");
            data.second = boost::lexical_cast<double>(*(value_iter));
            diamvec.push_back(data);

            if (++value_iter != stepData.end())
              throw std::runtime_error("Too many comma's");
          } catch (std::exception &except) {
            M_throw() << "Malformed step data, \"" << step << "\"\n"
                      << except.what();
          }
        }
      } else {
        diamvec.push_back(locpair(2.30, -0.06));
        diamvec.push_back(locpair(1.75, -0.22));
        diamvec.push_back(locpair(1.45, -0.55));
        diamvec.push_back(locpair(1.25, -0.98));
        diamvec.push_back(locpair(1.05, -0.47));
        diamvec.push_back(locpair(1.00, 0.76));
        diamvec.push_back(locpair(0.95, 3.81));
        diamvec.push_back(locpair(0.90, 10.95));
        diamvec.push_back(locpair(0.85, 27.55));
        diamvec.push_back(locpair(0.80, 66.74));
        diamvec.push_back(locpair(0.75, 1e300));
      }

      dout << "Building stepped potential" << std::endl;
      double oldr = std::numeric_limits<float>::infinity();
      for (locpair &p : diamvec) {
        dout << "Step r=" << p.first << ", E=" << p.second << std::endl;
        if (p.first > oldr)
          M_throw() << "Steps must be in descending order! r=" << p.first
                    << " is greater than old r=" << oldr;
        oldr = p.first;
      }

      potential.reset(new PotentialStepped(diamvec, false));
    } break;
    case 1: // Lennard-Jones potential
    {
      double kT = 1;
      if (energy_mode == PotentialLennardJones::VIRIAL) {
        if (!vm.count("thermostat"))
          M_throw()
              << "When using virial step placement, you must specify a "
                 "thermostat temperature using the --thermostat,-T option.";

        kT = vm["thermostat"].as<double>();
      }

      potential.reset(new PotentialLennardJones(1.0, 1.0, cutoff, energy_mode,
                                                placement_mode, Nsteps, kT));
    } break;
    }

    Sim->interactions.push_back(shared_ptr<Interaction>(new IStepped(
        Sim, potential, new IDPairRangeAll(), "Bulk", particleDiam, 1.0)));
    Sim->addSpecies(shared_ptr<Species>(
        new SpPoint(Sim, new IDRangeAll(Sim), 1.0, "Bulk", 0)));
    unsigned long nParticles = 0;
    Sim->particles.reserve(latticeSites.size());
    for (const Vector &position : latticeSites) {
      Sim->particles.push_back(Particle(
          position, getRandVelVec() * Sim->units.unitVelocity(), nParticles++));
      if (twoD)
        Sim->particles.back().getVelocity()[2] = 0;
    }
    break;
  }
  case 17:
  case 18: {
    M_throw() << "Option no longer supported";
  }
  case 19: {
    if (vm.count("help")) {
      std::cout << "Mode specific options:\n"
                   "  19: Oscillating plates bounding a system\n"
                   "       --b1 : Makes the particle collisions not affect the "
                   "plate\n"
                   "       --i1 : Picks the packing routine to use [0] "
                   "(0:FCC,1:BCC,2:SC)\n"
                   "       --i2 : Upper limit on the particles inserted [All]\n"
                   "       --f1 : Box to total particle mass ratio [2.93]\n"
                   "       --f2 : Length in particle radii [4]\n"
                   "       --f3 : Box frequency [1.23]\n"
                   "       --f4 : Initial displacement [10.7]\n"
                   "       --f5 : Particle-Particle inelasticity [0.75]\n"
                   "       --f6 : Particle-Wall inelasticity [0.76]\n"
                   "       --f7 : Cross section length [5.2]\n";
      exit(1);
    }

    double MassRatio = 2.93;
    if (vm.count("f1"))
      MassRatio = vm["f1"].as<double>();

    double L = 4.0;
    if (vm.count("f2"))
      L = vm["f2"].as<double>();
    // The available area to place the particle centers is actually -1 particle
    // diameter to the actual area
    L -= 1;

    double Omega0 = 1.23 * M_PI * 2.0;
    if (vm.count("f3"))
      Omega0 = vm["f3"].as<double>() * M_PI * 2.0;

    double Delta = 10.7;
    if (vm.count("f4"))
      Delta = vm["f4"].as<double>();

    double ParticleInelas = 0.75;
    if (vm.count("f5"))
      ParticleInelas = vm["f5"].as<double>();

    double PlateInelas = 0.76;
    if (vm.count("f6"))
      PlateInelas = vm["f6"].as<double>();

    double xy = 5.2;
    if (vm.count("f7"))
      xy = vm["f7"].as<double>();
    // The available area to place the particle centers is actually -1 particle
    // diameter to the actual area
    xy -= 1;

    // the  2.0 * L is to give an extra half box width on each side of the sim,
    // boxL is used as our unit length from now on.
    double boxL = 2.0 * L + 2.0 * Delta;
    double Aspect = xy / boxL;

    // Our simulation is set to be boxL X (1.1 * xy) X (1.1 * xy)
    Sim->primaryCellSize = Vector{1, 1.1 * Aspect, 1.1 * Aspect};

    //	Vector particleArea = Vector(0.5 * (L-2.0 * Sigma) / L ,
    //				     0.9 * Aspect, 0.9 * Aspect);
    //	//The minus one half spaces the particles off the wall
    //	Vector particleCOM = Vector(-(0.25 * (L - 2.0 * Sigma) + Delta - 0.5)/L,
    //				    0, 0);

    // The area in which we can place particle centers
    Vector particleArea = Vector{L / boxL, xy / boxL, xy / boxL};

    // The system starts at a full extention
    Vector particleCOM = Vector{Delta / boxL, 0, 0};

    UCell *sysPack;
    if (!vm.count("i1"))
      sysPack = new CUFCC(getCells(), particleArea, new UParticle());
    else
      switch (vm["i1"].as<size_t>()) {
      case 0: {
        sysPack = new CUFCC(getCells(), particleArea, new UParticle());
        break;
      }
      case 1: {
        sysPack = new CUBCC(getCells(), particleArea, new UParticle());
        break;
      }
      case 2: {
        sysPack = new CUSC(getCells(), particleArea, new UParticle());
        break;
      }
      default:
        M_throw() << "Not a valid packing type (--i1)";
      }

    std::unique_ptr<UCell> packptr(sysPack);
    packptr->initialise();

    std::vector<Vector> latticeSites(packptr->placeObjects(particleCOM));

    Sim->BCs = shared_ptr<BoundaryCondition>(new BCNone(Sim));

    // double simVol = 1.0;
    // for (size_t iDim = 0; iDim < NDIM; ++iDim)
    //   simVol *= Sim->primaryCellSize[iDim];

    double particleDiam = 1.0 / boxL;
    Sim->interactions.push_back(shared_ptr<Interaction>(new IHardSphere(
        Sim, particleDiam, ParticleInelas, new IDPairRangeAll(), "Bulk")));

    Sim->locals.push_back(shared_ptr<Local>(
        new LWall(Sim, PlateInelas, particleDiam, Vector{0, 0, 1},
                  Vector{0, 0, -0.5 * Aspect - 0.5 * particleDiam}, "Plate2",
                  new IDRangeAll(Sim))));

    Sim->locals.push_back(shared_ptr<Local>(
        new LWall(Sim, PlateInelas, particleDiam, Vector{0, 0, -1},
                  Vector{0, 0, +0.5 * Aspect + 0.5 * particleDiam}, "Plate3",
                  new IDRangeAll(Sim))));

    Sim->locals.push_back(shared_ptr<Local>(
        new LWall(Sim, PlateInelas, particleDiam, Vector{0, +1, 0},
                  Vector{0, -0.5 * Aspect - 0.5 * particleDiam, 0}, "Plate4",
                  new IDRangeAll(Sim))));

    Sim->locals.push_back(shared_ptr<Local>(
        new LWall(Sim, PlateInelas, particleDiam, Vector{0, -1, 0},
                  Vector{0, +0.5 * Aspect + 0.5 * particleDiam, 0}, "Plate5",
                  new IDRangeAll(Sim))));

    Sim->addSpecies(shared_ptr<Species>(
        new SpPoint(Sim, new IDRangeAll(Sim), 1.0, "Bulk", 0)));
    Sim->units.setUnitLength(particleDiam);

    size_t maxPart;
    if (vm.count("i2"))
      maxPart = vm["i2"].as<size_t>();
    else
      maxPart = latticeSites.size();

    unsigned long nParticles = 0;
    Sim->particles.reserve(maxPart);

    std::sort(latticeSites.begin(), latticeSites.end(), mySortPredictate);

    for (size_t i(0); i < maxPart; ++i)
      Sim->particles.push_back(
          Particle(latticeSites[i], getRandVelVec() * Sim->units.unitVelocity(),
                   nParticles++));

    bool strongPlate = false;
    if (vm.count("b1"))
      strongPlate = true;

    Sim->locals.push_back(shared_ptr<Local>(new LOscillatingPlate(
        Sim, Vector{0, 0, 0}, Vector{1, 0, 0}, Omega0, 0.5 * L / boxL,
        PlateInelas, Delta / boxL, MassRatio * nParticles, "Plate1",
        new IDRangeAll(Sim), 0.0, strongPlate)));
    break;
  }
  case 20: {
    // Pack of hard spheres then check overlaps against a set of triangles
    if (vm.count("help")) {
      std::cout << "Mode specific options:\n"
                   "  20: Load a set of triangles and plate it with spheres\n"
                   "       --i1 : Picks the packing routine to use [0] "
                   "(0:FCC,1:BCC,2:SC)\n"
                   "       --s1 : File name to load the triangles from\n"
                   "       --f1 : Size scale factor of the spheres when "
                   "checking for overlaps with triangles [1 = no scaling]\n";
      exit(1);
    }

    // Pack the system, determine the number of particles
    std::unique_ptr<UCell> packptr(standardPackingHelper(new UParticle()));
    size_t N = packptr->placeObjects(Vector{0, 0, 0}).size();
    Sim->primaryCellSize = packptr->systemDims();

    double simVol = 1.0;
    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      simVol *= Sim->primaryCellSize[iDim];

    double particleDiam =
        pow(simVol * vm["density"].as<double>() / N, double(1.0 / 3.0));
    double overlapDiameter = particleDiam;
    if (vm.count("f1"))
      overlapDiameter *= vm["f1"].as<double>();

    if (!vm.count("s1"))
      M_throw() << "No triangle file name specified";

    packptr.reset(
        new CUTriangleIntersect(standardPackingHelper(new UParticle()),
                                overlapDiameter, vm["s1"].as<std::string>()));

    packptr->initialise();

    std::vector<Vector> latticeSites(packptr->placeObjects(Vector{0, 0, 0}));

    Sim->interactions.push_back(shared_ptr<Interaction>(
        new IHardSphere(Sim, particleDiam, new IDPairRangeAll(), "Bulk")));
    Sim->addSpecies(shared_ptr<Species>(
        new SpPoint(Sim, new IDRangeAll(Sim), 1.0, "Bulk", 0)));

    Sim->units.setUnitLength(particleDiam);

    unsigned long nParticles = 0;
    Sim->particles.reserve(latticeSites.size());
    for (const Vector &position : latticeSites)
      Sim->particles.push_back(Particle(
          position, getRandVelVec() * Sim->units.unitVelocity(), nParticles++));
    break;
  }
  case 21: {
    // Pack of hard spheres in a cylinder
    if (vm.count("help")) {
      std::cout << "Mode specific options:\n"
                   "  21: Pack a cylinder with spheres\n"
                   "       --i1 : Picks the packing routine to use [0] "
                   "(0:FCC,1:BCC,2:SC)\n"
                   "       --f1 : Length over diameter of the cylinder\n";
      exit(1);
    }
    // Pack the system, determine the number of particles
    std::unique_ptr<UCell> packptr(standardPackingHelper(new UParticle()));
    packptr->initialise();

    std::vector<Vector> latticeSites(
        packptr->placeObjects(Vector{0.0, 0.0, 0.0}));

    double LoverD = 1;
    if (vm.count("f1"))
      LoverD = vm["f1"].as<double>();

    Sim->primaryCellSize = Vector{1.0, 1.0, 1.0};

    double boxlimit;
    double cylRad = 0.5;
    if (LoverD < 1) {
      // D is unity

      // Check if the cylinder limits the sim box
      boxlimit = LoverD;
      if ((1.0 / std::sqrt(2.0)) < LoverD)
        boxlimit = (1.0 / std::sqrt(2.0));

      Sim->primaryCellSize[0] = LoverD;
    } else {
      // L is unity
      Sim->primaryCellSize[1] = 1.0 / LoverD;
      Sim->primaryCellSize[2] = 1.0 / LoverD;

      boxlimit = 1.0;

      cylRad = 0.5 / LoverD;

      if ((1.0 / (LoverD * std::sqrt(2.0))) < 1.0)
        boxlimit = (1.0 / (LoverD * std::sqrt(2.0)));
    }

    // Shrink the box a little more
    boxlimit *= 0.9;

    Sim->BCs = shared_ptr<BoundaryCondition>(new BCPeriodicXOnly(Sim));

    double particleDiam = pow(vm["density"].as<double>() / latticeSites.size(),
                              double(1.0 / 3.0)) *
                          boxlimit;

    Sim->locals.push_back(shared_ptr<Local>(
        new LCylinder(Sim, 1.0, particleDiam, Vector{1, 0, 0}, Vector{0, 0, 0},
                      -cylRad, "Cylinder", new IDRangeAll(Sim))));

    Sim->interactions.push_back(shared_ptr<Interaction>(
        new IHardSphere(Sim, particleDiam, new IDPairRangeAll(), "Bulk")));

    Sim->addSpecies(shared_ptr<Species>(
        new SpPoint(Sim, new IDRangeAll(Sim), 1.0, "Bulk", 0)));

    Sim->units.setUnitLength(particleDiam);

    unsigned long nParticles = 0;
    Sim->particles.reserve(latticeSites.size());
    for (const Vector &position : latticeSites)
      Sim->particles.push_back(
          Particle(position * boxlimit,
                   getRandVelVec() * Sim->units.unitVelocity(), nParticles++));

    Sim->ensemble.reset(new dynamo::EnsembleNVE(Sim));
    break;
  }
  case 22: {
    // Pack of hard spheres on a plate
    if (vm.count("help")) {
      std::cout << "Mode specific options:\n"
                   "  22: Infinite system with spheres falling onto a plate "
                   "with gravity\n"
                   "       --i1 : Picks the packing routine to use [0] "
                   "(0:FCC,1:BCC,2:SC)\n";
      exit(1);
    }
    // Pack the system, determine the number of particles
    std::unique_ptr<UCell> packptr(standardPackingHelper(new UParticle()));
    packptr->initialise();

    std::vector<Vector> latticeSites(packptr->placeObjects(Vector{0, 0, 0}));

    Sim->primaryCellSize = packptr->systemDims();

    Sim->BCs = shared_ptr<BoundaryCondition>(new BCNone(Sim));

    double simVol = 1.0;

    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      simVol *= Sim->primaryCellSize[iDim];

    double particleDiam =
        pow(simVol * vm["density"].as<double>() / latticeSites.size(),
            double(1.0 / 3.0));

    Sim->units.setUnitLength(particleDiam);
    Sim->dynamics = shared_ptr<Dynamics>(
        new DynGravity(Sim, Vector{0, -Sim->units.unitAcceleration(), 0}));

    double elasticity = 1.0;

    if (vm.count("f1"))
      elasticity = vm["f1"].as<double>();

    Sim->interactions.push_back(shared_ptr<Interaction>(new IHardSphere(
        Sim, particleDiam, elasticity, new IDPairRangeAll(), "Bulk")));

    Sim->addSpecies(shared_ptr<Species>(
        new SpPoint(Sim, new IDRangeAll(Sim), 1.0, "Bulk", 0)));

    // We actually shrink our lattice length scale by 0.999 and our
    // wall spacing by 0.9995 to prevent particles being
    // initialised touching the wall and to insert the wall just
    // inside the primary image
    Sim->locals.push_back(shared_ptr<Local>(new LWall(
        Sim, 1.0, particleDiam, Vector{0, 1, 0},
        Vector{0, -0.5 * Sim->primaryCellSize[1] - 0.5 * particleDiam, 0},
        "GroundPlate", new IDRangeAll(Sim))));

    unsigned long nParticles = 0;
    Sim->particles.reserve(latticeSites.size());
    for (const Vector &position : latticeSites)
      Sim->particles.push_back(
          Particle(0.999 * position,
                   getRandVelVec() * Sim->units.unitVelocity(), nParticles++));
    break;
  }
  case 23: {
    // Pack of hard spheres to form a funnel
    if (vm.count("help")) {
      std::cout
          << "Mode specific options:\n"
             "  23: Funnel test for static spheres in gravity\n"
             "       --i1 : Number of rows to remove when making the cone hole "
             "[3]\n"
             "       --f1 : Height of the cone in particle diameters [10]\n"
             "       --f2 : Max radius of the cone in particle diameters "
             "[7.5]\n"
             "       --f3 : Elasticity of the particles [0.4]\n";
      exit(1);
    }

    double H = 10;
    if (vm.count("f1"))
      H = vm["f1"].as<double>();

    double R = 7.5;
    if (vm.count("f2"))
      R = vm["f2"].as<double>();

    size_t rowskip = 3;
    if (vm.count("i1"))
      rowskip = vm["i1"].as<size_t>();

    double elasticity = 0.4;
    if (vm.count("f3"))
      elasticity = vm["f3"].as<double>();

    double Sv = 1.0; // Vertical spacing
    double Sr = 1.0; // Radial spacing
    const double elasticV = 1.0;

    Sim->primaryCellSize = Vector{1, 1, 1};

    double particleDiam = std::min(1 / (2 * R + 1), 1 / (H + 1));

    Sim->units.setUnitLength(particleDiam);
    Sim->dynamics = shared_ptr<Dynamics>(
        new DynGravity(Sim, Vector{0, -Sim->units.unitAcceleration(), 0},
                       elasticV * Sim->units.unitVelocity()));

    Sim->interactions.push_back(shared_ptr<Interaction>(new IHardSphere(
        Sim, particleDiam, elasticity, new IDPairRangeAll(), "Bulk")));

    /// Now build our funnel, so we know how many particles it takes
    std::vector<Vector> funnelSites;
    size_t Nv =
        static_cast<size_t>(std::sqrt(H * H + R * R) / Sv); // Number of circles
    double deltaZ = H / Nv;
    for (size_t circle(rowskip); circle <= Nv; ++circle) {
      double r = R * circle / Nv;
      size_t Nr = static_cast<size_t>(M_PI / std::asin(Sr / (2 * r)));
      double deltaPhi = 2 * M_PI / Nr;

      for (size_t radialstep(0); radialstep < Nr; ++radialstep)
        funnelSites.push_back(particleDiam *
                                  Vector{r * std::sin(radialstep * deltaPhi),
                                         circle * deltaZ,
                                         r * std::cos(radialstep * deltaPhi)} -
                              Vector{0, 0.5, 0});
    }

    for (size_t circle(0);
         particleDiam * ((circle + 1) * Sv + Nv * deltaZ - 0.5) - 0.5 < 0.4;
         ++circle) {
      double r = R;
      size_t Nr = static_cast<size_t>(M_PI / std::asin(Sr / (2 * r)));
      double deltaPhi = 2 * M_PI / Nr;

      for (size_t radialstep(0); radialstep < Nr; ++radialstep)
        funnelSites.push_back(particleDiam *
                                  Vector{r * std::sin(radialstep * deltaPhi),
                                         (circle + 1) * Sv + Nv * deltaZ,
                                         r * std::cos(radialstep * deltaPhi)} -
                              Vector{0, 0.5, 0});
    }

    // Build a list of the dynamic particles
    std::vector<Vector> dynamicSites;
    Sr = std::max(1.1, Sr); // Increase the spacing to a min of 1.1
    Sv = std::max(1.1, Sv); //....
    for (size_t circle(0);
         particleDiam * ((circle + 1) * Sv + Nv * deltaZ - 0.5) - 0.5 < 0.4;
         ++circle)
      for (double r(R - Sr); r > 0; r -= Sr) {
        size_t Nr = static_cast<size_t>(M_PI / std::asin(Sr / (2 * r)));
        double deltaPhi = 2 * M_PI / Nr;

        for (size_t radialstep(0); radialstep < Nr; ++radialstep)
          dynamicSites.push_back(
              particleDiam * Vector{r * std::sin(radialstep * deltaPhi),
                                    (circle + 1) * Sv + Nv * deltaZ,
                                    r * std::cos(radialstep * deltaPhi)} -
              Vector{0, 0.5, 0});
      }

    Sim->addSpecies(shared_ptr<Species>(
        new SpFixedCollider(Sim, new IDRangeRange(0, funnelSites.size() - 1),
                            "FunnelParticles", 0)));
    Sim->addSpecies(shared_ptr<Species>(new SpPoint(
        Sim,
        new IDRangeRange(funnelSites.size(),
                         funnelSites.size() + dynamicSites.size() - 1),
        1.0, "Bulk", 0)));

    unsigned long nParticles = 0;
    Sim->particles.reserve(funnelSites.size() + dynamicSites.size());

    for (const Vector &position : funnelSites)
      Sim->particles.push_back(
          Particle(position, Vector{0, 0, 0}, nParticles++));

    for (const Vector &position : dynamicSites) {
      Vector vel = getRandVelVec() * Sim->units.unitVelocity();
      if (vel[1] > 0)
        vel[1] = -vel[1]; // So particles don't fly out of the hopper
      Sim->particles.push_back(Particle(position, vel, nParticles++));
    }
    break;
  }
  case 24: {
    // An isolated MJ model polymer
    if (vm.count("help")) {
      std::cout << "Mode specific options:\n"
                   "  24: Random walk of an isolated MJ model polymer\n"
                   "      "
                   "(DOI:10.1002/"
                   "(SICI)1097-0134(19990101)34:1<49::AID-PROT5>3.0.CO;2-L)\n"
                   "       --f1 : Diameter [1.6]\n"
                   "       --f2 : Well width factor [1.5]\n"
                   "       --f3 : Bond inner core [0.9]\n"
                   "       --f4 : Bond outer well [1.1]\n"
                   "       --s1 : Sequence to use [GVGTGSGRGQGVGTGSGRGQ]\n";
      exit(1);
    }
    std::string stringseq = "GVGTGSGRGQGVGTGSGRGQ";
    if (vm.count("s1"))
      stringseq = vm["s1"].as<std::string>();

    size_t chainlength = stringseq.size();

    double sigmin(0.9), sigmax(1.1), sigma(1.6), lambda(1.5);

    if (vm.count("f1"))
      sigma = vm["f1"].as<double>();

    if (vm.count("f2"))
      lambda = vm["f2"].as<double>();

    if (vm.count("f3"))
      sigmin = vm["f3"].as<double>();

    if (vm.count("f4"))
      sigmax = vm["f4"].as<double>();

    // We need a box which is at least 4 times the maximum
    // interaction distance for the neighbourlist to function.
    double length1 = 4 * std::max(sigma * lambda, sigmax);
    // We also would want at least a big enough neighbourlist to
    // encompass the fully stretched out chain without wrapping
    // it around
    double length2 = chainlength * sigmax + sigma * lambda;

    // And we double it just to be sure
    double diamScale = 1.0 / (2 * std::max(length1, length2));

    // Sit the particles 95% away of max distance from each other
    // to help with seriously overlapping wells
    CURandWalk sysPack(chainlength,
                       (sigmin + 0.95 * (sigmax - sigmin)) * diamScale,
                       sigma * diamScale, new UParticle());

    sysPack.initialise();

    // Drop them in the middle of the sim
    std::vector<Vector> latticeSites(sysPack.placeObjects(Vector{0, 0, 0}));

    Sim->interactions.push_back(shared_ptr<Interaction>(new ISquareBond(
        Sim, sigmin * diamScale, sigmax / sigmin, 1.0,
        new IDPairRangeChains(0, latticeSites.size() - 1, latticeSites.size()),
        "Bonds")));

    {
      std::vector<size_t> seq;

      seq.resize(chainlength, 0);

      // initialize MJ interaction matrix
      std::map<std::string, double> MJinter;

      MJinter["GW"] = -0.25;
      MJinter["GV"] = -0.15;
      MJinter["GT"] = -0.04;
      MJinter["GS"] = -0.01;
      MJinter["GR"] = 0.09;
      MJinter["GQ"] = 0.13;
      MJinter["GP"] = 0.02;
      MJinter["GY"] = -0.22;
      MJinter["GG"] = -0.29;
      MJinter["GF"] = -0.19;
      MJinter["GE"] = 0.32;
      MJinter["GD"] = 0.11;
      MJinter["GC"] = -0.31;
      MJinter["GA"] = -0.08;
      MJinter["GN"] = -0.01;
      MJinter["GM"] = -0.17;
      MJinter["GL"] = -0.16;
      MJinter["GK"] = 0.29;
      MJinter["GI"] = -0.13;
      MJinter["GH"] = 0.00;
      MJinter["EN"] = 0.12;
      MJinter["ME"] = 0.12;
      MJinter["MD"] = 0.30;
      MJinter["MG"] = -0.17;
      MJinter["MF"] = -0.83;
      MJinter["MA"] = -0.27;
      MJinter["MC"] = -0.61;
      MJinter["MM"] = -0.70;
      MJinter["ML"] = -0.70;
      MJinter["MN"] = 0.04;
      MJinter["MI"] = -0.66;
      MJinter["MH"] = -0.29;
      MJinter["MK"] = 0.29;
      MJinter["MT"] = -0.11;
      MJinter["MW"] = -0.73;
      MJinter["MV"] = -0.51;
      MJinter["MQ"] = -0.06;
      MJinter["MP"] = -0.13;
      MJinter["MS"] = 0.05;
      MJinter["MR"] = 0.03;
      MJinter["MY"] = -0.56;
      MJinter["FP"] = -0.19;
      MJinter["FQ"] = -0.11;
      MJinter["FR"] = -0.05;
      MJinter["FS"] = -0.12;
      MJinter["FT"] = -0.15;
      MJinter["FV"] = -0.67;
      MJinter["FW"] = -0.68;
      MJinter["FY"] = -0.58;
      MJinter["FA"] = -0.36;
      MJinter["FC"] = -0.67;
      MJinter["FD"] = 0.18;
      MJinter["FE"] = 0.14;
      MJinter["FF"] = -0.88;
      MJinter["FG"] = -0.19;
      MJinter["FH"] = -0.34;
      MJinter["FI"] = -0.73;
      MJinter["FK"] = 0.19;
      MJinter["FL"] = -0.80;
      MJinter["FM"] = -0.83;
      MJinter["FN"] = -0.01;
      MJinter["SY"] = -0.08;
      MJinter["SS"] = 0.05;
      MJinter["SR"] = 0.16;
      MJinter["SQ"] = 0.22;
      MJinter["SP"] = 0.20;
      MJinter["SW"] = -0.01;
      MJinter["SV"] = 0.04;
      MJinter["ST"] = 0.04;
      MJinter["SK"] = 0.36;
      MJinter["SI"] = 0.03;
      MJinter["SH"] = 0.04;
      MJinter["SN"] = 0.09;
      MJinter["SM"] = 0.05;
      MJinter["SL"] = -0.02;
      MJinter["SC"] = -0.13;
      MJinter["SA"] = 0.10;
      MJinter["SG"] = -0.01;
      MJinter["SF"] = -0.12;
      MJinter["SE"] = 0.18;
      MJinter["SD"] = 0.10;
      MJinter["YI"] = -0.49;
      MJinter["YH"] = -0.30;
      MJinter["YK"] = -0.05;
      MJinter["YM"] = -0.56;
      MJinter["YL"] = -0.55;
      MJinter["YN"] = -0.11;
      MJinter["YA"] = -0.20;
      MJinter["YC"] = -0.39;
      MJinter["YE"] = -0.08;
      MJinter["YD"] = -0.07;
      MJinter["YG"] = -0.22;
      MJinter["YF"] = -0.58;
      MJinter["YY"] = -0.45;
      MJinter["YQ"] = -0.14;
      MJinter["YP"] = -0.25;
      MJinter["YS"] = -0.08;
      MJinter["YR"] = -0.25;
      MJinter["YT"] = -0.09;
      MJinter["YW"] = -0.49;
      MJinter["YV"] = -0.38;
      MJinter["LF"] = -0.80;
      MJinter["LG"] = -0.16;
      MJinter["LD"] = 0.27;
      MJinter["LE"] = 0.17;
      MJinter["LC"] = -0.65;
      MJinter["LA"] = -0.38;
      MJinter["LN"] = 0.04;
      MJinter["LL"] = -0.84;
      MJinter["LM"] = -0.70;
      MJinter["LK"] = 0.22;
      MJinter["LH"] = -0.18;
      MJinter["LI"] = -0.81;
      MJinter["LV"] = -0.74;
      MJinter["LW"] = -0.62;
      MJinter["LT"] = -0.15;
      MJinter["LR"] = -0.04;
      MJinter["LS"] = -0.02;
      MJinter["LP"] = -0.12;
      MJinter["LQ"] = -0.04;
      MJinter["LY"] = -0.55;
      MJinter["RT"] = 0.11;
      MJinter["RV"] = 0.08;
      MJinter["RW"] = -0.21;
      MJinter["RP"] = 0.17;
      MJinter["RQ"] = 0.09;
      MJinter["RR"] = 0.19;
      MJinter["RS"] = 0.16;
      MJinter["RY"] = -0.25;
      MJinter["RD"] = -0.24;
      MJinter["RE"] = -0.22;
      MJinter["RF"] = -0.05;
      MJinter["RG"] = 0.09;
      MJinter["RA"] = 0.24;
      MJinter["RC"] = 0.08;
      MJinter["RL"] = -0.04;
      MJinter["RM"] = 0.03;
      MJinter["RN"] = 0.10;
      MJinter["RH"] = 0.05;
      MJinter["RI"] = 0.00;
      MJinter["RK"] = 0.66;
      MJinter["VH"] = -0.06;
      MJinter["VI"] = -0.67;
      MJinter["EM"] = 0.12;
      MJinter["EL"] = 0.17;
      MJinter["IR"] = 0.00;
      MJinter["EI"] = 0.17;
      MJinter["EH"] = 0.00;
      MJinter["EK"] = -0.06;
      MJinter["EE"] = 0.46;
      MJinter["ED"] = 0.44;
      MJinter["EG"] = 0.32;
      MJinter["EF"] = 0.14;
      MJinter["EA"] = 0.38;
      MJinter["EC"] = 0.20;
      MJinter["VM"] = -0.51;
      MJinter["EY"] = -0.08;
      MJinter["IW"] = -0.60;
      MJinter["ET"] = 0.16;
      MJinter["EW"] = -0.00;
      MJinter["EV"] = 0.26;
      MJinter["EQ"] = 0.27;
      MJinter["EP"] = 0.37;
      MJinter["ES"] = 0.18;
      MJinter["ER"] = -0.22;
      MJinter["II"] = -0.74;
      MJinter["IH"] = -0.13;
      MJinter["IK"] = 0.24;
      MJinter["IM"] = -0.66;
      MJinter["IN"] = 0.14;
      MJinter["KC"] = 0.33;
      MJinter["KA"] = 0.41;
      MJinter["KG"] = 0.29;
      MJinter["KF"] = 0.19;
      MJinter["KE"] = -0.06;
      MJinter["KD"] = -0.01;
      MJinter["KK"] = 0.76;
      MJinter["KI"] = 0.24;
      MJinter["KH"] = 0.38;
      MJinter["KN"] = 0.22;
      MJinter["KM"] = 0.29;
      MJinter["KL"] = 0.22;
      MJinter["KS"] = 0.36;
      MJinter["KR"] = 0.66;
      MJinter["KQ"] = 0.28;
      MJinter["KP"] = 0.47;
      MJinter["KW"] = 0.09;
      MJinter["KV"] = 0.29;
      MJinter["KT"] = 0.33;
      MJinter["KY"] = -0.05;
      MJinter["DN"] = 0.02;
      MJinter["DL"] = 0.27;
      MJinter["DM"] = 0.30;
      MJinter["DK"] = -0.01;
      MJinter["DH"] = -0.10;
      MJinter["DI"] = 0.22;
      MJinter["DF"] = 0.18;
      MJinter["DG"] = 0.11;
      MJinter["DD"] = 0.29;
      MJinter["DE"] = 0.44;
      MJinter["DC"] = 0.12;
      MJinter["DA"] = 0.27;
      MJinter["DY"] = -0.07;
      MJinter["DV"] = 0.36;
      MJinter["DW"] = 0.07;
      MJinter["DT"] = 0.11;
      MJinter["DR"] = -0.24;
      MJinter["DS"] = 0.10;
      MJinter["DP"] = 0.33;
      MJinter["DQ"] = 0.24;
      MJinter["QQ"] = 0.20;
      MJinter["QP"] = 0.17;
      MJinter["QS"] = 0.22;
      MJinter["QR"] = 0.09;
      MJinter["QT"] = 0.12;
      MJinter["QW"] = -0.02;
      MJinter["QV"] = 0.08;
      MJinter["QY"] = -0.14;
      MJinter["QA"] = 0.22;
      MJinter["QC"] = -0.07;
      MJinter["QE"] = 0.27;
      MJinter["QD"] = 0.24;
      MJinter["QG"] = 0.13;
      MJinter["QF"] = -0.11;
      MJinter["QI"] = -0.01;
      MJinter["QH"] = 0.15;
      MJinter["QK"] = 0.28;
      MJinter["QM"] = -0.06;
      MJinter["QL"] = -0.04;
      MJinter["QN"] = 0.06;
      MJinter["WG"] = -0.25;
      MJinter["WF"] = -0.68;
      MJinter["WE"] = -0.00;
      MJinter["WD"] = 0.07;
      MJinter["WC"] = -0.66;
      MJinter["WA"] = -0.27;
      MJinter["WN"] = -0.10;
      MJinter["WM"] = -0.73;
      MJinter["WL"] = -0.62;
      MJinter["WK"] = 0.09;
      MJinter["WI"] = -0.60;
      MJinter["WH"] = -0.37;
      MJinter["WW"] = -0.64;
      MJinter["WV"] = -0.51;
      MJinter["WT"] = -0.02;
      MJinter["WS"] = -0.01;
      MJinter["WR"] = -0.21;
      MJinter["WQ"] = -0.02;
      MJinter["WP"] = -0.37;
      MJinter["WY"] = -0.49;
      MJinter["PR"] = 0.17;
      MJinter["PS"] = 0.20;
      MJinter["PP"] = 0.11;
      MJinter["PQ"] = 0.17;
      MJinter["PV"] = -0.05;
      MJinter["PW"] = -0.37;
      MJinter["PT"] = 0.13;
      MJinter["PY"] = -0.25;
      MJinter["PC"] = -0.18;
      MJinter["PA"] = 0.15;
      MJinter["PF"] = -0.19;
      MJinter["PG"] = 0.02;
      MJinter["PD"] = 0.33;
      MJinter["PE"] = 0.37;
      MJinter["PK"] = 0.47;
      MJinter["PH"] = 0.01;
      MJinter["PI"] = -0.05;
      MJinter["PN"] = 0.18;
      MJinter["PL"] = -0.12;
      MJinter["PM"] = -0.13;
      MJinter["CK"] = 0.33;
      MJinter["CI"] = -0.64;
      MJinter["CH"] = -0.36;
      MJinter["CN"] = -0.01;
      MJinter["CM"] = -0.61;
      MJinter["CL"] = -0.65;
      MJinter["CC"] = -1.19;
      MJinter["CA"] = -0.33;
      MJinter["CG"] = -0.31;
      MJinter["CF"] = -0.67;
      MJinter["CE"] = 0.20;
      MJinter["CD"] = 0.12;
      MJinter["CY"] = -0.39;
      MJinter["CS"] = -0.13;
      MJinter["CR"] = 0.08;
      MJinter["CQ"] = -0.07;
      MJinter["CP"] = -0.18;
      MJinter["CW"] = -0.66;
      MJinter["CV"] = -0.59;
      MJinter["CT"] = -0.15;
      MJinter["IY"] = -0.49;
      MJinter["VA"] = -0.32;
      MJinter["VC"] = -0.59;
      MJinter["VD"] = 0.36;
      MJinter["VE"] = 0.26;
      MJinter["VF"] = -0.67;
      MJinter["VG"] = -0.15;
      MJinter["IQ"] = -0.01;
      MJinter["IP"] = -0.05;
      MJinter["IS"] = 0.03;
      MJinter["VK"] = 0.29;
      MJinter["VL"] = -0.74;
      MJinter["IT"] = -0.15;
      MJinter["VN"] = 0.12;
      MJinter["IV"] = -0.67;
      MJinter["VP"] = -0.05;
      MJinter["VQ"] = 0.08;
      MJinter["VR"] = 0.08;
      MJinter["VS"] = 0.04;
      MJinter["VT"] = -0.07;
      MJinter["IL"] = -0.81;
      MJinter["VV"] = -0.65;
      MJinter["VW"] = -0.51;
      MJinter["IA"] = -0.37;
      MJinter["VY"] = -0.38;
      MJinter["IC"] = -0.64;
      MJinter["IE"] = 0.17;
      MJinter["ID"] = 0.22;
      MJinter["IG"] = -0.13;
      MJinter["IF"] = -0.73;
      MJinter["HY"] = -0.30;
      MJinter["HR"] = 0.05;
      MJinter["HS"] = 0.04;
      MJinter["HP"] = 0.01;
      MJinter["HQ"] = 0.15;
      MJinter["HV"] = -0.06;
      MJinter["HW"] = -0.37;
      MJinter["HT"] = -0.03;
      MJinter["HK"] = 0.38;
      MJinter["HH"] = -0.40;
      MJinter["HI"] = -0.13;
      MJinter["HN"] = 0.00;
      MJinter["HL"] = -0.18;
      MJinter["HM"] = -0.29;
      MJinter["HC"] = -0.36;
      MJinter["HA"] = 0.07;
      MJinter["HF"] = -0.34;
      MJinter["HG"] = 0.00;
      MJinter["HD"] = -0.10;
      MJinter["HE"] = 0.00;
      MJinter["NH"] = 0.00;
      MJinter["NI"] = 0.14;
      MJinter["NK"] = 0.22;
      MJinter["NL"] = 0.04;
      MJinter["NM"] = 0.04;
      MJinter["NN"] = -0.06;
      MJinter["NA"] = 0.15;
      MJinter["NC"] = -0.01;
      MJinter["ND"] = 0.02;
      MJinter["NE"] = 0.12;
      MJinter["NF"] = -0.01;
      MJinter["NG"] = -0.01;
      MJinter["NY"] = -0.11;
      MJinter["NP"] = 0.18;
      MJinter["NQ"] = 0.06;
      MJinter["NR"] = 0.10;
      MJinter["NS"] = 0.09;
      MJinter["NT"] = 0.04;
      MJinter["NV"] = 0.12;
      MJinter["NW"] = -0.10;
      MJinter["TY"] = -0.09;
      MJinter["TV"] = -0.07;
      MJinter["TW"] = -0.02;
      MJinter["TT"] = 0.03;
      MJinter["TR"] = 0.11;
      MJinter["TS"] = 0.04;
      MJinter["TP"] = 0.13;
      MJinter["TQ"] = 0.12;
      MJinter["TN"] = 0.04;
      MJinter["TL"] = -0.15;
      MJinter["TM"] = -0.11;
      MJinter["TK"] = 0.33;
      MJinter["TH"] = -0.03;
      MJinter["TI"] = -0.15;
      MJinter["TF"] = -0.15;
      MJinter["TG"] = -0.04;
      MJinter["TD"] = 0.11;
      MJinter["TE"] = 0.16;
      MJinter["TC"] = -0.15;
      MJinter["TA"] = 0.04;
      MJinter["AA"] = -0.12;
      MJinter["AC"] = -0.33;
      MJinter["AE"] = 0.38;
      MJinter["AD"] = 0.27;
      MJinter["AG"] = -0.08;
      MJinter["AF"] = -0.36;
      MJinter["AI"] = -0.37;
      MJinter["AH"] = 0.07;
      MJinter["AK"] = 0.41;
      MJinter["AM"] = -0.27;
      MJinter["AL"] = -0.38;
      MJinter["AN"] = 0.15;
      MJinter["AQ"] = 0.22;
      MJinter["AP"] = 0.15;
      MJinter["AS"] = 0.10;
      MJinter["AR"] = 0.24;
      MJinter["AT"] = 0.04;
      MJinter["AW"] = -0.27;
      MJinter["AV"] = -0.32;
      MJinter["AY"] = -0.20;

      // Transcribe the sequence
      std::cout << std::endl;
      std::cout << "chainlength=" << stringseq.size() << std::endl;
      size_t type_int = 0;
      std::string type_string;
      std::map<std::string, size_t> mapping;
      std::map<std::string, size_t>::iterator it;

      // translate letters to numbers
      for (size_t i = 0; i < chainlength; ++i) {
        type_string = stringseq.at(i);
        it = mapping.find(type_string);
        if (it == mapping.end()) {
          mapping[type_string] = type_int;
          ++type_int;
        }
        seq[i] = mapping[type_string];
      }
      for (const auto &e1 : mapping) {
        std::cout << e1.first << "  " << e1.second << std::endl;
      }
      std::cout << "protein sequence:" << std::endl;
      for (size_t i = 0; i < chainlength; ++i) {
        std::cout << i << "  " << seq[i] << std::endl;
      }

      Sim->interactions.push_back(shared_ptr<Interaction>(
          new ISWSequence(Sim, sigma * diamScale, lambda, 1.0, seq,
                          new IDPairRangeAll(), "Bulk")));

      ISWSequence &interaction(
          dynamic_cast<ISWSequence &>(*Sim->interactions["Bulk"]));

      //	    interaction.getAlphabet().at(0).at(0) = 1.0;
      //	    interaction.getAlphabet().at(1).at(0) = 0.5;
      //	    interaction.getAlphabet().at(0).at(1) = 0.5;

      // set interaction matrix
      std::map<std::string, double>::iterator it_MJ;
      std::string pair;
      for (const auto &e1 : mapping) {
        for (const auto &e2 : mapping) {
          pair = e1.first + e2.first;
          it_MJ = MJinter.find(pair);
          std::cout << e1.second << "  " << e2.second << "  " << e1.first
                    << e2.first << "  "
                    << pair
                    //<< MJinter[pair]
                    << std::endl;
          if (it_MJ == MJinter.end()) {
            M_throw() << "Entered a monomer not in the database.";
          }
          interaction.getAlphabet().at(e1.second).at(e2.second) =
              -MJinter[pair];
          pair = e2.first + e1.first;
          interaction.getAlphabet().at(e2.second).at(e1.second) =
              -MJinter[pair];
        }
      }
    }

    Sim->addSpecies(shared_ptr<Species>(
        new SpPoint(Sim, new IDRangeAll(Sim), 1.0, "Bulk", 0)));

    Sim->units.setUnitLength(diamScale);
    // Set the unit energy to 1 (assuming the unit of mass is 1);
    Sim->units.setUnitTime(diamScale);

    Sim->topology.push_back(
        shared_ptr<Topology>(new TChain(Sim, 1, "HelixPolymer")));

    Sim->topology.back()->addMolecule(new IDRangeAll(Sim));

    unsigned long nParticles = 0;

    Sim->BCs = shared_ptr<BoundaryCondition>(new BCNone(Sim));

    Sim->particles.reserve(latticeSites.size());
    for (const Vector &position : latticeSites)
      Sim->particles.push_back(Particle(
          position, getRandVelVec() * Sim->units.unitVelocity(), nParticles++));
    break;
  }
  case 25: {
    // Pack of hard spheres to form a funnel, slide and cup
    if (vm.count("help")) {
      std::cout << "Mode specific options:\n"
                   "  25: Funnel and cup simulation (with sleepy particles)\n"
                   "       --f1 : Elasticity [0.4]\n"
                   "       --f2 : Elastic Velocity [Disabled]\n"
                   "       --f3 : Sleep velocity [Disabled]\n"
                   "       --f4 : tc model time [0.1] (0=off)\n"
                   "       --f5 : If using a sleep velocity, this sets the "
                   "periodic wake up time [Disabled]\n";
      exit(1);
    }

    double elasticity = 0.4;
    if (vm.count("f1"))
      elasticity = vm["f1"].as<double>();

    double elasticV = 0;
    if (vm.count("f2"))
      elasticV = vm["f2"].as<double>();

    double sleepV = 0;
    if (vm.count("f3"))
      sleepV = vm["f3"].as<double>();

    double tc = 0.04;
    if (vm.count("f4"))
      tc = vm["f4"].as<double>();

    if (!tc)
      tc = -std::numeric_limits<float>::infinity();

    double wakeTime = 0;
    if (vm.count("f5"))
      wakeTime = vm["f5"].as<double>();

    Sim->primaryCellSize = Vector{1, 1, 1};

    double Rmax = 0.01999;
    double l = 4;
    double particleDiam = (2 * Rmax) / l;

    Sim->units.setUnitLength(particleDiam);

    Sim->dynamics = shared_ptr<Dynamics>(new DynGravity(
        Sim, Vector{0, -Sim->units.unitAcceleration(), 0},
        elasticV * Sim->units.unitVelocity(), tc * Sim->units.unitTime()));

    /// Now build our funnel, so we know how many particles it takes
    std::vector<Vector> funnelSites;
    Vector move({0, 0, -0.1});
    double factor = particleDiam / Rmax;
    double x, y, z;

    double spacing = 2.01 * particleDiam / factor;

    /// Build our funnel, so we know how many particles it takes
    double R = 0.3;
    double H = 0.34;
    size_t Nv = static_cast<size_t>(std::sqrt(H * H + R * R) /
                                    (spacing)); // Number of circles
    double deltaZ = H / Nv;
    for (size_t circle(3); circle <= Nv; ++circle) {
      double r = R * circle / Nv;
      size_t Nr = static_cast<size_t>(M_PI / std::asin(spacing / (2 * r)));
      double deltaPhi = 2 * M_PI / Nr;

      for (size_t radialstep(0); radialstep < Nr; ++radialstep)
        funnelSites.push_back(factor *
                                  Vector{r * std::sin(radialstep * deltaPhi),
                                         circle * deltaZ + 0.0052,
                                         r * std::cos(radialstep * deltaPhi)} -
                              move);
    }

    spacing = 2.1 * particleDiam / factor;

    // Slide
    for (int k = -1; k < 13; k++) {
      for (int i = -1; i < 10; i++) {
        x = cos(i * 2 * M_PI / 16) * 0.11;
        y = -sin(i * 2 * M_PI / 16) * 0.11 - k * 0.02 + 0.05;
        z = k * 0.02 * 2 - 0.1;
        funnelSites.push_back(factor * Vector{x, y, z} - move);
      }
    }

    // wall
    for (int k = -2; k < 2; k++) { // Slide blocking Wall
      for (int i = -2; i < 1; i++) {
        x = k * 2 * Rmax + 0.02;
        y = i * 2 * Rmax + 0.08;
        z = -0.18;
        funnelSites.push_back(factor * Vector{x, y, z} - move);
      }
    }

    // End container
    double r = 0.26;
    size_t Nr = static_cast<size_t>(M_PI / std::asin(spacing / (2 * r)));
    for (int k = 0; k < 10; k++) { // Box Walls
      {
        for (size_t i = 0; i < Nr; i++) {
          x = sin(i * 2 * M_PI / Nr) * r;
          y = -0.68 + k * 0.04001;
          z = -cos(i * 2 * M_PI / Nr) * r + 0.46;
          funnelSites.push_back(factor * Vector{x, y, z} - move);
        }
      }
    }

    for (int k = 10; k < 19; k++) { // Box Deflection wall, 1/4 is missing
      for (size_t i = Nr / 8; i < (Nr * 7) / 8; i++) {
        x = sin(i * 2 * M_PI / Nr) * r;
        y = -0.68 + k * 0.04001;
        z = -cos(i * 2 * M_PI / Nr) * r + 0.46;
        funnelSites.push_back(factor * Vector{x, y, z} - move);
      }
    }

    // Box bottom
    y = -0.72;
    funnelSites.push_back(factor * Vector{0, y, 0.46} - move);
    for (double lr = spacing; lr < r + spacing; lr += spacing) {
      size_t Nr = static_cast<size_t>(M_PI / std::asin(spacing / (2 * lr)));
      for (size_t i = 0; i < Nr; i++) {
        x = cos(i * 2 * M_PI / Nr) * lr;
        z = sin(i * 2 * M_PI / Nr) * lr + 0.46;
        funnelSites.push_back(factor * Vector{x, y, z} - move);
      }
    }

    for (int j = 1; j < 10; ++j)
      for (int i = 0; i < 46; i++) { // Funnel circular walls
        x = cos(i * 2 * M_PI / 46) * 0.30;
        z = sin(i * 2 * M_PI / 46) * 0.30;
        y = 0.35 + spacing * j;
        funnelSites.push_back(factor * Vector{x, y, z} - move);
      }

    // Clear out overlapping funnel particles
    for (std::vector<Vector>::iterator iPtr = funnelSites.begin();
         iPtr != funnelSites.end(); ++iPtr) {
      bool overlapping = false;
      for (std::vector<Vector>::iterator jPtr = iPtr + 1;
           jPtr != funnelSites.end(); ++jPtr) {
        Vector rij = (*iPtr) - (*jPtr);
        if (rij.nrm() < 2.0 * particleDiam) {
          overlapping = true;
          break;
        }
      }

      if (overlapping) {
        iPtr = funnelSites.erase(iPtr);
        if (iPtr == funnelSites.end())
          break;
      }
    }

    // Build a list of the dynamic particles
    std::vector<Vector> dynamicSites;

    for (double r = 0.30 - spacing; r > spacing; r -= spacing) {
      size_t Nr = static_cast<size_t>(M_PI / std::asin(spacing / (2 * r)));
      for (double y = 0.35 + spacing; y < 0.65; y += spacing)
        for (size_t i = 0; i < Nr; i++) {
          x = cos(i * 2 * M_PI / Nr) * r;
          z = sin(i * 2 * M_PI / Nr) * r;
          dynamicSites.push_back(factor * Vector{x, y, z} - move);
        }
    }
    //	for(int i=0;i<dim;i++){
    //	  for(int k=0;k<1.5*dim;k++){
    //	    for(int l=0;l<1.5*dim;l++){
    //	      double x=0.041*k-0.15;
    //	      double y=0.30+0.041*i;
    //	      double z=0.041*l-0.15;
    //	      dynamicSites.push_back(factor * Vector(x,y,z) - move);
    //	    }
    //	  }
    //	}

    Sim->interactions.push_back(shared_ptr<Interaction>(new IHardSphere(
        Sim, particleDiam * 2.0, elasticity, new IDPairRangeAll(), "Bulk")));

    Sim->addSpecies(shared_ptr<Species>(
        new SpFixedCollider(Sim, new IDRangeRange(0, funnelSites.size() - 1),
                            "FunnelParticles", 0)));
    Sim->addSpecies(shared_ptr<Species>(new SpPoint(
        Sim,
        new IDRangeRange(funnelSites.size(),
                         funnelSites.size() + dynamicSites.size() - 1),
        1.0, "Bulk", 0)));

    if (sleepV) {
      Sim->systems.push_back(shared_ptr<System>(new SSleep(
          Sim, "Sleeper",
          new IDRangeRange(funnelSites.size(),
                           funnelSites.size() + dynamicSites.size() - 1),
          sleepV * Sim->units.unitVelocity())));

      if (wakeTime)
        Sim->globals.push_back(shared_ptr<Global>(new GWaker(
            Sim, "Waker",
            new IDRangeRange(funnelSites.size(),
                             funnelSites.size() + dynamicSites.size() - 1),
            wakeTime * Sim->units.unitTime(),
            0.5 * sleepV * Sim->units.unitVelocity(), "SchedulerNBList")));
    }

    unsigned long nParticles = 0;
    Sim->particles.reserve(funnelSites.size() + dynamicSites.size());

    for (const Vector &position : funnelSites)
      Sim->particles.push_back(
          Particle(position, Vector{0, 0, 0}, nParticles++));

    for (const Vector &position : dynamicSites) {
      Vector vel = 0.001 * getRandVelVec() * Sim->units.unitVelocity();
      if (vel[1] > 0)
        vel[1] = -vel[1]; // So particles don't fly out of the hopper
      Sim->particles.push_back(Particle(position, vel, nParticles++));
    }
    break;
  }
  case 26: {
    // Pack of polydisperse sheared hard spheres
    if (vm.count("help")) {
      std::cout
          << "Mode specific options:\n"
             "  26: Polydisperse (Gaussian) hard spheres in LEBC (shearing)\n"
             "      Note: Generated particle diameters are restricted to the "
             "range (0,1].\n"
             "            Mass is distributed according to volume (constant "
             "density).\n"
             "            A particle with diameter of 1 has a mass of 1.\n"
             "       --i1 : Picks the packing routine to use [0] "
             "(0:FCC,1:BCC,2:SC)\n"
             "       --f1 : Inelasticity [1.0]\n"
             "       --f2 : Mean size [0.5]\n"
             "       --f3 : Standard deviation [0.1]\n";
      exit(1);
    }

    double mean = 0.5;
    if (vm.count("f2"))
      mean = vm["f2"].as<double>();

    double variance = 0.1;
    if (vm.count("f3"))
      variance = vm["f3"].as<double>();

    // FCC simple cubic pack of hard spheres with inelasticity and shearing
    // Pack the system, determine the number of particles
    std::unique_ptr<UCell> packptr(standardPackingHelper(new UParticle()));
    packptr->initialise();

    std::vector<Vector> latticeSites(packptr->placeObjects(Vector{0, 0, 0}));

    Sim->primaryCellSize = packptr->systemDims();

    double simVol = 1.0;

    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      simVol *= Sim->primaryCellSize[iDim];

    double particleDiam =
        pow(simVol * vm["density"].as<double>() / latticeSites.size(),
            double(1.0 / 3.0));

    bool twoD = false;
    if (vm.count("rectangular-box") &&
        (vm.count("i1") && vm["i1"].as<size_t>() == 2)) {
      std::array<long, 3> cells = getCells();
      if ((cells[0] == 1) || (cells[1] == 1) || (cells[2] == 1)) {
        twoD = true;
        derr << "Warning! Now assuming that you're trying to set up a 2D "
                "simulation!\n"
                "I'm going to temporarily calculate the density by the 2D "
                "definition!"
             << std::endl;

        size_t dimension;
        if (cells[0] == 1)
          dimension = 0;
        if (cells[1] == 1)
          dimension = 1;
        if (cells[2] == 1)
          dimension = 2;

        particleDiam =
            std::sqrt(simVol * vm["density"].as<double>() /
                      (Sim->primaryCellSize[dimension] * latticeSites.size()));

        dout << "I'm changing what looks like the unused box dimension ("
             << dimension
             << ") to the smallest value allowed by the neighbourlist "
                "implementation (slightly more than 4 particle diameters)"
             << std::endl;

        Sim->primaryCellSize[dimension] = 4.0000001 * particleDiam;
      }
    }

    double elasticity = 1.0;

    if (vm.count("f1"))
      elasticity = vm["f1"].as<double>();

    Sim->BCs = shared_ptr<BoundaryCondition>(new BCLeesEdwards(Sim));
    const double shearRate(1);

    shared_ptr<ParticleProperty> D(new ParticleProperty(
        latticeSites.size(), Property::Units::Length(), "D", particleDiam));

    shared_ptr<ParticleProperty> M(new ParticleProperty(
        latticeSites.size(), Property::Units::Mass(), "M", 1.0));
    Sim->_properties.push(D);
    Sim->_properties.push(M);

    std::normal_distribution<> normal_dist(mean, variance);

    for (size_t i(0); i < latticeSites.size(); ++i) {
      double diameter = normal_dist(Sim->ranGenerator);
      for (size_t attempt(0);
           ((diameter <= 0) || (diameter > 1)) && (attempt < 100); ++attempt)
        diameter = normal_dist(Sim->ranGenerator);

      if ((diameter <= 0) || (diameter > 1))
        M_throw() << "After 100 attempts, not a single valid particle diameter "
                     "could be generated."
                  << "Please recheck the distribution parameters";

      D->getProperty(i) = diameter * particleDiam;

      // A particle with unit diameter has unit mass
      double mass = diameter * diameter * diameter;

      M->getProperty(i) = mass;
    }

    Sim->interactions.push_back(shared_ptr<Interaction>(
        new IHardSphere(Sim, "D", elasticity, new IDPairRangeAll(), "Bulk")));

    Sim->addSpecies(shared_ptr<Species>(
        new SpPoint(Sim, new IDRangeAll(Sim), "M", "Bulk", 0)));

    Sim->units.setUnitLength(particleDiam);

    unsigned long nParticles = 0;
    Sim->particles.reserve(latticeSites.size());
    for (const Vector &position : latticeSites) {
      Sim->particles.push_back(Particle(
          position, getRandVelVec() * Sim->units.unitVelocity(), nParticles++));
      if (twoD)
        Sim->particles.back().getVelocity()[2] = 0;
    }

    // Insert a linear profile, zero momentum then add a vel gradient
    Sim->setCOMVelocity();
    for (Particle &part : Sim->particles)
      part.getVelocity()[0] += part.getPosition()[1] * shearRate;
    break;
  }
  case 27: {
    if (vm.count("help")) {
      std::cout << "  27: Crystal pack of snowmen molecules\n"
                   "       --i1 : Picks the packing routine to use [0] "
                   "(0:FCC,1:BCC,2:SC)\n"
                   "       --f1 : Inelasticity [1.0]\n"
                   "       --f2 : Size ratio [1.0]\n"
                   "       --f3 : Mass ratio [size_ratio^3]\n";
      exit(1);
    }

    const double elasticity = (vm.count("f1")) ? vm["f1"].as<double>() : 1.0;
    const double sizeratio = (vm.count("f2")) ? vm["f2"].as<double>() : 1.0;
    const double massratio = (vm.count("f3"))
                                 ? vm["f3"].as<double>()
                                 : sizeratio * sizeratio * sizeratio;

    std::unique_ptr<UCell> packptr(standardPackingHelper(new UParticle()));
    packptr->initialise();

    std::vector<Vector> latticeSites(packptr->placeObjects(Vector{0, 0, 0}));

    Sim->primaryCellSize = packptr->systemDims();

    double simVol = 1.0;

    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      simVol *= Sim->primaryCellSize[iDim];

    double sigmaA =
        pow(simVol * vm["density"].as<double>() / latticeSites.size(),
            double(1.0 / 3.0));

    bool twoD = false;
    size_t unusedDimension = 0;
    if (vm.count("rectangular-box") &&
        (vm.count("i1") && vm["i1"].as<size_t>() == 2)) {
      std::array<long, 3> cells = getCells();
      if ((cells[0] == 1) || (cells[1] == 1) || (cells[2] == 1)) {
        derr << "Warning! Now assuming that you're trying to set up a 2D "
                "simulation!\n"
                "I'm going to temporarily calculate the density by the 2D "
                "definition!"
             << std::endl;

        if (cells[0] == 1)
          unusedDimension = 0;
        else if (cells[1] == 1)
          unusedDimension = 1;
        else if (cells[2] == 1)
          unusedDimension = 2;
        else
          M_throw() << "Continuity error!";

        sigmaA = std::sqrt(
            simVol * vm["density"].as<double>() /
            (Sim->primaryCellSize[unusedDimension] * latticeSites.size()));

        dout << "I'm changing what looks like the unused box dimension ("
             << unusedDimension
             << ") to the smallest value allowed by the neighbourlist "
                "implementation (slightly more than 4 particle diameters)"
             << std::endl;
        Sim->primaryCellSize[unusedDimension] =
            10.0000001 * std::max(sigmaA, sigmaA * sizeratio);
        twoD = true;
      }
    }

    const double sigmaB = sigmaA * sizeratio;
    const double LB = (sigmaA + sigmaB) * 0.5 / (1 + massratio);
    const double LA = massratio * LB;
    const double mA = 1 / (1 + massratio);
    const double mB = mA * massratio;
    const double I = (mA * sigmaA * sigmaA + mB * sigmaB * sigmaB) * 0.1 +
                     mA * LA * LA + mB * LB * LB;

    shared_ptr<IDumbbells> interaction(
        new IDumbbells(Sim, elasticity, new IDPairRangeAll(), "Bulk"));
    Sim->interactions.push_back(interaction);
    interaction->addSphere(LA * Vector{0.0, 0.0, 1.0}, sigmaA);
    interaction->addSphere(-LB * Vector{0.0, 0.0, 1.0}, sigmaA * sizeratio);

    Sim->addSpecies(shared_ptr<Species>(
        new SpSphericalTop(Sim, new IDRangeAll(Sim), 1.0, "Bulk", 0, I)));

    Sim->units.setUnitLength(sigmaA);

    size_t nParticles = 0;
    Sim->particles.reserve(latticeSites.size());
    for (const Vector &position : latticeSites)
      Sim->particles.push_back(Particle(
          position, getRandVelVec() * Sim->units.unitVelocity(), nParticles++));

    Sim->dynamics->initOrientations();
    if (twoD) {
      Vector rotationAxis{0, 0, 0};
      rotationAxis[unusedDimension] = 1;
      std::normal_distribution<> dist(0, 1);
      for (size_t i(0); i < Sim->particles.size(); ++i) {
        Sim->particles[i].getVelocity()[2] = 0;
        auto &data = Sim->dynamics->getRotData(i);
        Vector orientation{0, 0, 0};
        orientation[(unusedDimension + 1) % 3] = dist(Sim->ranGenerator);
        orientation[(unusedDimension + 2) % 3] = dist(Sim->ranGenerator);
        data.orientation =
            magnet::math::Quaternion::fromToVector(orientation.normal());
        data.angularVelocity = rotationAxis * dist(Sim->ranGenerator);
      }
      static_pointer_cast<IDumbbells>(Sim->interactions["Bulk"])
          ->setUnusedDimension(unusedDimension);
    }
    break;
  }
  case 28: {
    if (vm.count("help")) {
      std::cout
          << "Mode specific options:\n"
             "  28: Rotating drum made out of particles.\n"
             "       --i1 : Depth of the drum in particle diameters [5]\n"
             "       --f1 : Radius of the drum in particle diameters (from "
             "centre to boundary particle centre) [7.5]\n"
             "       --f2 : Elasticity of the particles [0.4]\n"
             "       --f3 : Spacing of the particles in particle diameters "
             "[3]\n"
             "       --f4 : Incline of the system in degrees [6]\n"
             "       --f5 : Rotations per unit time [0.001]\n"
             "       --f6 : \"Steps\" per rotation [360]\n"
             "       --f7 : Elastic velocity [0.5]\n"
             "       --f8 : Tangential restitution coefficient [disabled]\n";
      exit(1);
    }

    size_t depth = 5;
    if (vm.count("i1"))
      depth = vm["i1"].as<size_t>();

    double R = 7.5;
    if (vm.count("f1"))
      R = vm["f1"].as<double>();

    double elasticity = 0.4;
    if (vm.count("f2"))
      elasticity = vm["f2"].as<double>();

    double dynamicSpacing = 3;
    if (vm.count("f3"))
      dynamicSpacing = vm["f3"].as<double>();

    double incline = 6;
    if (vm.count("f4"))
      incline = vm["f4"].as<double>();

    double RPT = 0.001;
    if (vm.count("f5"))
      RPT = vm["f5"].as<double>();

    double steps_per_rotation = 360;
    if (vm.count("f6"))
      steps_per_rotation = vm["f6"].as<double>();

    double elasticV = 0.5;
    if (vm.count("f7"))
      elasticV = vm["f7"].as<double>();

    double et = 1.0;
    if (vm.count("f8"))
      et = vm["f8"].as<double>();

    const double diameter = 1.0;
    const double g = 1.0;

    Sim->units.setUnitLength(diameter);

    Sim->primaryCellSize = Vector{2 * R + 1, 2 * R + 1, double(depth)};

    // Set up a standard simulation
    Sim->scheduler = shared_ptr<SNeighbourList>(
        new SNeighbourList(Sim, new CBTFEL<HeapPEL>()));

    incline *= M_PI / 180.0;
    Sim->dynamics = shared_ptr<Dynamics>(new DynGravity(
        Sim, g * Vector{0, -cos(incline), sin(incline)}, elasticV));

    if (et == 1.0)
      Sim->interactions.push_back(shared_ptr<Interaction>(new IHardSphere(
          Sim, diameter, elasticity, new IDPairRangeAll(), "Bulk")));
    else
      Sim->interactions.push_back(shared_ptr<Interaction>(new IHardSphere(
          Sim, diameter, elasticity, et, new IDPairRangeAll(), "Bulk")));

    Sim->systems.push_back(shared_ptr<System>(new SysRotateGravity(
        Sim, "GravityRotator", 1.0 / (RPT * steps_per_rotation), 2 * M_PI * RPT,
        Vector{0, 0, 1})));

    /// Now build our funnel, so we know how many particles it takes
    std::vector<Vector> funnelSites;

    for (size_t circle(0); circle < depth; ++circle) {
      const size_t Nr =
          static_cast<size_t>(M_PI / std::asin(diameter / (2 * R)));
      const double deltaPhi = 2 * M_PI / Nr;

      for (size_t radialstep(0); radialstep < Nr; ++radialstep)
        funnelSites.push_back(Vector{R * std::sin(radialstep * deltaPhi),
                                     R * std::cos(radialstep * deltaPhi),
                                     circle * diameter});
    }

    // Build a list of the dynamic particles
    std::vector<Vector> dynamicSites;
    for (double circlePos(0); circlePos < depth * diameter;
         circlePos += dynamicSpacing * diameter)
      for (double circleR = R - dynamicSpacing * diameter; circleR > diameter;
           circleR -= dynamicSpacing * diameter) {
        const size_t Nr =
            static_cast<size_t>(M_PI / std::asin(diameter / (2 * circleR)));
        const double deltaPhi = 2 * M_PI / Nr;

        for (size_t radialstep(0); radialstep < Nr; ++radialstep)
          dynamicSites.push_back(
              Vector{circleR * std::sin(radialstep * deltaPhi),
                     circleR * std::cos(radialstep * deltaPhi), circlePos});
      }

    Sim->addSpecies(shared_ptr<Species>(
        new SpFixedCollider(Sim, new IDRangeRange(0, funnelSites.size() - 1),
                            "FunnelParticles", 0)));

    if (et == 1.0)
      Sim->addSpecies(shared_ptr<Species>(new SpPoint(
          Sim,
          new IDRangeRange(funnelSites.size(),
                           funnelSites.size() + dynamicSites.size() - 1),
          1.0, "Bulk", 0)));
    else
      Sim->addSpecies(shared_ptr<Species>(new SpSphericalTop(
          Sim,
          new IDRangeRange(funnelSites.size(),
                           funnelSites.size() + dynamicSites.size() - 1),
          1.0, "Bulk", 0, diameter * diameter / 10.0)));

    unsigned long nParticles = 0;
    Sim->particles.reserve(funnelSites.size() + dynamicSites.size());

    for (const Vector &position : funnelSites)
      Sim->particles.push_back(
          Particle(position, Vector{0, 0, 0}, nParticles++));

    for (const Vector &position : dynamicSites) {
      Vector vel = getRandVelVec() * Sim->units.unitVelocity();
      if (vel[1] > 0)
        vel[1] = -vel[1]; // So particles don't fly out of the hopper
      Sim->particles.push_back(Particle(position, vel, nParticles++));
    }

    if (et != 1.0)
      Sim->dynamics->initOrientations(1);
    break;
  }
  default:
    M_throw() << "Did not recognise the packer mode you wanted";
  }

  Sim->ensemble = dynamo::Ensemble::loadEnsemble(*Sim);
}

Vector IPPacker::getNormalisedCellDimensions() {
  std::array<long, 3> cells = getCells();
  size_t maxdim = 0;

  // Determine the biggest dimension
  for (size_t iDim = 1; iDim < NDIM; ++iDim)
    if (cells[iDim] > cells[maxdim])
      maxdim = iDim;

  Vector retval;

  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    retval[iDim] =
        static_cast<double>(cells[iDim]) / static_cast<double>(cells[maxdim]);

  return retval;
}

UCell *IPPacker::standardPackingHelper(UCell *tmpPtr, bool forceRectangular) {
  UCell *sysPack;

  Vector boxDimensions{1, 1, 1};

  if (vm.count("rectangular-box") || forceRectangular) {
    boxDimensions = getNormalisedCellDimensions();
  }

  if (!vm.count("i1"))
    sysPack = new CUFCC(getCells(), boxDimensions, tmpPtr);
  else
    switch (vm["i1"].as<size_t>()) {
    case 0: {
      sysPack = new CUFCC(getCells(), boxDimensions, tmpPtr);
      break;
    }
    case 1: {
      sysPack = new CUBCC(getCells(), boxDimensions, tmpPtr);
      break;
    }
    case 2: {
      sysPack = new CUSC(getCells(), boxDimensions, tmpPtr);
      break;
    }
    case 3: {
      sysPack = new CUHCP(getCells(), boxDimensions, tmpPtr);
      break;
    }
    default:
      M_throw() << "Not a valid packing type (--i1)";
    }

  return sysPack;
}

std::array<long, 3> IPPacker::getCells() {
  long NCells = vm["NCells"].as<unsigned long>();
  std::array<long, 3> cells = {{NCells, NCells, NCells}};

  if (vm.count("xcell"))
    cells[0] = vm["xcell"].as<unsigned long>();

  if (vm.count("ycell"))
    cells[1] = vm["ycell"].as<unsigned long>();

  if (vm.count("zcell"))
    cells[2] = vm["zcell"].as<unsigned long>();

  return cells;
}

Vector IPPacker::getRandVelVec() {
  // See http://mathworld.wolfram.com/SpherePointPicking.html
  std::normal_distribution<> normal_dist(0.0, (1.0 / sqrt(double(NDIM))));

  Vector tmpVec;
  for (size_t iDim = 0; iDim < NDIM; iDim++)
    tmpVec[iDim] = normal_dist(Sim->ranGenerator);

  return tmpVec;
}
} // namespace dynamo
