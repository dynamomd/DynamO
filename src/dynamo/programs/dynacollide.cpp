/*  dynamo:- Event driven molecular dynamics simulator
    http://www.dynamomd.org
    Copyright (C) 2013 Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
/*! \file dynarun.cpp

  \brief Contains the main() function for dynarun

  Although this contains the main() function, most of the behaviour peculiar
  to dynarun is carried out by the Coordinator class.
*/

#include <dynamo/BC/include.hpp>
#include <dynamo/dynamics/newtonian.hpp>
#include <dynamo/inputplugins/cells/include.hpp>
#include <dynamo/inputplugins/compression.hpp>
#include <dynamo/inputplugins/include.hpp>
#include <dynamo/interactions/hardsphere.hpp>
#include <dynamo/outputplugins/misc.hpp>
#include <dynamo/outputplugins/msd.hpp>
#include <dynamo/ranges/include.hpp>
#include <dynamo/schedulers/include.hpp>
#include <dynamo/schedulers/sorters/include.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/species/point.hpp>
#include <magnet/arg_share.hpp>
#include <magnet/stream/formattedostream.hpp>
#include <random>

/*! \brief Starting point for the dynarun program.

  This merely boots the Coordinator class. The Coordinator class is
  the true "main" function for the simulations.

  \param argc The number of command line arguments.
  \param argv A pointer to the array of command line arguments.
*/
int main(int argc, char *argv[]) {
  // Output the program licence
  std::cout << "dynacollide Copyright (C) 2014 Marcus N Campbell Bannerman\n"
            << "This program comes with ABSOLUTELY NO WARRANTY.\n"
            << "This is free software, and you are welcome to redistribute it\n"
            << "under certain conditions. See the licence you obtained with\n"
            << "the code\n";

  // Run the simulation
  try {
    magnet::ArgShare::getInstance().setArgs(argc, argv);
    dynamo::Simulation Sim;

    // Strip off all arguments after our special -GLGTK marker
    int args = argc;
    for (int i(0); i < argc; ++i)
      if (!strcmp("-GLGTK", argv[i])) {
        args = i;
        break;
      }

    std::mt19937 RNG;
    RNG.seed(std::random_device()());
    Sim.ranGenerator.seed(std::random_device()());
    Sim.dynamics =
        dynamo::shared_ptr<dynamo::Dynamics>(new dynamo::DynNewtonian(&Sim));
    Sim.BCs =
        dynamo::shared_ptr<dynamo::BoundaryCondition>(new dynamo::BCNone(&Sim));
    Sim.ptrScheduler =
        dynamo::shared_ptr<dynamo::SNeighbourList>(new dynamo::SNeighbourList(
            &Sim, new dynamo::FELBoundedPQ<dynamo::PELMinMax<3>>()));
    Sim.primaryCellSize = dynamo::Vector(10, 10, 10);

    const double D = 1.0, e = 1.0, M = 1.0;
    Sim.interactions.push_back(
        dynamo::shared_ptr<dynamo::Interaction>(new dynamo::IHardSphere(
            &Sim, D, e, new dynamo::IDPairRangeAll(), "Bulk")));
    Sim.addSpecies(dynamo::shared_ptr<dynamo::Species>(
        new dynamo::SpPoint(&Sim, new dynamo::IDRangeAll(&Sim), M, "Bulk", 0)));

    Sim.particles.reserve(2);
    Sim.particles.push_back(
        dynamo::Particle(dynamo::Vector(0, 0, 0), dynamo::Vector(0, 0, 0), 0));
    Sim.particles.push_back(
        dynamo::Particle(dynamo::Vector(4, 0, 0), dynamo::Vector(-1, 0, 0), 1));
    Sim.ensemble = dynamo::Ensemble::loadEnsemble(Sim);

    Sim.endEventCount = 1;
    Sim.initialise();
    while (Sim.runSimulationStep()) {
    }

    return 0;
  } catch (std::exception &cep) {
    std::cout.flush();
    magnet::stream::FormattedOStream os(
        std::cerr, magnet::console::bold() + magnet::console::red_fg() +
                       "Main(): " + magnet::console::reset());
    os << cep.what() << std::endl;
#ifndef DYNAMO_DEBUG
    os << "For a stack trace please run the debugging executables."
       << std::endl;
#endif
    return 1;
  }
}
