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

#include <dynamo/coordinator/coordinator.hpp>
#include <dynamo/coordinator/engine/single.hpp>
#include <dynamo/systems/snapshot.hpp>
#include <dynamo/systems/visualizer.hpp>
#include <stdio.h>

namespace dynamo {
ESingleSimulation::ESingleSimulation(
    const boost::program_options::variables_map &nVM,
    magnet::thread::ThreadPool &tp)
    : Engine(nVM, "config.out.xml", "output.xml", tp) {}

void ESingleSimulation::runSimulation() {
  try {
    while (true) {
      if (!simulation.runSimulationStep())
        break;
      if (_SIGINT) {
        // Clear the writes to screen
        std::cout.flush();
        std::cerr << "\n<S>hutdown or <P>eek at data output:";

        char c;
        // Clear the input buffer
        std::cin.clear();
        setvbuf(stdin, NULL, _IONBF, 0);
        c = getchar();
        setvbuf(stdin, NULL, _IOLBF, BUFSIZ);
        switch (c) {
        case 's':
        case 'S':
          simulation.simShutdown();
          break;
        case 'p':
        case 'P':
          simulation.outputData("peek.data.xml.bz2");
          break;
        }

        _SIGINT = false;
        Coordinator::setup_signal_handler();
      }
      if (_SIGTERM) {
        _SIGTERM = false;
        simulation.simShutdown();
      }
    }
  } catch (std::exception &cep) {
    try {
      std::cerr << "\nEngine: Trying to output config to config.error.xml.bz2, "
                   "and output to output.error.xml.bz2";
      simulation.writeXMLfile("config.error.xml.bz2", !vm.count("unwrapped"));
      simulation.outputData("output.error.xml.bz2");
    } catch (...) {
      std::cerr
          << "\nEngine: Could not write out config/output in error state!";
    }
    throw;
  }
}

void ESingleSimulation::initialisation() {
  preSimInit();

  if (!(vm.count("config-file")) ||
      (vm["config-file"].as<std::vector<std::string>>().size() != 1))
    M_throw() << "You must only provide one input file in single mode";

  setupSim(simulation, vm["config-file"].as<std::vector<std::string>>()[0]);

#ifdef DYNAMO_visualizer
  if (_loadVisualizer)
    simulation.systems.push_back(shared_ptr<System>(new SVisualizer(
        &simulation, vm["config-file"].as<std::vector<std::string>>()[0],
        simulation.lastRunMFT)));
#endif

  if (vm.count("snapshot"))
    simulation.systems.push_back(shared_ptr<System>(
        new SysSnapshot(&simulation, vm["snapshot"].as<double>(),
                        "SnapshotTimer", "%COUNT", !vm.count("unwrapped"))));

  if (vm.count("snapshot-events"))
    simulation.systems.push_back(shared_ptr<System>(new SysSnapshot(
        &simulation, vm["snapshot-events"].as<size_t>(), "SnapshotEventTimer",
        "%COUNTe", !vm.count("unwrapped"))));

  simulation.initialise();

  postSimInit(simulation);

  if (vm.count("ticker-period"))
    simulation.setTickerPeriod(vm["ticker-period"].as<double>());
}

void ESingleSimulation::outputData() {
  simulation.outputData(outputFormat.c_str());
}

void ESingleSimulation::outputConfigs() {
  simulation.writeXMLfile(configFormat.c_str(), !vm.count("unwrapped"));
}
} // namespace dynamo
