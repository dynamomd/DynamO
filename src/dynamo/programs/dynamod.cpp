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

#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <dynamo/BC/include.hpp>
#include <dynamo/inputplugins/include.hpp>
#include <dynamo/schedulers/include.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/systems/andersenThermostat.hpp>
#include <magnet/exception.hpp>

namespace po = boost::program_options;

int main(int argc, char *argv[]) {
  std::cout << "dynamod  Copyright (C) 2013  Marcus N Campbell Bannerman\n"
            << "This program comes with ABSOLUTELY NO WARRANTY.\n"
            << "This is free software, and you are welcome to redistribute it\n"
            << "under certain conditions. See the licence you obtained with\n"
            << "the code\n";

  dynamo::Simulation sim;

  ////////////////////////PROGRAM OPTIONS!!!!!!!!!!!!!!!!!!!!!!!
  try {
    po::options_description allopts("General Options"),
        loadopts("Load Config File Options"), hiddenopts, helpOpts;

#ifdef DYNAMO_bzip2_support
    std::string extension(".bz2");
#else
    std::string extension("");
#endif

    allopts.add_options()(
        "help,h", "Produces this message OR if --pack-mode/-m is set, it lists "
                  "the specific options available for that packer mode.")(
        "out-config-file,o",
        po::value<string>()->default_value("config.out.xml" + extension),
        "Configuration output file.")(
        "random-seed,s", po::value<unsigned int>(),
        "Seed value for the random number generator.")(
        "rescale-T,r", po::value<double>(),
        "Rescales the kinetic temperature of the input/generated config to "
        "this value.")(
        "thermostat,T", po::value<double>(),
        "Change or add a thermostatt with the temperature provided. A "
        "temperature of zero will remove the thermostatt.")(
        "zero-momentum,Z",
        "Zeros the total momentum of the input/generated config.")(
        "zero-com", "Zeros the centre of mass of the input/generated config.")(
        "zero-vel", po::value<size_t>(),
        "Sets the velocity in the [arg=0,1,or 2] dimension of each particle to "
        "zero.")(
        "set-com-vel", po::value<std::string>(),
        "Sets the velocity of the COM of the system (format x,y,z no spaces).")(
        "mirror-system,M", po::value<unsigned int>(),
        "Mirrors the particle co-ordinates and velocities. Argument is "
        "dimension to reverse/mirror.")(
        "round", "Output the XML config file with one less digit of accuracy "
                 "to remove rounding errors (used in the test harness).")(
        "unwrapped", "Don't apply the boundary conditions of the system when "
                     "writing out the particle positions.")(
        "check", "Runs tests on the configuration to ensure the system is not "
                 "in an invalid state.");

    loadopts.add_options()("config-file", po::value<string>(),
                           "Config file to initialise from (Non packer mode).");

    helpOpts.add(allopts);
    helpOpts.add(dynamo::IPPacker::getOptions());

    allopts.add(loadopts);
    allopts.add(dynamo::IPPacker::getOptions());

    hiddenopts.add_options()("b1", "boolean option one.")(
        "b2", "boolean option two.")("i1", po::value<size_t>(),
                                     "integer option one.")(
        "i2", po::value<size_t>(), "integer option two.")(
        "i3", po::value<size_t>(), "integer option three.")(
        "i4", po::value<size_t>(), "integer option four.")(
        "s1", po::value<std::string>(), "string option one.")(
        "s2", po::value<std::string>(),
        "string option two.")("f1", po::value<double>(), "double option one.")(
        "f2", po::value<double>(), "double option two.")(
        "f3", po::value<double>(), "double option three.")(
        "f4", po::value<double>(), "double option four.")(
        "f5", po::value<double>(),
        "double option five.")("f6", po::value<double>(), "double option six.")(
        "f7", po::value<double>(), "double option seven.")(
        "f8", po::value<double>(), "double option eight.")(
        "f9", po::value<double>(), "double option nine.")(
        "f10", po::value<double>(), "double option ten.")(
        "NCells,C", po::value<unsigned long>()->default_value(7),
        "Default number of unit cells per dimension, used for crystal packing "
        "of particles.")("xcell,x", po::value<unsigned long>(),
                         "Number of unit cells in the x dimension.")(
        "ycell,y", po::value<unsigned long>(),
        "Number of unit cells in the y dimension.")(
        "zcell,z", po::value<unsigned long>(),
        "Number of unit cells in the z dimension.")(
        "rectangular-box",
        "Force the simulation box to be deformed so "
        "that the x,y,z cells also specify the box aspect ratio.")(
        "density,d", po::value<double>()->default_value(0.5),
        "System number density.");

    allopts.add(hiddenopts);

    po::positional_options_description p;
    p.add("config-file", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv)
                  .options(allopts)
                  .positional(p)
                  .run(),
              vm);
    po::notify(vm);

    if (vm.count("random-seed"))
      sim.ranGenerator.seed(vm["random-seed"].as<unsigned int>());

    if (!vm.count("pack-mode") &&
        (vm.count("help") || !vm.count("config-file"))) {
      cout << "Usage : dynamod <OPTIONS>...[CONFIG FILE]\n"
           << " Either modifies a config file (if a file name is passed as an "
              "argument) OR generates a new config file depending on the "
              "packing mode (if --pack-mode/-m is used).\n"
           << helpOpts;
      return 1;
    }

    if (vm.count("pack-mode") && vm.count("config-file")) {
      cout << "You cannot specify a packing mode and pass a configuration file "
              "as an argument"
           << std::endl;
      return 1;
    }

    ////////////////////////Simulation Initialisation!!!!!!!!!!!!!
    // Now load the config
    if (vm.count("pack-mode")) {
      dynamo::IPPacker plug(vm, &sim);
      plug.initialise();

      // We don't zero momentum and rescale for certain packer modes
      if ((vm["pack-mode"].as<size_t>() != 23) &&
          (vm["pack-mode"].as<size_t>() != 25) &&
          (vm["pack-mode"].as<size_t>() != 28)) {
        dynamo::InputPlugin(&sim, "Rescaler").zeroMomentum();
        dynamo::InputPlugin(&sim, "Rescaler").rescaleVels(1.0);
      }
    } else
      sim.loadXMLfile(vm["config-file"].as<string>());

    sim.endEventCount = 0;

    if (vm.count("thermostat")) {
      if (vm["thermostat"].as<double>() == 0.0) {
        auto thermostat_base_ptr = sim.systems.find("Thermostat");
        if (thermostat_base_ptr == sim.systems.end())
          M_throw() << "Could not locate thermostat to disable";

        if (!std::dynamic_pointer_cast<dynamo::SysAndersen>(
                *thermostat_base_ptr))
          M_throw() << "Could not upcast System event named \"Thermostat\" to "
                       "Thermostat type";

        sim.systems.erase(thermostat_base_ptr);
        sim.ensemble = dynamo::Ensemble::loadEnsemble(sim);
      } else {
        if (sim.systems.find("Thermostat") == sim.systems.end())
          // Should use reduced temperature here, but we set it later anyway
          sim.systems.push_back(dynamo::shared_ptr<dynamo::System>(
              new dynamo::SysAndersen(&sim, 1.0 / sim.N(), 1.0, "Thermostat")));

        auto thermostat_base_ptr = sim.systems.find("Thermostat");
        auto thermostat_ptr = std::dynamic_pointer_cast<dynamo::SysAndersen>(
            *thermostat_base_ptr);

        if (!thermostat_ptr)
          M_throw() << "Could not upcast System event named \"Thermostat\" to "
                       "SysAndersen";

        thermostat_ptr->setReducedTemperature(vm["thermostat"].as<double>());
        sim.ensemble = dynamo::Ensemble::loadEnsemble(sim);
      }
    }

    sim.initialise();

    // Here we modify the sim accordingly
    if (vm.count("zero-momentum"))
      dynamo::InputPlugin(&sim, "MomentumZeroer").zeroMomentum();

    if (vm.count("check"))
      sim.checkSystem();

    if (vm.count("zero-com"))
      dynamo::InputPlugin(&sim, "CentreOfMassZeroer").zeroCentreOfMass();

    if (vm.count("rescale-T"))
      dynamo::InputPlugin(&sim, "Rescaler")
          .rescaleVels(vm["rescale-T"].as<double>());

    if (vm.count("mirror-system"))
      dynamo::InputPlugin(&sim, "Mirrorer")
          .mirrorDirection(vm["mirror-system"].as<unsigned int>());

    if (vm.count("set-com-vel")) {
      boost::tokenizer<boost::char_separator<char>> tokens(
          vm["set-com-vel"].as<std::string>(),
          boost::char_separator<char>(","));

      boost::tokenizer<boost::char_separator<char>>::iterator details_iter =
          tokens.begin();

      dynamo::Vector vel{0, 0, 0};

      if (details_iter == tokens.end())
        M_throw() << "set-com-vel requires 3 components";
      vel[0] = boost::lexical_cast<double>(*(details_iter++));
      if (details_iter == tokens.end())
        M_throw() << "set-com-vel requires 3 components";
      vel[1] = boost::lexical_cast<double>(*(details_iter++));
      if (details_iter == tokens.end())
        M_throw() << "set-com-vel requires 3 components";
      vel[2] = boost::lexical_cast<double>(*(details_iter));

      dynamo::InputPlugin(&sim, "velSetter").setCOMVelocity(vel);
    }

    if (vm.count("zero-vel"))
      dynamo::InputPlugin(&sim, "Vel-Component-Zeroer")
          .zeroVelComp(vm["zero-vel"].as<size_t>());

    sim.writeXMLfile(vm["out-config-file"].as<string>(), !vm.count("unwrapped"),
                     vm.count("round"));
  } catch (std::exception &cep) {
    std::cout.flush();
    magnet::stream::FormattedOStream os(std::cerr, "Main(): ");
    os << cep.what() << std::endl;
#ifndef DYNAMO_DEBUG
    os << "Try using the debugging executable for more information on the "
          "error."
       << std::endl;
#endif

    return 1;
  }

  return 0;
}
