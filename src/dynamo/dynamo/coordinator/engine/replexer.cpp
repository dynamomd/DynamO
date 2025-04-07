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
#include <dynamo/coordinator/engine/replexer.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/systems/andersenThermostat.hpp>
#include <dynamo/systems/snapshot.hpp>
#include <dynamo/systems/tHalt.hpp>
#include <fstream>
#include <limits>
#include <magnet/string/searchreplace.hpp>
#include <magnet/thread/threadpool.hpp>

namespace dynamo {
void EReplicaExchangeSimulation::getOptions(
    boost::program_options::options_description &opts) {
  boost::program_options::options_description ropts(
      "REplica EXchange Engine Options (--engine=2)");

  ropts.add_options()(
      "replex-interval,i",
      boost::program_options::value<double>()->default_value(1.0),
      "Interval between attempting swaps on the coldest temperature. Every"
      "other systems exchange interval is scaled by (T_cold/T_i)^{1/2} to try"
      "to keep the simulation calculation times approximately"
      "constant. Otherwise the high temperature system would consume all the"
      "calculation time.")(
      "replex-swap-mode",
      boost::program_options::value<unsigned int>()->default_value(1),
      "System Swap Mode:\n"
      " Values:\n"
      "  0: \tDisable swapping (For debugging or 1 system)\n"
      "  1: \tAlternating sets of pairs (~Nsims/2 attempts per swap event)\n"
      "  2: \tRandom pair per swap\n"
      "  3: \t5 * Nsim random pairs per swap\n"
      "  4: \tRandom selection of the above methods");

  opts.add(ropts);
}

EReplicaExchangeSimulation::EReplicaExchangeSimulation(
    const boost::program_options::variables_map &nVm,
    magnet::thread::ThreadPool &tp)
    : Engine(nVm, "config.%ID.end.xml", "output.%ID.xml", tp),
      replicaEndTime(0), ReplexMode(RandomSelection), replexSwapCalls(0),
      round_trips(0), SeqSelect(false), nSims(0) {
  if (vm["events"].as<size_t>() != std::numeric_limits<size_t>::max())
    M_throw() << "You cannot use collisions to control a replica exchange "
                 "simulation\n"
              << "See the following DynamO issue: "
                 "https://github.com/toastedcrumpets/DynamO/issues/7\n";
}

void EReplicaExchangeSimulation::initialisation() {
  preSimInit();
  for (unsigned int i = 0; i < nSims; i++) {
    setupSim(Simulations[i],
             vm["config-file"].as<std::vector<std::string>>()[i]);

    if (vm.count("snapshot"))
      Simulations[i].systems.push_back(shared_ptr<System>(new SysSnapshot(
          &(Simulations[i]), vm["snapshot"].as<double>(), "SnapshotTimer",
          "ID%ID.%COUNT", !vm.count("unwrapped"))));

    if (vm.count("snapshot-events"))
      Simulations[i].systems.push_back(shared_ptr<System>(new SysSnapshot(
          &(Simulations[i]), vm["snapshot-events"].as<size_t>(),
          "SnapshotEventTimer", "%COUNTe", !vm.count("unwrapped"))));

    Simulations[i].initialise();

    postSimInit(Simulations[i]);
  }

  // Ensure we are in the right ensemble for all simulations
  for (size_t i = nSims; i != 0;)
    if (dynamic_cast<const dynamo::EnsembleNVT *>(
            Simulations[--i].ensemble.get()) == NULL)
      M_throw() << vm["config-file"].as<std::vector<std::string>>()[i]
                << " does not have an NVT ensemble";

  // Ensure the types of the simulation Dynamics match
  const Dynamics &dyn0 = *Simulations[0].dynamics;
  const std::type_info &dyntype0 = typeid(dyn0);

  for (size_t i(1); i < nSims; ++i) {
    const Dynamics &dyni = *Simulations[i].dynamics;
    if (typeid(dyni) != dyntype0)
      M_throw() << vm["config-file"].as<std::vector<std::string>>()[i]
                << " does not have the same Dynamics type as "
                << vm["config-file"].as<std::vector<std::string>>()[0];
  }

  // Test a thermostat is available
  for (size_t i = 0; i < nSims; i++)
    try {
      Simulations[i].systems["Thermostat"];
    } catch (...) {
      M_throw() << "Could not find the Thermostat for system " << i
                << "\nFilename "
                << vm["config-file"].as<std::vector<std::string>>()[i];
    }

  // Set up the replex organisation
  temperatureList.clear();

  for (unsigned int i = 1; i < nSims; i++)
    if (Simulations[0].N() != Simulations[i].N())
      M_throw() << "Every replica configuration file must have the same number "
                   "of particles!";

  for (unsigned int i = 0; i < nSims; i++) {
    shared_ptr<System> &sysPtr1 = Simulations[i].systems["Thermostat"];
    if (dynamic_cast<SysAndersen *>(sysPtr1.get()) == NULL)
      M_throw() << "Found a System event called \"Thermostat\" but could not "
                   "convert it to an Andersen Thermostat";

    temperatureList.push_back(replexPair(
        Simulations[i].ensemble->getEnsembleVals()[2],
        simData(i, Simulations[i].ensemble->getReducedEnsembleVals()[2])));
  }

  std::sort(temperatureList.begin(), temperatureList.end());

  SimDirection.resize(temperatureList.size(), 0);
  roundtrip.resize(temperatureList.size(), false);

  SimDirection[temperatureList.front().second.simID] = 1; // Going up
  SimDirection[temperatureList.back().second.simID] = -1; // Going down

  for (size_t id(0); id < nSims; ++id)
    Simulations[id].stateID = id;

  // If a system ticker is set we scale the ticker time such that the
  // number of ticks in all systems is equal.
  if (vm.count("ticker-period"))
    for (size_t i = 0; i < nSims; ++i) {
      double tFactor =
          std::sqrt(temperatureList.begin()->second.realTemperature /
                    Simulations[i].ensemble->getReducedEnsembleVals()[2]);
      Simulations[i].setTickerPeriod(vm["ticker-period"].as<double>() *
                                     tFactor);
    }

  if (vm.count("snapshot"))
    for (size_t i = 0; i < nSims; ++i) {
      double tFactor =
          std::sqrt(temperatureList.begin()->second.realTemperature /
                    Simulations[i].ensemble->getReducedEnsembleVals()[2]);
      dynamic_cast<SysSnapshot &>(*Simulations[i].systems["SnapshotTimer"])
          .setTickerPeriod(vm["snapshot"].as<double>() * tFactor);
    }
}

void EReplicaExchangeSimulation::preSimInit() {
  Engine::preSimInit();

  ReplexMode =
      static_cast<Replex_Mode_Type>(vm["replex-swap-mode"].as<unsigned int>());

  nSims = vm["config-file"].as<std::vector<std::string>>().size();

  replicaEndTime = vm["sim-end-time"].as<double>();

  if (nSims < 2 && vm.count("replex")) {
    std::cout << "\nTurning off replica exchange as you have Nsystems < 2";
    ReplexMode = NoSwapping;
  }

  if (configFormat.find("%ID") == configFormat.npos)
    M_throw() << "Replex mode, but format string for config file output"
                 " doesnt contain %ID";

  if (outputFormat.find("%ID") == outputFormat.npos)
    M_throw() << "Multiple configs loaded, but format string for output"
                 " file doesnt contain %ID";

  Simulations.reset(new Simulation[nSims]);

  // We set this straight away
  for (size_t id(0); id < nSims; ++id)
    Simulations[id].simID = id;
}

void EReplicaExchangeSimulation::setupSim(Simulation &Sim,
                                          const std::string filename) {
  Engine::setupSim(Sim, filename);
  // Add a SystHalt event to allow us to periodically halt simulations for
  // replica exchange
  Sim.systems.push_back(
      shared_ptr<System>(new SystHalt(&Sim, 0, "ReplexHalt")));
}

void EReplicaExchangeSimulation::ReplexSwap(Replex_Mode_Type localMode) {
  if (temperatureList.size() < 2)
    return;

  switch (localMode) {
  case NoSwapping:
    break;
  case SinglePair: {
    if (temperatureList.size() == 2)
      AttemptSwap(0, 1);
    else {
      // Select a image to mess with
      std::uniform_int_distribution<size_t> tmpDist(0,
                                                    temperatureList.size() - 2);
      size_t ID = tmpDist(Simulations[0].ranGenerator);
      AttemptSwap(ID, ID + 1);
    }
  } break;
  case AlternatingSequence: {
    for (size_t i = (SeqSelect) ? 0 : 1; i < (nSims - 1); i += 2)
      AttemptSwap(i, i + 1);

    SeqSelect = !SeqSelect;
  } break;
  case RandomPairs: {
    std::uniform_int_distribution<size_t> tmpDist(0,
                                                  temperatureList.size() - 1);
    size_t amount = temperatureList.size() * 5;
    for (size_t i = 0; i < amount; ++i) {
      size_t ID1(tmpDist(Simulations[0].ranGenerator)),
          ID2(tmpDist(Simulations[0].ranGenerator));

      while (ID2 == ID1)
        ID2 = tmpDist(Simulations[0].ranGenerator);

      AttemptSwap(ID1, ID2);
    }
  } break;
  case RandomSelection: {
    std::uniform_int_distribution<size_t> tmpDist(0, 1);
    switch (tmpDist(Simulations[0].ranGenerator)) {
    case 0:
      ReplexSwap(AlternatingSequence);
      break;
    case 1:
      ReplexSwap(RandomPairs);
      break;
    default:
      M_throw() << "Error, randomly picked a replex move that doesn't exist";
    }
  } break;
  }
}

void EReplicaExchangeSimulation::ReplexSwapTicker() {
  // Update the counters indicating the replexSwap count
  ++replexSwapCalls;

  for (size_t i(0); i < nSims; ++i)
    ++(Simulations[i].replexExchangeNumber);

  // Now update the histogramming
  for (replexPair &dat : temperatureList)
    if (SimDirection[dat.second.simID]) {
      if (SimDirection[dat.second.simID] > 0)
        ++dat.second.upSims;
      else
        ++dat.second.downSims;
    }

  if (SimDirection[temperatureList.front().second.simID] == -1) {
    if (roundtrip[temperatureList.front().second.simID])
      ++round_trips;

    roundtrip[temperatureList.front().second.simID] = true;
  }

  if (SimDirection[temperatureList.back().second.simID] == 1) {
    if (roundtrip[temperatureList.back().second.simID])
      ++round_trips;

    roundtrip[temperatureList.back().second.simID] = true;
  }

  SimDirection[temperatureList.front().second.simID] = 1; // Going up
  SimDirection[temperatureList.back().second.simID] = -1; // Going down
}

void EReplicaExchangeSimulation::AttemptSwap(const unsigned int sim1ID,
                                             const unsigned int sim2ID) {
  Simulation &sim1 = Simulations[temperatureList[sim1ID].second.simID];
  Simulation &sim2 = Simulations[temperatureList[sim2ID].second.simID];

  temperatureList[sim1ID].second.attempts++;
  temperatureList[sim2ID].second.attempts++;

  std::uniform_real_distribution<> uniform_dist;
  // No need to check sign, it will just accept the move anyway due to
  // the [0,1) limits of the random number generator
  if (sim1.ensemble->exchangeProbability(*sim2.ensemble) >
      uniform_dist(sim1.ranGenerator)) {
    sim1.replexerSwap(sim2);

    // Swapping the sort data
    std::swap(temperatureList[sim1ID].second.simID,
              temperatureList[sim2ID].second.simID);

    ++(temperatureList[sim1ID].second.swaps);
    ++(temperatureList[sim2ID].second.swaps);
  }
}

void EReplicaExchangeSimulation::outputData() {
  {
    std::fstream replexof("replex.dat", std::ios::out | std::ios::trunc);

    for (const replexPair &myPair : temperatureList)
      replexof << myPair.second.realTemperature << " " << myPair.second.swaps
               << " "
               << (static_cast<double>(myPair.second.swaps) /
                   static_cast<double>(myPair.second.attempts))
               << " " << myPair.second.upSims << " " << myPair.second.downSims
               << "\n";

    replexof.close();
  }

  {
    std::fstream replexof("replex.stats", std::ios::out | std::ios::trunc);

    replexof
        << "Number_of_replex_cycles " << replexSwapCalls
        << "\nTime_spent_replexing "
        << std::chrono::duration<double>(_end_time - _start_time).count() << "s"
        << "\nReplex Rate "
        << static_cast<double>(replexSwapCalls) /
               std::chrono::duration<double>(_end_time - _start_time).count()
        << "\n";

    replexof.close();
  }

  int i = 0;

  for (replexPair p1 : temperatureList)
    Simulations[p1.second.simID].outputData(
        (magnet::string::search_replace(outputFormat, "%ID",
                                        boost::lexical_cast<std::string>(i++)))
            .c_str());
}

void EReplicaExchangeSimulation::runSimulation() {
  _start_time = std::chrono::system_clock::now();

  while (((Simulations[temperatureList.front().second.simID].systemTime /
           Simulations[temperatureList.front().second.simID].units.unitTime()) <
          replicaEndTime) &&
         (Simulations[0].eventCount < vm["events"].as<size_t>())) {
    if (_SIGTERM) {
      replicaEndTime = 0.0;
      for (unsigned int i = 0; i < nSims; i++)
        Simulations[i].simShutdown();
      _SIGTERM = false;
      continue;
    }

    if (_SIGINT) {
      // Clear the writes to screen
      std::cout.flush();
      std::cerr << "\n<S>hutdown, <D>ata or <P>eek at data output:";

      char c;
      // Clear the input buffer
      std::cin.clear();
      setvbuf(stdin, NULL, _IONBF, 0);
      c = getchar();
      setvbuf(stdin, NULL, _IOLBF, 0);
      _SIGINT = false;

      switch (c) {
      case 's':
      case 'S': {
        replicaEndTime = 0.0;
        for (unsigned int i = 0; i < nSims; i++)
          Simulations[i].simShutdown();
        continue;
        break;
      }
      case 'p':
      case 'P': {
        _end_time = std::chrono::system_clock::now();

        size_t i = 0;
        for (replexPair p1 : temperatureList) {
          Simulations[p1.second.simID].endEventCount =
              vm["events"].as<size_t>();
#ifdef DYNAMO_bzip2_support
          Simulations[p1.second.simID].outputData(
              (magnet::string::search_replace(
                  std::string("peek.data.%ID.xml.bz2"), "%ID",
                  boost::lexical_cast<std::string>(i++))));
#else
          Simulations[p1.second.simID].outputData(
              (magnet::string::search_replace(
                  std::string("peek.data.%ID.xml"), "%ID",
                  boost::lexical_cast<std::string>(i++))));
#endif
        }

        {
          std::fstream replexof("replex.dat", std::ios::out | std::ios::trunc);

          for (const replexPair &myPair : temperatureList)
            replexof << myPair.second.realTemperature << " "
                     << myPair.second.swaps << " "
                     << (static_cast<double>(myPair.second.swaps) /
                         static_cast<double>(myPair.second.attempts))
                     << " " << myPair.second.upSims << " "
                     << myPair.second.downSims << "\n";

          replexof.close();
        }

        {
          std::fstream replexof("replex.stats",
                                std::ios::out | std::ios::trunc);

          replexof
              << "Number_of_replex_cycles " << replexSwapCalls
              << "\nTime_spent_replexing "
              << std::chrono::duration<double>(_end_time - _start_time).count()
              << "s"
              << "\nReplex Rate "
              << static_cast<double>(replexSwapCalls) /
                     std::chrono::duration<double>(_end_time - _start_time)
                         .count()
              << "\n";

          replexof.close();
        }
        break;
      }
      case 'd':
      case 'D': {
        std::cout << "Replica Exchange, ReplexSwap No." << replexSwapCalls
                  << ", Round Trips " << round_trips
                  << "\n        T   ID     NColl   A-Ratio     Swaps    UpSims "
                     "    DownSims\n";

        for (const replexPair &dat : temperatureList) {
          std::cout << std::setw(9)
                    << Simulations[dat.second.simID]
                           .ensemble->getReducedEnsembleVals()[2]
                    << " " << std::setw(4) << dat.second.simID << " "
                    << std::setw(8)
                    << Simulations[dat.second.simID].eventCount / 1000 << "k"
                    << " " << std::setw(9)
                    << (static_cast<double>(dat.second.swaps) /
                        dat.second.attempts)
                    << " " << std::setw(9) << dat.second.swaps << " "
                    << std::setw(9) << dat.second.upSims << " "
                    << (SimDirection[dat.second.simID] > 0 ? "/\\" : "  ")
                    << " " << std::setw(9) << dat.second.downSims << " "
                    << (SimDirection[dat.second.simID] < 0 ? "\\/" : "  ")
                    << "\n";
        }
        break;
      }
      }
      Coordinator::setup_signal_handler();
    }
    {
      // Reset the stop events
      for (size_t i = nSims; i != 0;) {
        // Reset the stop event
        shared_ptr<SystHalt> tmpRef = std::dynamic_pointer_cast<SystHalt>(
            Simulations[--i].systems["ReplexHalt"]);

#ifdef DYNAMO_DEBUG
        if (!tmpRef)
          M_throw() << "Could not find the time halt event error";
#endif
        // Each simulations exchange time is inversly proportional to its
        // temperature
        double tFactor =
            std::sqrt(temperatureList.begin()->second.realTemperature /
                      Simulations[i].ensemble->getReducedEnsembleVals()[2]);

        tmpRef->increasedt(vm["replex-interval"].as<double>() * tFactor);

        Simulations[i].scheduler->rebuildSystemEvents();

        // Reset the max collisions
        Simulations[i].endEventCount = vm["events"].as<size_t>();
      }

      // Run the simulations. We also generate all tasks at once
      // and submit them all at once to minimise lock contention.
      std::vector<std::function<void()>> tasks;
      tasks.reserve(nSims);

      for (size_t i(0); i < nSims; ++i)
        tasks.push_back(std::bind(&Simulation::runSimulation,
                                  &static_cast<Simulation &>(Simulations[i]),
                                  true));

      threads.queueTasks(tasks);
      try {
        threads.wait(); // This syncs the systems for the replica exchange
      } catch (std::exception &e) {
        int i = 0;
        std::cerr << e.what() << std::endl;
        std::cerr << "Attempting to write out configurations at the error."
                  << std::endl;
        for (replexPair p1 : temperatureList) {
          Simulations[p1.second.simID].endEventCount =
              vm["events"].as<size_t>();
          Simulations[p1.second.simID].writeXMLfile(
              magnet::string::search_replace(
                  "config.%ID.error.xml", "%ID",
                  boost::lexical_cast<std::string>(i++)),
              !vm.count("unwrapped"));
        }
        M_throw() << "Exception caught while performing simulations";
      }

      // Swap calculation
      ReplexSwap(ReplexMode);

      ReplexSwapTicker();

      double duration = std::chrono::duration<double>(
                            std::chrono::system_clock::now() - _start_time)
                            .count();

      double fractionComplete =
          (Simulations[temperatureList.front().second.simID].systemTime /
           Simulations[temperatureList.front().second.simID].units.unitTime()) /
          replicaEndTime;
      double seconds_remaining_double = duration * (1 / fractionComplete - 1);
      size_t seconds_remaining = seconds_remaining_double;

      if (seconds_remaining < std::numeric_limits<size_t>::max()) {
        size_t ETA_hours = seconds_remaining / 3600;
        size_t ETA_mins = (seconds_remaining / 60) % 60;
        size_t ETA_secs = seconds_remaining % 60;

        std::cout << "\rReplica Exchange No." << replexSwapCalls << ", ETA ";
        if (ETA_hours)
          std::cout << ETA_hours << "hr ";

        if (ETA_mins)
          std::cout << ETA_mins << "min ";

        std::cout << ETA_secs << "s        ";
        std::cout.flush();
      }
    }
  }
  _end_time = std::chrono::system_clock::now();
}

void EReplicaExchangeSimulation::outputConfigs() {
  std::fstream TtoID("TtoID.dat", std::ios::out | std::ios::trunc);

  int i = 0;
  for (replexPair p1 : temperatureList) {
    TtoID << p1.second.realTemperature << " " << i << "\n";
    Simulations[p1.second.simID].endEventCount = vm["events"].as<size_t>();
    Simulations[p1.second.simID].writeXMLfile(
        magnet::string::search_replace(configFormat, "%ID",
                                       boost::lexical_cast<std::string>(i++)),
        !vm.count("unwrapped"));
  }
}
} // namespace dynamo
