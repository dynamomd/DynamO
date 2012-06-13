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

#include <dynamo/coordinator/engine/replexer.hpp>
#include <dynamo/systems/tHalt.hpp>
#include <dynamo/systems/andersenThermostat.hpp>
#include <dynamo/liouvillean/liouvillean.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/outputplugins/1partproperty/uenergy.hpp>
#include <magnet/thread/threadpool.hpp>
#include <magnet/string/searchreplace.hpp>
#include <boost/random/uniform_int.hpp>
#include <fstream>
#include <limits>

namespace dynamo {
  void
  EReplicaExchangeSimulation::getOptions(boost::program_options::options_description& opts)
  {
    boost::program_options::options_description 
      ropts("REplica EXchange Engine Options (--engine=2)");

    ropts.add_options()
      ("replex-interval,i", boost::program_options::value<double>()->default_value(1.0), 
       "Interval between attempting swaps on the coldest temperature. Every"
       "other systems exchange interval is scaled by (T_cold/T_i)^{1/2} to try"
       "to keep the simulation calculation times approximately"
       "constant. Otherwise the high temperature system would consume all the"
       "calculation time.")
      ("replex-swap-mode", boost::program_options::value<unsigned int>()->default_value(1), 
       "System Swap Mode:\n"
       " Values:\n"
       "  0: \tDisable swapping (For debugging or 1 system)\n"
       "  1: \tAlternating sets of pairs (~Nsims/2 attempts per swap event)\n"
       "  2: \tRandom pair per swap\n"
       "  3: \t5 * Nsim random pairs per swap\n"
       "  4: \tRandom selection of the above methods")
      ;
  
    opts.add(ropts);
  }

  EReplicaExchangeSimulation::EReplicaExchangeSimulation(const boost::program_options::variables_map& nVm,
							 magnet::thread::ThreadPool& tp):
    Engine(nVm, "config.%ID.end.xml.bz2", "output.%ID.xml.bz2", tp),
    replicaEndTime(0),
    ReplexMode(RandomSelection),
    replexSwapCalls(0),
    round_trips(0),
    SeqSelect(false),
    nSims(0),
    peekMode(false)
  {
    if (vm["events"].as<unsigned long long>() != std::numeric_limits<unsigned long long>::max())
      M_throw() << "You cannot use collisions to control a replica exchange simulation\n"
		<< "See the following DynamO issue: https://github.com/toastedcrumpets/DynamO/issues/7\n";
  }


  void
  EReplicaExchangeSimulation::initialisation()
  {
    preSimInit();

    for (unsigned int i = 0; i < nSims; i++)
      {
	setupSim(Simulations[i], 
		 vm["config-file"].as<std::vector<std::string> >()[i]);

	Simulations[i].initialise();

	postSimInit(Simulations[i]);
      }

    //Ensure we are in the right ensemble for all simulations
    for (size_t i = nSims; i != 0;)
      if (dynamic_cast<const dynamo::EnsembleNVT* >(Simulations[--i].ensemble.get()) == NULL)
	M_throw() << vm["config-file"].as<std::vector<std::string> >()[i]
		  << " does not have an NVT ensemble";

    //Ensure the types of the simulation Liouvilleans match
    for (size_t i(1); i < nSims; ++i)
      if (typeid(*Simulations[i].liouvillean)
	  != typeid(*Simulations[0].liouvillean))
	M_throw() << vm["config-file"].as<std::vector<std::string> >()[i]
		  << " does not have the same Liouvillean type as "
		  << vm["config-file"].as<std::vector<std::string> >()[0];

    //Test a thermostat is available
    for (size_t i = 0; i < nSims; i++)
      try {
	Simulations[i].systems["Thermostat"];
      } catch (...) {
	M_throw() << "Could not find the Thermostat for system " << i 
		  << "\nFilename " << vm["config-file"].as<std::vector<std::string> >()[i];
      }
  
    //Set up the replex organisation
    temperatureList.clear();

    for (unsigned int i = 1; i < nSims; i++)
      if (Simulations[0].N != Simulations[i].N)
	M_throw() << "Every replica configuration file must have the same number of particles!";
  
    for (unsigned int i = 0; i < nSims; i++)
      {
	bool didWork = false;
	BOOST_FOREACH(shared_ptr<System>& sysPtr1, Simulations[i].systems)
	  if (sysPtr1->getName() == "Thermostat")
	    {
	      if (dynamic_cast<SysAndersen*>(sysPtr1.get()) == NULL)
		M_throw() << "Could not upcast thermostat to Andersens";
	    
	      temperatureList.push_back
		(replexPair
		 (Simulations[i].ensemble->getEnsembleVals()[2], 
		  simData(i,Simulations[i].ensemble->getReducedEnsembleVals()[2])));
	    
	      didWork = true;
	      break;
	    }
      
#ifdef DYNAMO_DEBUG
	if (!didWork)
	  {	 
	    std::cout << "Could not find thermostat system event";
	    exit(1);
	  }
#endif
      
      }
  
    std::sort(temperatureList.begin(), temperatureList.end());  
  
    SimDirection.resize(temperatureList.size(), 0);
    roundtrip.resize(temperatureList.size(), false);
  
    SimDirection[temperatureList.front().second.simID] = 1; //Going up
    SimDirection[temperatureList.back().second.simID] = -1; //Going down
  
    //If a system ticker is set we scale the ticker time such that the
    //number of ticks in all systems is equal.
    if (vm.count("ticker-period"))
      for (size_t i = 0; i < nSims; ++i)
	{
	  double tFactor 
	    = std::sqrt(temperatureList.begin()->second.realTemperature
			/ Simulations[i].ensemble->getReducedEnsembleVals()[2]); 
	  Simulations[i].setTickerPeriod(vm["ticker-period"].as<double>() * tFactor);
	}
  }

  void 
  EReplicaExchangeSimulation::outputData()
  {
    {
      std::fstream replexof("replex.dat",std::ios::out | std::ios::trunc);
    
      BOOST_FOREACH(const replexPair& myPair, temperatureList)
	replexof << myPair.second.realTemperature << " " 
		 << myPair.second.swaps << " " 
		 << (static_cast<double>(myPair.second.swaps) 
		     / static_cast<double>(myPair.second.attempts))  << " "
		 << myPair.second.upSims << " "
		 << myPair.second.downSims
		 << "\n";
    
      replexof.close();      
    }
  
    {      
      std::fstream replexof("replex.stats", std::ios::out | std::ios::trunc);
    
      replexof << "Number_of_replex_cycles " << replexSwapCalls
	       << "\nTime_spent_replexing " <<  boost::posix_time::to_simple_string(end_Time - start_Time)
	       << "\nReplex Rate " << static_cast<double>(replexSwapCalls) / static_cast<double>((end_Time - start_Time).total_seconds())
	       << "\n";	
    
      replexof.close();
    }    
  
    int i = 0;
  
    BOOST_FOREACH(replexPair p1, temperatureList)
      Simulations[p1.second.simID].outputData
      ((magnet::string::search_replace(outputFormat, "%ID", boost::lexical_cast<std::string>(i++))).c_str());
  }

  void
  EReplicaExchangeSimulation::preSimInit()
  {
    Engine::preSimInit();

    ReplexMode = static_cast<Replex_Mode_Type>(vm["replex-swap-mode"].as<unsigned int>());
  
    nSims = vm["config-file"].as<std::vector<std::string> >().size();
  
    replicaEndTime = vm["sim-end-time"].as<double>();
  
    if (nSims < 2 && vm.count("replex"))
      {
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

    //We set this straight away
    for (size_t id(0); id < nSims; ++id)
      Simulations[id].simID = id;
  }

  void 
  EReplicaExchangeSimulation::forceShutdown()
  {
    replicaEndTime = 0.0;
    for (unsigned int i = 0; i < nSims; i++)
      Simulations[i].simShutdown();
  }

  void 
  EReplicaExchangeSimulation::peekData()
  {
    peekMode = true;
  
    for (unsigned int i = 0; i < nSims; i++)
      Simulations[i].simShutdown();
  }

  void 
  EReplicaExchangeSimulation::setupSim(Simulation & Sim, const std::string filename)
  {
    Engine::setupSim(Sim, filename);

    //Add the halt time, set to zero so a replica exchange occurrs immediately
    Sim.systems.push_back(shared_ptr<System>(new SystHalt(&Sim, 0, "ReplexHalt")));

    Sim.addOutputPlugin("UEnergy");
  }

  void 
  EReplicaExchangeSimulation::printStatus()
  { 
    std::cout << "Replica Exchange, ReplexSwap No." << replexSwapCalls 
	      << ", Round Trips " << round_trips
	      << "\n        T   ID     NColl   A-Ratio     Swaps    UpSims     DownSims\n";

    BOOST_FOREACH(const replexPair& dat, temperatureList)
      {       
	std::cout << std::setw(9)
		  << Simulations[dat.second.simID].ensemble->getReducedEnsembleVals()[2] 
		  << " " << std::setw(4)
		  << dat.second.simID
		  << " " << std::setw(8)
		  << Simulations[dat.second.simID].eventCount/1000 << "k" 
		  << " " << std::setw(9)
		  << ( static_cast<double>(dat.second.swaps) / dat.second.attempts)
		  << " " << std::setw(9)
		  << dat.second.swaps 
		  << " " << std::setw(9)
		  << dat.second.upSims
		  << " "
		  << (SimDirection[dat.second.simID] > 0 ? "/\\" : "  ")
		  << " " << std::setw(9)
		  << dat.second.downSims
		  << " "
		  << (SimDirection[dat.second.simID] < 0 ? "\\/" : "  ")
		  << "\n";
      }
  }

  void 
  EReplicaExchangeSimulation::ReplexSwap(Replex_Mode_Type localMode)
  {
    if (temperatureList.size() < 2) return;

    switch (localMode)
      {
      case NoSwapping:
	break;
      case SinglePair:
	{
	  if (temperatureList.size() == 2)
	    AttemptSwap(0, 1);
	  else  
	    {
	      //Select a image to mess with
	      boost::uniform_int<unsigned int> tmpDist(0, temperatureList.size()-2);
	      size_t ID = boost::variate_generator<dynamo::baseRNG&,
		boost::uniform_int<unsigned int> >
		(Simulations[0].ranGenerator, tmpDist)();
	      AttemptSwap(ID, ID+1);
	    }
	}
	break;
      case AlternatingSequence:
	{
	  for (size_t i = (SeqSelect) ? 0 : 1; i < (nSims -1); i +=2)
	    AttemptSwap(i, i+1);
	
	  SeqSelect = !SeqSelect;
	}
	break;
      case RandomPairs:
	{
	  boost::variate_generator<dynamo::baseRNG&,
	    boost::uniform_int<unsigned int> >
	    rPID(Simulations[0].ranGenerator, boost::uniform_int<unsigned int>(0, temperatureList.size()-1));
	
	  size_t amount = temperatureList.size() * 5;
	  for (size_t i = 0; i < amount; ++i)
	    {
	      size_t ID1(rPID()), ID2(rPID());

	      while (ID2 == ID1)
		ID2 = rPID();

	      AttemptSwap(ID1,ID2);
	    }

	}
	break;
      case RandomSelection:
	{
	  boost::variate_generator<dynamo::baseRNG&,
	    boost::uniform_int<> >
	    rPID(Simulations[0].ranGenerator, boost::uniform_int<>(0, 1));
	
	  switch (rPID())
	    {
	    case 0:
	      ReplexSwap(AlternatingSequence);
	      break;
	    case 1:
	      ReplexSwap(RandomPairs);
	      break;
	    default:
	      M_throw() << "Error, randomly picked a replex move that doesn't exist";
	    }
	}
	break;
      }

  }

  void 
  EReplicaExchangeSimulation::ReplexSwapTicker()
  {
    //Update the counters indicating the replexSwap count
    ++replexSwapCalls;

    for (size_t i(0); i < nSims; ++i)
      ++(Simulations[i].replexExchangeNumber);

    //Now update the histogramming
    BOOST_FOREACH(replexPair& dat, temperatureList)
      if (SimDirection[dat.second.simID])
	{
	  if (SimDirection[dat.second.simID] > 0)
	    ++dat.second.upSims;
	  else
	    ++dat.second.downSims;
	}

    if (SimDirection[temperatureList.front().second.simID] == -1)
      {
	if (roundtrip[temperatureList.front().second.simID])
	  ++round_trips;
	
	roundtrip[temperatureList.front().second.simID] = true;
      }
 
    if (SimDirection[temperatureList.back().second.simID] == 1)
      {
	if (roundtrip[temperatureList.back().second.simID])
	  ++round_trips;

	roundtrip[temperatureList.back().second.simID] = true;
      }

    SimDirection[temperatureList.front().second.simID] = 1; //Going up
    SimDirection[temperatureList.back().second.simID] = -1; //Going down
  }

  void 
  EReplicaExchangeSimulation::AttemptSwap(const unsigned int sim1ID, const unsigned int sim2ID)
  {
    Simulation& sim1 = Simulations[temperatureList[sim1ID].second.simID];
    Simulation& sim2 = Simulations[temperatureList[sim2ID].second.simID];

    temperatureList[sim1ID].second.attempts++;
    temperatureList[sim2ID].second.attempts++;
    
    //No need to check sign, it will just accept the move anyway due to
    //the [0,1) limits of the random number generator
    if (sim1.ensemble->exchangeProbability(*sim2.ensemble)
	> boost::uniform_01<dynamo::baseRNG, double>(sim1.ranGenerator)())
      {
	sim1.replexerSwap(sim2);

	//Swapping the sort data
	std::swap(temperatureList[sim1ID].second.simID, temperatureList[sim2ID].second.simID);

	++(temperatureList[sim1ID].second.swaps);
	++(temperatureList[sim2ID].second.swaps);
      }
  }

  void EReplicaExchangeSimulation::runSimulation()
  {
    start_Time = boost::posix_time::second_clock::local_time();

    while ((Simulations[0].getSysTime() < replicaEndTime) && (Simulations[0].eventCount < vm["events"].as<unsigned long long>()))
      {
	if (peekMode)
	  {
	    end_Time = boost::posix_time::second_clock::local_time();
		  
	    size_t i = 0;
	    BOOST_FOREACH(replexPair p1, temperatureList)
	      {
		Simulations[p1.second.simID].endEventCount = vm["events"].as<unsigned long long>();
		Simulations[p1.second.simID].outputData((magnet::string::search_replace(std::string("peek.data.%ID.xml.bz2"), 
											"%ID", boost::lexical_cast<std::string>(i++))));
	      }
		  
	    peekMode = false;
		  
	    {
	      std::fstream replexof("replex.dat",std::ios::out | std::ios::trunc);
		    
	      BOOST_FOREACH(const replexPair& myPair, temperatureList)
		replexof << myPair.second.realTemperature << " " 
			 << myPair.second.swaps << " " 
			 << (static_cast<double>(myPair.second.swaps) 
			     / static_cast<double>(myPair.second.attempts))  << " "
			 << myPair.second.upSims << " "
			 << myPair.second.downSims
			 << "\n";
		    
	      replexof.close();      
	    }
		  
	    {      
	      std::fstream replexof("replex.stats", std::ios::out | std::ios::trunc);
		    
	      replexof << "Number_of_replex_cycles " << replexSwapCalls
		       << "\nTime_spent_replexing " <<  boost::posix_time::to_simple_string(end_Time - start_Time)
		       << "\nReplex Rate " << static_cast<double>(replexSwapCalls) / static_cast<double>((end_Time - start_Time).total_seconds())
		       << "\n";	
		    
	      replexof.close();
	    }
		  
	  }
	else 
	  {
	    //Run the simulations. We also generate all tasks at once
	    //and submit them all at once to minimise lock contention.
	    std::vector<magnet::function::Task*> tasks(nSims, NULL);
	  
	    for (size_t i(0); i < nSims; ++i)
	      tasks[i] = magnet::function::Task::makeTask(&Simulation::runSimulation, 
							  &(Simulations[i]), true);

	    threads.queueTasks(tasks);
	    threads.wait();//This syncs the systems for the replica exchange
		  
	    //Swap calculation
	    ReplexSwap(ReplexMode);
		  
	    ReplexSwapTicker();
		  
	    //Reset the stop events
	    for (size_t i = nSims; i != 0;)
	      {
		//Reset the stop event
		SystHalt* tmpRef = dynamic_cast<SystHalt*>
		  (&Simulations[--i].systems["ReplexHalt"]);
		      
#ifdef DYNAMO_DEBUG
		if (tmpRef == NULL)
		  M_throw() << "Could not find the time halt event error";
#endif			
		//Each simulations exchange time is inversly proportional to its temperature
		double tFactor 
		  = std::sqrt(temperatureList.begin()->second.realTemperature
			      / Simulations[i].ensemble->getReducedEnsembleVals()[2]); 

		tmpRef->increasedt(vm["replex-interval"].as<double>() * tFactor);

		Simulations[i].ptrScheduler->rebuildSystemEvents();

		//Reset the max collisions
		Simulations[i].endEventCount = vm["events"].as<unsigned long long>();
	      }
	  }
      }
    end_Time = boost::posix_time::second_clock::local_time();
  }

  void 
  EReplicaExchangeSimulation::outputConfigs()
  {
    std::fstream TtoID("TtoID.dat",std::ios::out | std::ios::trunc);
  
    int i = 0;
    BOOST_FOREACH(replexPair p1, temperatureList)
      {
	TtoID << p1.second.realTemperature << " " << i << "\n";
	Simulations[p1.second.simID].endEventCount = vm["events"].as<unsigned long long>();
	Simulations[p1.second.simID].writeXMLfile(magnet::string::search_replace(configFormat, "%ID", boost::lexical_cast<std::string>(i++)), 
						  !vm.count("unwrapped"));
      }
  }
}
