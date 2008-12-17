/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "replexer.hpp"
#include "../../dynamics/systems/tHalt.hpp"
#include "../../dynamics/systems/ghost.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"
#include "../../schedulers/scheduler.hpp"
#include "../../outputplugins/1partproperty/uenergy.hpp"
#include "../../extcode/threadpool.hpp"
#include <fstream>
#include <boost/random/uniform_int.hpp>
#include <limits>

void
CEReplexer::getOptions(boost::program_options::options_description& opts)
{
  boost::program_options::options_description 
    ropts("REplica EXchange Engine Options");

  ropts.add_options()
    ("sim-end-time,f", boost::program_options::value<Iflt>()->default_value(std::numeric_limits<Iflt>::max()), 
     "Simulation end time")
    ("replex-interval,i", boost::program_options::value<Iflt>()->default_value(1.0), 
     "Interval between attempting swaps")
    ("replex-swap-mode", boost::program_options::value<unsigned int>()->default_value(4), 
     "System Swap Mode:\n"
     " Values:\n"
     "  0: \tDisable swapping (For debugging or 1 system)\n"
     "  1: \tAlternating sets of pairs (~Nsims/2 attempts per swap event)\n"
     "  2: \tRandom pair per swap\n"
     "  3: \t5 * Nsim random pairs per swap\n"
     "  4: \tRandom selection of the most effective methods")
    ;
  
  opts.add(ropts);
}

CEReplexer::CEReplexer(const boost::program_options::variables_map& nVm,
		       CThreadPool& tp):
  CEngine(nVm, "config.%ID.end.xml.bz2", "output.%ID.xml.bz2", tp),
  replicaEndTime(0),
  ReplexMode(RandomSelection),
  replexSwapCalls(0),
  round_trips(0),
  SeqSelect(false),
  nSims(0),
  peekMode(false)
{}

void
CEReplexer::initialisation()
{
  preSimInit();

  for (unsigned int i = 0; i < nSims; i++)
    {
      setupSim(Simulations[i], 
	       vm["config-file"].as<std::vector<std::string> >()[i]);

      Simulations[i].initialise();

      postSimInit(Simulations[i]);

      Simulations[i].initPlugins();

      if (vm.count("ticker-period"))
	Simulations[i].setTickerPeriod(vm["ticker-period"].as<Iflt>());

      if (vm.count("scale-ticker"))
	Simulations[i].scaleTickerPeriod(vm["scale-ticker"].as<Iflt>());

    }

  //Ensure we are in the right ensemble for all simulations
  try {
    for (size_t i = nSims; i != 0;)
      dynamic_cast<const DYNAMO::CENVT&>(*(Simulations[--i].getEnsemble()));
  }
  catch (std::bad_cast&)
    { D_throw() << "One of the systems does not have an NVT ensemble"; }
  
  //Test a thermostat is available
  for (unsigned int i = 0; i < nSims; i++)
    if (Simulations[i].getSystem("Thermostat") == NULL)
      D_throw() << "Could not find the Thermostat for system " << i 
		<< "\nFilename " << vm["config-file"].as<std::vector<std::string> >()[i];
  
  //Set up the replex organisation
  temperatureList.clear();
  
  for (unsigned int i = 0; i < nSims; i++)
    {
      bool didWork = false;
      BOOST_FOREACH(smrtPlugPtr<CSystem>& sysPtr1, Simulations[i].Dynamics.getSystemEvents())
	if (sysPtr1->getName() == "Thermostat")
	  {
	    if (dynamic_cast<CSysGhost*>(sysPtr1.get_ptr()) == NULL)
	      {
		std::cout << "Could not upcast thermostat to Andersens";
		exit(1);
	      }
	    
	    temperatureList.push_back
	      (replexPair
	       (Simulations[i].getEnsemble()->getEnsembleVals()[2], 
		simData(i,Simulations[i].getEnsemble()->getReducedEnsembleVals()[2])));
	    
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
}

void 
CEReplexer::outputData()
{
  {
    std::fstream replexof("replex.dat",std::ios::out | std::ios::trunc);
    
    BOOST_FOREACH(const replexPair& myPair, temperatureList)
      replexof << myPair.second.realTemperature << " " 
	       << myPair.second.swaps << " " 
	       << (static_cast<Iflt>(myPair.second.swaps) 
		   / static_cast<Iflt>(myPair.second.attempts))  << " "
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
    ((DYNAMO::searchReplace(outputFormat, "%ID", boost::lexical_cast<std::string>(i++))).c_str());
}

void
CEReplexer::preSimInit()
{
  CEngine::preSimInit();

  ReplexMode = static_cast<Replex_Mode_Type>(vm["replex-swap-mode"].as<unsigned int>());
  
  nSims = vm["config-file"].as<std::vector<std::string> >().size();
  
  replicaEndTime = vm["sim-end-time"].as<Iflt>();
  
  if (nSims < 2 && vm.count("replex"))
    {
      std::cout << "\nTurning off replica exchange as you have Nsystems < 2";
      ReplexMode = NoSwapping;
    }
  
  if (configFormat.find("%ID") == configFormat.npos)
    D_throw() << "Replex mode, but format string for config file output"
      " doesnt contain %ID";
  
  if (outputFormat.find("%ID") == outputFormat.npos)
    D_throw() << "Multiple configs loaded, but format string for output"
      " file doesnt contain %ID";  
  
  Simulations.reset(new CSimulation[nSims]);

  //We set this straight away
  for (size_t id(0); id < nSims; ++id)
    Simulations[id].simID = id;
}


void 
CEReplexer::postSimInit(CSimulation& Sim)
{
  CEngine::postSimInit(Sim);
  Sim.addOutputPlugin("UEnergy");
}

void 
CEReplexer::forceShutdown()
{
  replicaEndTime = 0.0;
  for (unsigned int i = 0; i < nSims; i++)
    Simulations[i].simShutdown();
}

void 
CEReplexer::peekData()
{
  peekMode = true;
  
  for (unsigned int i = 0; i < nSims; i++)
    Simulations[i].simShutdown();
}

void 
CEReplexer::setupSim(CSimulation & Sim, const std::string filename)
{
  CEngine::setupSim(Sim, filename);

  Sim.addSystem(new CStHalt(&Sim, vm["replex-interval"].as<Iflt>(), 
			    "ReplexHalt"));
}

void 
CEReplexer::printStatus()
{ 
  std::cout << "Replica Exchange, ReplexSwap No." << replexSwapCalls 
	    << ", Round Trips " << round_trips
	    << "\n        T   ID     NColl   A-Ratio     Swaps    UpSims     DownSims\n";
  
  BOOST_FOREACH(const replexPair& dat, temperatureList)
    {       
      std::cout << std::setw(9)
		<< Simulations[dat.second.simID].getEnsemble()->getReducedEnsembleVals()[2] 
		<< " " << std::setw(4)
		<< dat.second.simID;
      
      std::cout << " " << std::setw(8)
		<< Simulations[dat.second.simID].getnColl()/1000 << "k" 
		<< " " << std::setw(9)
		<< ( static_cast<float>(dat.second.swaps) / dat.second.attempts)
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
CEReplexer::ReplexSwap(Replex_Mode_Type localMode)
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
	    size_t ID = boost::variate_generator<DYNAMO::baseRNG&,
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
	boost::variate_generator<DYNAMO::baseRNG&,
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
	boost::variate_generator<DYNAMO::baseRNG&,
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
	    D_throw() << "Error, randomly picked a replex move that doesn't exist";
	  }
      }
      break;
    }

}

void 
CEReplexer::ReplexSwapTicker()
{
  //Now update the histogramming
  ++replexSwapCalls;

  BOOST_FOREACH(replexPair& dat, temperatureList)
    {
      if (SimDirection[dat.second.simID])
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
CEReplexer::AttemptSwap(const unsigned int sim1ID, const unsigned int sim2ID)
{
  //D_throw() << "I broke this when this was added as a commit, check the gitk loh\nIt's something to do with improperly handled system events on exchange which I broke during the system event revamp.\n The last working commit was bf52ad0782fcf46e7690e961cc9fa51af57fb08e";
  CSimulation& sim1 = Simulations[temperatureList[sim1ID].second.simID];
  CSimulation& sim2 = Simulations[temperatureList[sim2ID].second.simID];

  temperatureList[sim1ID].second.attempts++;
  temperatureList[sim2ID].second.attempts++;
    
  //No need to check sign, it will just accept the move anyway due to
  //the [0,1) limits of the random number generator
  if (exp(sim1.getEnsemble()->exchangeProbability(*sim2.getEnsemble()))
      > boost::uniform_01<DYNAMO::baseRNG, Iflt>(sim1.ranGenerator)())
    {
      //Get all particles up to date and zero the pecTimes
      sim1.Dynamics.Liouvillean().updateAllParticles();
      sim2.Dynamics.Liouvillean().updateAllParticles();
      
      std::swap(sim1.dSysTime, sim2.dSysTime);
      std::swap(sim1.lNColl, sim2.lNColl);
      
      sim1.Dynamics.getSystemEvents().swap(sim2.Dynamics.getSystemEvents());

      BOOST_FOREACH(smrtPlugPtr<CSystem>& aPtr, sim1.Dynamics.getSystemEvents())
	aPtr->changeSystem(&sim1);

      BOOST_FOREACH(smrtPlugPtr<CSystem>& aPtr, sim2.Dynamics.getSystemEvents())
	aPtr->changeSystem(&sim2);
      
      //Rescale the velocities     
      Iflt scale1(sqrt(sim2.getEnsemble()->getEnsembleVals()[2]
		       / sim1.getEnsemble()->getEnsembleVals()[2]));

      BOOST_FOREACH(CParticle& part, sim1.vParticleList)
	part.scaleVelocity(scale1);

      sim2.ptrScheduler->rescaleTimes(scale1);
      
      Iflt scale2(1.0 / scale1);
      
      BOOST_FOREACH(CParticle& part, sim2.vParticleList)
	part.scaleVelocity(scale2);
      
      sim1.ptrScheduler->rescaleTimes(scale2);
      
      sim1.ptrScheduler->rebuildSystemEvents();
      sim2.ptrScheduler->rebuildSystemEvents();
      //If this is broken use this
      //sim1.ptrScheduler->rebuildList();
    
      //Globals?
#ifdef DYNAMO_DEBUG
      if (sim1.outputPlugins.size() != sim2.outputPlugins.size())
	std::cerr << "Error, could not swap output plugin lists as they are not equal in size";
#endif

      sim1.outputPlugins.swap(sim2.outputPlugins);      

      {
	std::vector<smrtPlugPtr<COutputPlugin> >::iterator iPtr1 = sim1.outputPlugins.begin(), 
	  iPtr2 = sim2.outputPlugins.begin();

	while (iPtr1 != sim1.outputPlugins.end())
	  {
#ifdef DYNAMO_DEBUG
	    if (typeid(*(*iPtr1)) != typeid(*(*iPtr2)))
	      D_throw() << "Output plugin mismatch while replexing! lists not sorted the same perhaps?";
#endif

	    (*iPtr1)->changeSystem(iPtr2->get_ptr());

	    (*iPtr1)->temperatureRescale(scale1);
	    (*iPtr2)->temperatureRescale(scale2);

	    ++iPtr1; 
	    ++iPtr2;
	  }
      }

      sim1.getEnsemble()->exchange(*sim2.getEnsemble());

      //Swapping the sort data
      std::swap(temperatureList[sim1ID].second.simID, temperatureList[sim2ID].second.simID);

      ++(temperatureList[sim1ID].second.swaps);
      ++(temperatureList[sim2ID].second.swaps);
    }
}

void CEReplexer::runSimulation()
{
  start_Time = boost::posix_time::second_clock::local_time();

  while ((Simulations[0].getSysTime() < replicaEndTime) && (Simulations[0].getnColl() < vm["ncoll"].as<unsigned long long>()))
    {
      if (peekMode)
	{
	  end_Time = boost::posix_time::second_clock::local_time();
		  
	  size_t i = 0;
	  BOOST_FOREACH(replexPair p1, temperatureList)
	    {
	      Simulations[p1.second.simID].setTrajectoryLength(vm["ncoll"].as<unsigned long long>());
	      Simulations[p1.second.simID].outputData((DYNAMO::searchReplace(std::string("peek.data.%ID.xml.bz2"), 
									   "%ID", boost::lexical_cast<std::string>(i++))).c_str());	    
	    }
		  
	  peekMode = false;
		  
	  {
	    std::fstream replexof("replex.dat",std::ios::out | std::ios::trunc);
		    
	    BOOST_FOREACH(const replexPair& myPair, temperatureList)
	      replexof << myPair.second.realTemperature << " " 
		       << myPair.second.swaps << " " 
		       << (static_cast<Iflt>(myPair.second.swaps) 
			   / static_cast<Iflt>(myPair.second.attempts))  << " "
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
	  //Run the simulations
	  //This is reversed as the high temperature sims generally run longer
	  for (unsigned int i = nSims; i != 0;)
	    threads.invoke(CThreadPool::task_noarg<CSimulation>
			   (Simulations[--i], &CSimulation::runSilentSimulation));
		  
	  threads.wait();//This syncs the systems for the replica exchange
		  
	  //Swap calculation
	  ReplexSwap(ReplexMode);
		  
	  ReplexSwapTicker();
		  
	  //Reset the stop events
	  for (size_t i = nSims; i != 0;)
	    {
	      //Reset the stop event
	      CStHalt* tmpRef = dynamic_cast<CStHalt*>
		(Simulations[--i].getSystem("ReplexHalt"));
		      
#ifdef DYNAMO_DEBUG
	      if (tmpRef == NULL)
		D_throw() << "Could not find the time halt event error";
#endif			
	      tmpRef->increasedt(vm["replex-interval"].as<Iflt>());

	      Simulations[i].ptrScheduler->rebuildSystemEvents();

	      //Reset the max collisions
	      Simulations[i].setTrajectoryLength(vm["ncoll"].as<unsigned long long>());
	    }
	}
    }
  end_Time = boost::posix_time::second_clock::local_time();
}

void 
CEReplexer::outputConfigs()
{
  std::fstream TtoID("TtoID.dat",std::ios::out | std::ios::trunc);
  
  int i = 0;
  BOOST_FOREACH(replexPair p1, temperatureList)
    {
      TtoID << p1.second.realTemperature << " " << i << "\n";
      Simulations[p1.second.simID].setTrajectoryLength(vm["ncoll"].as<unsigned long long>());
      Simulations[p1.second.simID].writeXMLfile((DYNAMO::searchReplace(configFormat, "%ID", boost::lexical_cast<std::string>(i++))).c_str());
    }
}
