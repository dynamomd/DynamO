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

#include "simulation.hpp"
#include <iomanip>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/chain.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include "../dynamics/include.hpp"
#include "../schedulers/scheduler.hpp"
#include "../base/is_exception.hpp"
#include "../dynamics/include.hpp"
#include "../dynamics/interactions/intEvent.hpp"
#include "../outputplugins/outputplugin.hpp"
#include "../outputplugins/0partproperty/XMLconfig.hpp"
#include "../inputplugins/XMLconfig.hpp"
#include "../extcode/xmlParser.h"
#include "../dynamics/liouvillean/liouvillean.hpp"
#include "../dynamics/systems/system.hpp"
#include "../dynamics/NparticleEventData.hpp"
#include "../outputplugins/tickerproperty/ticker.hpp"
#include "../dynamics/systems/sysTicker.hpp"

CSimulation::CSimulation():
  Base_Class("Simulation",IC_green),
  rebuildnColl(-1),
  localeps(0)
{
  //Reasonable precision for periodic output
  std::cout << std::setprecision(std::numeric_limits<float>::digits10);
}

CSimulation::~CSimulation()
{
  if (double frac = static_cast<double>(lReverseEvents) 
      / static_cast<double>(lNColl) > 0.01 && lNColl)
    I_cerr() << "EEEEEEEK! Number of reverse events > 1 \% , =" << frac;
  
  I_cout() << lReverseEvents << " reverse events found";
}


void 
CSimulation::setTickerPeriod(Iflt nP)
{
   CSTicker* ptr = dynamic_cast<CSTicker*>(getSystem("SystemTicker"));
  if (ptr == NULL)
    D_throw() << "Could not find system ticker (maybe not required?)";

  ptr->setTickerPeriod(nP*Dynamics.units().unitTime());
}

CSystem* 
CSimulation::getSystem(std::string name)
{
  BOOST_FOREACH(smrtPlugPtr<CSystem>& sysPtr, Dynamics.getSystemEvents())
    if (sysPtr->getName() == name)
      return sysPtr.get_ptr();
  
  return NULL;
}

void 
CSimulation::addGlobal(CGlobal* tmp)
{
  if (tmp == NULL)
    D_throw() << "Adding a NULL global";

  if (status != CONFIG_LOADED)
    D_throw() << "Cannot add global events now its initialised";

  Dynamics.addGlobal(tmp);
}

void 
CSimulation::addSystem(CSystem* tmp)
{
  if (tmp == NULL)
    D_throw() << "Adding a NULL systemEvent";

  if (status != CONFIG_LOADED)
    D_throw() << "Cannot add system events now it is initialised";

  Dynamics.addSystem(tmp);
}

void 
CSimulation::addOutputPlugin(std::string Name)
{
  if (status != INITIALISED)
    D_throw() << "Cannot add plugins now";
  
  smrtPlugPtr<COutputPlugin> tempPlug(COutputPlugin::getPlugin(Name, this));
  outputPlugins.push_back(tempPlug);
}

void 
CSimulation::setRandSeed(unsigned int x)
{ ranGenerator.seed(x); }

void 
CSimulation::setnPrint(unsigned long long newnPrint)
{ 
  I_cout() << "Periodic output length set to " << newnPrint << " collisions";
  lNPrint = newnPrint; 
}

void 
CSimulation::simShutdown()
{ lPrintLimiter = lMaxNColl = lNColl; }

void 
CSimulation::setTrajectoryLength(unsigned long long newMaxColl)
{ 
  //I_cout() << "Trajectory length set to " << newMaxColl << " collisions";
  lMaxNColl = newMaxColl; 
}

void
CSimulation::initialise()
{
  if (status != CONFIG_LOADED)
    D_throw() << "Sim initialised at wrong time";

  localeps = -eps * Dynamics.units().unitTime();

  lN = vParticleList.size();
  
  I_cout() << "Initialising Simulation";  

  if (ptrScheduler == NULL)
    D_throw() << "The scheduler has not been set!";      
  
  I_cout() << "Initialising the Dynamics";
  Dynamics.initialise();

  Ensemble->initialise();
    
  I_cout() << "Initialising the scheduler";
  fflush(stdout);

  if (lMaxNColl) //Only initialise the scheduler if we're simulating
    ptrScheduler->initialise();
  
  status = INITIALISED;
}

void
CSimulation::runSimulation(bool silentMode)
{
  if (status != INITIALISED && status != PRODUCTION)
    D_throw() << "Bad state for runSimulation()";
  status = PRODUCTION;

  if (silentMode)
    try
      {
	for (; lNColl < lMaxNColl; ++lNColl)
	  executeEvent();
      }
    catch (std::exception &cep)
      {
	D_throw() << "\nWhile executing collision "
		  << lNColl << cep.what();
      }
  else
    for (lPrintLimiter = lNColl + lNPrint; lNColl < lMaxNColl; 
	 lPrintLimiter += lNPrint)
      try
	{
	  for (; lNColl < lPrintLimiter; ++lNColl)
	    executeEvent();
	  
	  //Periodic work
	  if (outputPlugins.size())
	    std::cout << "\n";
	  //Print the screen data plugins
	  BOOST_FOREACH( smrtPlugPtr<COutputPlugin> & Ptr, outputPlugins)
	    Ptr->periodicOutput();
	  
	  fflush(stdout);
	}
      catch (std::exception &cep)
	{
	  D_throw() << "\nWhile executing collision "
		    << lNColl << cep.what();
	}
}

void 
CSimulation::executeEvent()
{
#ifdef DYNAMO_DEBUG 
  if (status != PRODUCTION)
    D_throw() << "Attempted Collision execution at improper stage";
#endif

  switch (ptrScheduler->nextEventType())
    {
    case Interaction:
      executeIntEvent();
      break;
    case Global:
      executeGlobEvent();
      break;
    case System:
#ifdef DYNAMO_DEBUG
	if (Dynamics.getSystemEvents().empty()) 
		D_throw() << "A system event has been scheduled yet there are no system events";
#endif
      executeSysEvent();
      break;
    }
}

void 
CSimulation::executeIntEvent()
{
  CIntEvent iEvent = ptrScheduler->earliestIntEvent();

  if (iEvent.getType() == NONE)
    {
      I_cerr() << "A glancing or tenuous collision may have become invalid due"
	"\nto free streaming inaccuracies"
	"\nOtherwise the simulation has run out of events!"
	"\nThis occured when confirming the event with the scheduler"
	"\nIgnoring this NONE event below\n"
	       << iEvent.stringData(this);

      //Now we're past the event, update the scheduler and plugins
      ptrScheduler->update(iEvent.getParticle1(), iEvent.getParticle2());
      
      return;
    }
    
  if (iEvent.getdt() < localeps)
    ++lReverseEvents;
  
#ifdef DYNAMO_DEBUG 
  if (isnan(iEvent.getdt()))
    D_throw() << "A NAN Interaction collision time has been found"
	      << iEvent.stringData(this);
  
  if (iEvent.getdt() == HUGE_VAL)
    D_throw() << "An infinite Interaction (not marked as NONE) collision time has been found\n"
	      << iEvent.stringData(this);
#endif

  //Debug section
#ifdef DYNAMO_CollDebug
  if (iEvent.getParticle1().getID() < iEvent.getParticle2().getID())
    std::cerr << "\nsysdt " << iEvent.getdt() + dSysTime
	      << "  ID1 " << iEvent.getParticle1().getID() 
	      << "  ID2 " << iEvent.getParticle2().getID()
	      << "  dt " << iEvent.getdt()
	      << "  Type " << CIntEvent::getCollEnumName(iEvent.getType());
  else
    std::cerr << "\nsysdt " << iEvent.getdt() + dSysTime
	      << "  ID1 " << iEvent.getParticle2().getID() 
	      << "  ID2 " << iEvent.getParticle1().getID()
	      << "  dt " << iEvent.getdt()
	      << "  Type " << CIntEvent::getCollEnumName(iEvent.getType());
#endif
  
  dSysTime += iEvent.getdt();
    
  ptrScheduler->stream(iEvent.getdt());
  
  //dynamics must be updated first
  Dynamics.stream(iEvent.getdt());
  
  //Run the collision and catch the data
  C2ParticleData EDat = Dynamics.runEvent(iEvent);
  
  //Now we're past the event, update the scheduler and plugins
  ptrScheduler->update(iEvent.getParticle1(), iEvent.getParticle2());
  
  BOOST_FOREACH( smrtPlugPtr<COutputPlugin> & Ptr, outputPlugins)
    Ptr->eventUpdate(iEvent,EDat);
}

void 
CSimulation::executeGlobEvent()
{
  CGlobEvent iEvent = ptrScheduler->earliestGlobEvent();
  
  if (iEvent.getType() == NONE)
    D_throw() << "No global collision found\n"
	      << iEvent.stringData(this);
  
  if (iEvent.getdt() < localeps)
    ++lReverseEvents;

  /*if (iEvent.getdt() < localeps)
    D_throw() << "Reverse Time Global Collision!\n" 
    << iEvent.stringData(this); */
  
#ifdef DYNAMO_DEBUG 
  if (isnan(iEvent.getdt()))
    D_throw() << "A NAN Global collision time has been found\n"
	      << iEvent.stringData(this);
  
  if (iEvent.getdt() == HUGE_VAL)
    D_throw() << "An infinite (not marked as NONE) Global collision time has been found\n"
	      << iEvent.stringData(this);
#endif
  
  dSysTime += iEvent.getdt();
  
  ptrScheduler->stream(iEvent.getdt());
  
  //dynamics must be updated first
  Dynamics.stream(iEvent.getdt());
  
  //Run the collision and catch the data
  CNParticleData EDat = Dynamics.runEvent(iEvent);
  
  //Now we're past the event update the scheduler and plugins
  ptrScheduler->update(iEvent.getParticle());
  
  BOOST_FOREACH( smrtPlugPtr<COutputPlugin> & Ptr, outputPlugins)
    Ptr->eventUpdate(iEvent,EDat);
}

void 
CSimulation::executeSysEvent()
{
  CSystem& sEvent = **min_element(Dynamics.getSystemEvents().begin(), Dynamics.getSystemEvents().end());
  //Must have dt copied here as we dont have systemEvent classes
  Iflt dt = sEvent.getdt();
  
  if (dt < localeps)
    D_throw() << "Reverse Time system event!\n"; 
  
#ifdef DYNAMO_DEBUG 
  if (isnan(dt))
    D_throw() << "A NAN system event time has been found";
#endif

    
  dSysTime += dt;
    
  ptrScheduler->stream(dt);
  
  //dynamics must be updated first
  Dynamics.stream(dt);
  
  //Run the collision and catch the data
  CNParticleData SDat = sEvent.runEvent();

#ifdef DYNAMO_CollDebug
    std::cerr << "\nSYSTEM sysdt " << dt + dSysTime
	      << "  dt " << dt;
#endif
  
  //Now we're past the event update the scheduler and plugins
  BOOST_FOREACH(const C1ParticleData& pData, SDat.L1partChanges)
    {
#ifdef DYNAMO_CollDebug
      std::cerr << "\nPart Single =" << pData.getParticle().getID();
#endif      
      ptrScheduler->update(pData.getParticle());
    }

  BOOST_FOREACH(const C2ParticleData& pData, SDat.L2partChanges)
    {
#ifdef DYNAMO_CollDebug
      std::cerr << "\nPart Pair 1 =" <<pData.particle1_.getParticle().getID()
		<< "  Pair 2 = " << pData.particle2_.getParticle().getID();
#endif

      ptrScheduler->update(pData.particle1_.getParticle());
      ptrScheduler->update(pData.particle2_.getParticle());
    }
  
  BOOST_FOREACH(smrtPlugPtr<COutputPlugin>& Ptr, outputPlugins)
    {
      Ptr->eventUpdate(sEvent, SDat, dt);
    }
}

void 
CSimulation::configLoaded()
{
  //Handled by an input plugin
  if (status != START)
    D_throw() << "Loading config at wrong time";
  
  status = CONFIG_LOADED;
}

void
CSimulation::loadXMLfile(const char *fileName)
{
  //Handled by an input plugin
  if (status != START)
    D_throw() << "Loading config at wrong time";
  
  CIPConfig XMLconfig(fileName,this);
  XMLconfig.initialise();
	  
  status = CONFIG_LOADED;
}

void
CSimulation::writeXMLfile(const char *fileName)
{
  if (status < INITIALISED || status == ERROR)
    D_throw() << "Cannot write out configuration in this state";
  
  //Particle data output handled by an output plugin
  COPConfig XMLconfig(this);
  XMLconfig.fileOutput(fileName);
}

void 
CSimulation::initPlugins()
{
  I_cout() << "Sort and init the Output Plugins";
  std::sort(outputPlugins.begin(), outputPlugins.end());
  
  bool needTicker = false;
  
  BOOST_FOREACH(smrtPlugPtr<COutputPlugin> & Ptr, outputPlugins)
    {
      Ptr->initialise();
      if (!needTicker && dynamic_cast<COPTicker*>(Ptr.get_ptr()) != NULL)
	needTicker = true; 
    }

  if (needTicker)
    Dynamics.addSystemTicker();
}

void
CSimulation::loadPlugins(std::string pluginFileName)
{
  XMLNode xMainNode;

  if (!boost::filesystem::exists(pluginFileName))
    D_throw() << "Plugin file \"" << pluginFileName << "\" doesn't exist";


  if (std::string(pluginFileName.end()-4, pluginFileName.end()) == ".xml")
    {
      xMainNode=XMLNode::openFileHelper(pluginFileName.c_str(), "Plugins");
      smrtPlugPtr<COutputPlugin> tmpPlug(NULL);
      for (int i = 0; i < xMainNode.nChildNode("Plugin"); i++)
	{
	  tmpPlug.set_ptr(COutputPlugin::getPlugin(xMainNode.getChildNode("Plugin", i), this));
	  outputPlugins.push_back(tmpPlug);
	}
    }
  else
    D_throw() << "plugin filename should end in .xml and be xml";
}


void
CSimulation::outputData(const char* filename)
{
  if (status < INITIALISED || status == ERROR)
    D_throw() << "Cannot output data when not initialised!";
  
  namespace io = boost::iostreams;
  
  io::filtering_ostream coutputFile;
  coutputFile.push(io::bzip2_compressor());
  coutputFile.push(io::file_sink(filename));
  
  xmlw::XmlStream XML(coutputFile);
  
  XML << std::setprecision(std::numeric_limits<Iflt>::digits10)
      << xmlw::prolog() << xmlw::tag("OutputData");
  
  //Output the data and delete the outputplugins
  BOOST_FOREACH( smrtPlugPtr<COutputPlugin> & Ptr, outputPlugins)
    Ptr->output(XML);

  XML << xmlw::endtag("OutputData");
}

long double 
CSimulation::getSysTime()
{ return dSysTime / Dynamics.units().unitTime(); }



/* Obsolete

void 
CSimulation::setNThreads(size_t nT)
{
  I_cout() << "Thread count set to " << nT;
  threadPool.setMaxThreads(nT);
}

void 
CSimulation::streamTaskInit()
{
  //Handle a small or unset nStreamTasks
  if (!(threadPool.getMaxThreads()))
    {
      nStreamTasks = 0;
      I_cout() << "Not threaded, using a single stream task";
      streamTasks.push_back(CStreamTask(vParticleList.begin(), 
					vParticleList.end(),
					this));
      return;
    }
  else if (nStreamTasks == 0)
    nStreamTasks = vParticleList.size()/10000;
  
  if (nStreamTasks == 0)
    {
      I_cout() << "Single threaded streaming enabled due to system size";
      streamTasks.push_back(CStreamTask(vParticleList.begin(), 
					vParticleList.end(),
					this));
      return;
    }
  
  unsigned long step = vParticleList.size()/nStreamTasks;
  
  std::vector<CParticle>::iterator start = vParticleList.begin();

  for (int i = 1; i < nStreamTasks; i++)
    {
      streamTasks.push_back(CStreamTask(start,start+step,this));
      start += step;
    }
  
  I_cout() << "Using " << nStreamTasks <<" streaming thread tasks";

  streamTasks.push_back(CStreamTask(start,vParticleList.end(),this));
}


  //Wait before proceeding
  //threadPool.wait();

  I_cout() << "Generating thread streaming tasks";
  streamTaskInit();


void 
CSimulation::setnStreamTasks(unsigned int x)
{
  I_cout() << "Number of StreamTasks set to " << x;
  nStreamTasks = x; 
}
*/
