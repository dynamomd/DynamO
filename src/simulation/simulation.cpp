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
  Base_Class("Simulation",IC_green)
{}

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
	for (; lNColl < lMaxNColl;)
	  ptrScheduler->runNextEvent();
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
	  for (; lNColl < lPrintLimiter; )
	    ptrScheduler->runNextEvent();
	  
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
      for (int i = 0; i < xMainNode.nChildNode("Plugin"); ++i)
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
