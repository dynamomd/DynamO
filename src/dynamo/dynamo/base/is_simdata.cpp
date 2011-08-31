/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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

#include <dynamo/base/is_simdata.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/dynamics/liouvillean/liouvillean.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/dynamics/systems/system.hpp>
#include <dynamo/outputplugins/0partproperty/misc.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/chain.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/copy.hpp>
#include <iomanip>

//! The configuration file version, a version mismatch prevents an XML file load.
const char configFileVersion[] = "1.4.0";

namespace dynamo
{
  SimData::SimData():
    Base("Simulation"),
    ensemble(NULL),
    dSysTime(0.0),
    freestreamAcc(0.0),
    eventCount(0),
    endEventCount(100000),
    eventPrintInterval(50000),
    nextPrintEvent(0),
    N(0),
    ptrScheduler(NULL),
    dynamics(this),
    primaryCellSize(1,1,1),
    ranGenerator(static_cast<unsigned>(std::time(0))),
    normal_sampler(ranGenerator, boost::normal_distribution<double>()),
    uniform_sampler(ranGenerator),
    lastRunMFT(0.0),
    simID(0),
    replexExchangeNumber(0),
    status(START)
  {
  }

  SimData::~SimData()
  {
    if (ptrScheduler != NULL) delete ptrScheduler;
  }

  void
  SimData::loadXMLfile(std::string fileName)
  {
    if (status != START)
      M_throw() << "Loading config at wrong time, status = " << status;

    using namespace magnet::xml;
    Document doc;
    
    namespace io = boost::iostreams;
    
    if (!boost::filesystem::exists(fileName))
      M_throw() << "Could not find the XML file named " << fileName
		<< "\nPlease check the file exists.";
    { //This scopes out the file objects
      
      //We use the boost iostreams library to load the file into a
      //string which may be compressed.
      
      //We make our filtering iostream
      io::filtering_istream inputFile;
      
      //Now check if we should add a decompressor filter
      if (std::string(fileName.end()-8, fileName.end()) == ".xml.bz2")
	inputFile.push(io::bzip2_decompressor());
      else if (!(std::string(fileName.end()-4, fileName.end()) == ".xml"))
	M_throw() << "Unrecognized extension for xml file";

      //Finally, add the file as a source
      inputFile.push(io::file_source(fileName));
	  
      io::copy(inputFile, io::back_inserter(doc.getStoredXMLData()));
    }

    doc.parseData();

    Node mainNode = doc.getNode("DynamOconfig");

    {
      std::string version(mainNode.getAttribute("version"));
      if (version != configFileVersion)
	M_throw() << "This version of the config file is obsolete"
		  << "\nThe current version is " << configFileVersion
		  << "\nPlease look at the XMLFILE.VERSION file in the root directory of the dynamo source."
	  ;
    }

    Node subNode= mainNode.getNode("Simulation");
  
    //Don't fail if the MFT is not valid
    try {
      if (subNode.getNode("Trajectory").hasAttribute("lastMFT"))
	lastRunMFT = subNode.getNode("Trajectory").getAttribute("lastMFT").as<double>();
    } catch (std::exception&)
      {}

    ssHistory << subNode.getNode("History");

    ensemble.reset(dynamo::Ensemble::getClass(subNode.getNode("Ensemble"), this));

    _properties << mainNode;
    dynamics << mainNode;
    ptrScheduler 
      = Scheduler::getClass(subNode.getNode("Scheduler"), this);

    dynamics.getLiouvillean().loadParticleXMLData(mainNode);
  
    //Fixes or conversions once system is loaded
    lastRunMFT *= dynamics.units().unitTime();
    //Scale the loaded properties to the simulation units
    _properties.rescaleUnit(Property::Units::L, 
			    dynamics.units().unitLength());

    _properties.rescaleUnit(Property::Units::T, 
			    dynamics.units().unitTime());

    _properties.rescaleUnit(Property::Units::M, 
			    dynamics.units().unitMass());

    status = CONFIG_LOADED;
  }

  void
  SimData::writeXMLfile(std::string fileName, bool applyBC, bool round)
  {
    if (status < INITIALISED || status == ERROR)
      M_throw() << "Cannot write out configuration in this state";
  
    namespace io = boost::iostreams;
    io::filtering_ostream coutputFile;

    if (std::string(fileName.end()-4, fileName.end()) == ".bz2")
      coutputFile.push(io::bzip2_compressor());
  
    coutputFile.push(io::file_sink(fileName));
  
    magnet::xml::XmlStream XML(coutputFile);
    XML.setFormatXML(true);

    dynamics.getLiouvillean().updateAllParticles();

    //Rescale the properties to the configuration file units
    _properties.rescaleUnit(Property::Units::L, 
			    1.0 / dynamics.units().unitLength());

    _properties.rescaleUnit(Property::Units::T, 
			    1.0 / dynamics.units().unitTime());

    _properties.rescaleUnit(Property::Units::M, 
			    1.0 / dynamics.units().unitMass());

    XML << std::scientific
      //This has a minus one due to the digit in front of the decimal
      //An extra one is added if we're rounding
	<< std::setprecision(std::numeric_limits<double>::digits10 - 1 - round)
	<< magnet::xml::prolog() << magnet::xml::tag("DynamOconfig")
	<< magnet::xml::attr("version") << configFileVersion
	<< magnet::xml::tag("Simulation")
	<< magnet::xml::tag("Trajectory")
	<< magnet::xml::attr("Coll") << endEventCount
	<< magnet::xml::attr("nCollPrint") << eventPrintInterval;

    //Allow this block to fail if need be
    try {
      double mft = getOutputPlugin<OPMisc>()->getMFT();
      if (!std::isinf(mft))
	XML << magnet::xml::attr("lastMFT")
	    << mft;
    }
    catch (std::exception&)
      {}

    XML << magnet::xml::endtag("Trajectory")
	<< *ensemble
	<< magnet::xml::tag("Scheduler")
	<< *ptrScheduler
	<< magnet::xml::endtag("Scheduler")
	<< magnet::xml::tag("History") 
	<< magnet::xml::chardata()
	<< ssHistory.str()
	<< "\nRun for " << eventCount << " collisions"
	<< magnet::xml::endtag("History") << magnet::xml::endtag("Simulation")
	<< dynamics
	<< _properties;

    dynamics.getLiouvillean().outputParticleXMLData(XML, applyBC);

    XML << magnet::xml::endtag("DynamOconfig");

    dout << "Config written to " << fileName << std::endl;

    //Rescale the properties back to the simulation units
    _properties.rescaleUnit(Property::Units::L, 
			    dynamics.units().unitLength());

    _properties.rescaleUnit(Property::Units::T, 
			    dynamics.units().unitTime());

    _properties.rescaleUnit(Property::Units::M, 
			    dynamics.units().unitMass());
  }
  
  void 
  SimData::signalParticleUpdate
  (const NEventData& pdat) const
  {
    BOOST_FOREACH(const particleUpdateFunc& func, _particleUpdateNotify)
      func(pdat);
  }

  void 
  SimData::replexerSwap(SimData& other)
  {
    //Get all particles up to date and zero the pecTimes
    dynamics.getLiouvillean().updateAllParticles();
    other.dynamics.getLiouvillean().updateAllParticles();
      
    std::swap(dSysTime, other.dSysTime);
    std::swap(eventCount, other.eventCount);
    std::swap(_particleUpdateNotify, other._particleUpdateNotify);
    
    dynamics.getSystemEvents().swap(other.dynamics.getSystemEvents());

    BOOST_FOREACH(std::tr1::shared_ptr<System>& aPtr, dynamics.getSystemEvents())
      aPtr->changeSystem(this);

    BOOST_FOREACH(std::tr1::shared_ptr<System>& aPtr, other.dynamics.getSystemEvents())
      aPtr->changeSystem(&other);

    dynamics.getLiouvillean().swapSystem(other.dynamics.getLiouvillean());

    //Rescale the velocities 
    double scale1(sqrt(other.ensemble->getEnsembleVals()[2]
		       / ensemble->getEnsembleVals()[2]));
    
    BOOST_FOREACH(Particle& part, particleList)
      part.getVelocity() *= scale1;
    
    other.ptrScheduler->rescaleTimes(scale1);
    
    double scale2(1.0 / scale1);

    BOOST_FOREACH(Particle& part, other.particleList)
      part.getVelocity() *= scale2;
    
    ptrScheduler->rescaleTimes(scale2);
    
    ptrScheduler->rebuildSystemEvents();
    other.ptrScheduler->rebuildSystemEvents();    

    //Globals?
#ifdef DYNAMO_DEBUG
    if (outputPlugins.size() != other.outputPlugins.size())
      std::cerr << "Error, could not swap output plugin lists as they are not equal in size";
#endif

    outputPlugins.swap(other.outputPlugins);      
    
    {
      std::vector<magnet::ClonePtr<OutputPlugin> >::iterator iPtr1 = outputPlugins.begin(), 
	iPtr2 = other.outputPlugins.begin();
      
      while (iPtr1 != outputPlugins.end())
	{
#ifdef DYNAMO_DEBUG
	  if (typeid(*(*iPtr1)) != typeid(*(*iPtr2)))
	    M_throw() << "Output plugin mismatch while replexing! lists not sorted the same perhaps?";
#endif
	  
	  (*iPtr1)->changeSystem(iPtr2->get_ptr());
	  
	  (*iPtr1)->temperatureRescale(scale1 * scale1);
	  (*iPtr2)->temperatureRescale(scale2 * scale2);
	  
	  ++iPtr1; 
	  ++iPtr2;
	}
    }

    //This is swapped last as things need it for calcs
    ensemble->swap(*other.ensemble);
  }
}
