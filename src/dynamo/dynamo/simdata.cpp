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

#include <dynamo/simdata.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/liouvillean/liouvillean.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/systems/sysTicker.hpp>
#include <dynamo/locals/local.hpp>
#include <dynamo/species/species.hpp>
#include <dynamo/topology/topology.hpp>
#include <dynamo/globals/global.hpp>
#include <dynamo/interactions/interaction.hpp>
#include <dynamo/outputplugins/0partproperty/misc.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/chain.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/copy.hpp>
#include <dynamo/BC/BC.hpp>
#include <iomanip>

//! The configuration file version, a version mismatch prevents an XML file load.
static const std::string configFileVersion("1.5.0");

namespace dynamo
{
  SimData::SimData():
    Base("Simulation"),
    dSysTime(0.0),
    freestreamAcc(0.0),
    eventCount(0),
    endEventCount(100000),
    eventPrintInterval(50000),
    nextPrintEvent(0),
    N(0),
    primaryCellSize(1,1,1),
    ranGenerator(static_cast<unsigned>(std::time(0))),
    normal_sampler(ranGenerator, boost::normal_distribution<double>()),
    uniform_sampler(ranGenerator, boost::uniform_01<double>()),
    lastRunMFT(0.0),
    simID(0),
    replexExchangeNumber(0),
    status(START)
  {}

  void
  SimData::initialise()
  {
    BOOST_FOREACH(shared_ptr<Species>& ptr, species)
      ptr->initialise();

    unsigned int count = 0;
    //Now confirm that every species has only one species type!
    BOOST_FOREACH(const Particle& part, particleList)
      {
	BOOST_FOREACH(shared_ptr<Species>& ptr, species)
	  if (ptr->isSpecies(part)) { count++; break; }
	
	if (count < 1)
	  M_throw() << "Particle ID=" << part.getID() << " has no species";
	
	if (count > 1)
	  M_throw() << "Particle ID=" << part.getID() << " has more than one species";
	count = 0;
      }
    
    //Now confirm that there are not more counts from each species
    //than there are particles
    {
      unsigned long tot = 0;
      BOOST_FOREACH(shared_ptr<Species>& ptr, species)
	tot += ptr->getCount();
    
      if (tot < N)
	M_throw() << "The particle count according to the species definition is too low\n"
		  << "discrepancy = " << tot - N
		  << "\nN = " << N;
    
      if (tot > N)
	M_throw() << "The particle count according to the species definition is too high\n"
		  << "discrepancy = " << tot - N
		  << "\nN = " << N;
    }

    liouvillean->initialise();

    {
      size_t ID=0;
      
      BOOST_FOREACH(shared_ptr<Interaction>& ptr, interactions)
	ptr->initialise(ID++);
    }

    {
      size_t ID=0;
      //Must be initialised before globals. Neighbour lists are
      //implemented as globals and must initialise where locals are and their ID.
      BOOST_FOREACH(shared_ptr<Local>& ptr, locals)
	ptr->initialise(ID++);
    }

    {
      size_t ID=0;
      
      BOOST_FOREACH(shared_ptr<Global>& ptr, globals)
	ptr->initialise(ID++);
    }

    {
      size_t ID=0;
      
      BOOST_FOREACH(shared_ptr<System>& ptr, systems)
	ptr->initialise(ID++);
    }
  }

  IntEvent 
  SimData::getEvent(const Particle& p1, const Particle& p2) const
  {
    BOOST_FOREACH(const shared_ptr<Interaction>& ptr, interactions)
      if (ptr->isInteraction(p1,p2))
	{
#ifdef dynamo_UpdateCollDebug
	  std::cerr << "\nGOT INTERACTION P1 = " << p1.getID() << " P2 = " 
		    << p2.getID() << " NAME = " << typeid(*(ptr.get())).name();
#endif
	  return ptr->getEvent(p1,p2);
	}
    
    M_throw() << "Could not find the right interaction to test for";
  }

  void 
  SimData::stream(const double dt)
  {
    BCs->update(dt);

    liouvillean->stream(dt);

    BOOST_FOREACH(shared_ptr<System>& ptr, systems)
      ptr->stream(dt);
  }

  double 
  SimData::getLongestInteraction() const
  {
    double maxval = 0.0;

    BOOST_FOREACH(const shared_ptr<Interaction>& ptr, interactions)
      if (ptr->maxIntDist() > maxval)
	maxval = ptr->maxIntDist();

    return maxval;
  }

  const shared_ptr<Interaction>&
  SimData::getInteraction(const Particle& p1, const Particle& p2) const 
  {
    BOOST_FOREACH(const shared_ptr<Interaction>& ptr, interactions)
      if (ptr->isInteraction(p1,p2))
	return ptr;
  
    M_throw() << "Could not find the interaction requested";
  }

  const Species& 
  SimData::SpeciesContainer::operator[](const Particle& p1) const 
  {
    BOOST_FOREACH(const shared_ptr<Species>& ptr, *this)
      if (ptr->isSpecies(p1)) return *ptr;
    
    M_throw() << "Could not find the species corresponding to particle ID=" 
	      << p1.getID(); 
  }

  void SimData::addSpecies(shared_ptr<Species> sp)
  {
    if (status >= INITIALISED)
      M_throw() << "Cannot add species after simulation initialisation";

    species.push_back(sp);

    BOOST_FOREACH(shared_ptr<Interaction>& intPtr, interactions)
      if (intPtr->isInteraction(*sp))
	{
	  sp->setIntPtr(intPtr.get());
	  return;
	}

    M_throw() << "Could not find the interaction for the species \"" 
	      << sp->getName() << "\"";
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

    Node simNode= mainNode.getNode("Simulation");
  
    //Don't fail if the MFT is not valid
    try {
      if (simNode.hasAttribute("lastMFT"))
	lastRunMFT = simNode.getAttribute("lastMFT").as<double>();
    } catch (std::exception&)
      {}

    ensemble = dynamo::Ensemble::getClass(simNode.getNode("Ensemble"), this);

    _properties << mainNode;

    //Load the Primary cell's size
    primaryCellSize << simNode.getNode("SimulationSize");
    primaryCellSize /= units.unitLength();

    { 
      size_t i(0);
      for (magnet::xml::Node node = simNode.getNode("Genus").fastGetNode("Species");
	   node.valid(); ++node, ++i)
	species.push_back(Species::getClass(node, this, i));
    }
    
    BCs = BoundaryCondition::getClass(simNode.getNode("BC"), this);

    liouvillean = Liouvillean::getClass(simNode.getNode("Dynamics"), this);

    if (simNode.hasNode("Topology"))
      {
	size_t i(0);
	for (magnet::xml::Node node = simNode.getNode("Topology").fastGetNode("Structure");
	     node.valid(); ++node, ++i)
	  topology.push_back(Topology::getClass(node, this, i));
      }

    for (magnet::xml::Node node = simNode.getNode("Interactions").fastGetNode("Interaction");
	 node.valid(); ++node)
      interactions.push_back(Interaction::getClass(node, this));

    //Link the species and interactions
    BOOST_FOREACH(shared_ptr<Species>& sp , species)
      BOOST_FOREACH(shared_ptr<Interaction>& intPtr, interactions)
      if (intPtr->isInteraction(*sp))
	{
	  sp->setIntPtr(intPtr.get());
	  break;
	}

    if (simNode.hasNode("Locals"))
      for (magnet::xml::Node node = simNode.getNode("Locals").fastGetNode("Local"); 
	   node.valid(); ++node)
	locals.push_back(Local::getClass(node, this));

    if (simNode.hasNode("Globals"))
      for (magnet::xml::Node node = simNode.getNode("Globals").fastGetNode("Global"); 
	   node.valid(); ++node)
	globals.push_back(Global::getClass(node, this));

    if (simNode.hasNode("SystemEvents"))
      for (magnet::xml::Node node = simNode.getNode("SystemEvents").fastGetNode("System"); 
	   node.valid(); ++node)
	systems.push_back(System::getClass(node, this));

    ptrScheduler = Scheduler::getClass(simNode.getNode("Scheduler"), this);

    liouvillean->loadParticleXMLData(mainNode);
  
    //Fixes or conversions once system is loaded
    lastRunMFT *= units.unitTime();
    //Scale the loaded properties to the simulation units
    _properties.rescaleUnit(Property::Units::L, 
			    units.unitLength());

    _properties.rescaleUnit(Property::Units::T, 
			    units.unitTime());

    _properties.rescaleUnit(Property::Units::M, 
			    units.unitMass());
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

    liouvillean->updateAllParticles();

    //Rescale the properties to the configuration file units
    _properties.rescaleUnit(Property::Units::L, 
			    1.0 / units.unitLength());

    _properties.rescaleUnit(Property::Units::T, 
			    1.0 / units.unitTime());

    _properties.rescaleUnit(Property::Units::M, 
			    1.0 / units.unitMass());

    XML << std::scientific
      //This has a minus one due to the digit in front of the decimal
      //An extra one is added if we're rounding
	<< std::setprecision(std::numeric_limits<double>::digits10 - 1 - round)
	<< magnet::xml::prolog() << magnet::xml::tag("DynamOconfig")
	<< magnet::xml::attr("version") << configFileVersion
	<< magnet::xml::tag("Simulation");
    
    //Allow this block to fail if need be
    if (getOutputPlugin<OPMisc>())
      {
	double mft = getOutputPlugin<OPMisc>()->getMFT();
	if (!std::isinf(mft))
	  XML << magnet::xml::attr("lastMFT")
	      << mft;
      }

    XML	<< *ensemble
	<< magnet::xml::tag("Scheduler")
	<< *ptrScheduler
	<< magnet::xml::endtag("Scheduler")
	<< magnet::xml::tag("SimulationSize")
	<< primaryCellSize / units.unitLength()
	<< magnet::xml::endtag("SimulationSize")
      	<< magnet::xml::tag("Genus");
  
    BOOST_FOREACH(const shared_ptr<Species>& ptr, species)
      XML << magnet::xml::tag("Species")
	  << *ptr
	  << magnet::xml::endtag("Species");
  
    XML << magnet::xml::endtag("Genus")
      	<< magnet::xml::tag("BC")
	<< *BCs
	<< magnet::xml::endtag("BC")
	<< magnet::xml::tag("Topology");
  
    BOOST_FOREACH(const shared_ptr<Topology>& ptr, topology)
      XML << magnet::xml::tag("Structure")
	  << *ptr
	  << magnet::xml::endtag("Structure");
    
    XML << magnet::xml::endtag("Topology")
	<< magnet::xml::tag("Interactions");
  
    BOOST_FOREACH(const shared_ptr<Interaction>& ptr, interactions)
      XML << magnet::xml::tag("Interaction")
	  << *ptr
	  << magnet::xml::endtag("Interaction");
  
    XML << magnet::xml::endtag("Interactions")
	<< magnet::xml::tag("Locals");
    
    BOOST_FOREACH(const shared_ptr<Local>& ptr, locals)
      XML << magnet::xml::tag("Local")
	  << *ptr
	  << magnet::xml::endtag("Local");
    
    XML << magnet::xml::endtag("Locals")
      	<< magnet::xml::tag("Globals");
    
    BOOST_FOREACH(const shared_ptr<Global>& ptr, globals)
      XML << *ptr;
    
    XML << magnet::xml::endtag("Globals")
	<< magnet::xml::tag("SystemEvents");
    
    BOOST_FOREACH(const shared_ptr<System>& ptr, systems)
      XML << *ptr;
  
    XML << magnet::xml::endtag("SystemEvents")
      	<< magnet::xml::tag("Dynamics")
	<< *liouvillean
	<< magnet::xml::endtag("Dynamics")
	<< magnet::xml::endtag("Simulation")
	<< _properties;

    liouvillean->outputParticleXMLData(XML, applyBC);

    XML << magnet::xml::endtag("DynamOconfig");

    dout << "Config written to " << fileName << std::endl;

    //Rescale the properties back to the simulation units
    _properties.rescaleUnit(Property::Units::L, 
			    units.unitLength());

    _properties.rescaleUnit(Property::Units::T, 
			    units.unitTime());

    _properties.rescaleUnit(Property::Units::M, 
			    units.unitMass());
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
    liouvillean->updateAllParticles();
    other.liouvillean->updateAllParticles();
      
    std::swap(dSysTime, other.dSysTime);
    std::swap(eventCount, other.eventCount);
    std::swap(_particleUpdateNotify, other._particleUpdateNotify);
    
    systems.swap(other.systems);

    BOOST_FOREACH(shared_ptr<System>& aPtr, systems)
      aPtr->changeSystem(this);

    BOOST_FOREACH(shared_ptr<System>& aPtr, other.systems)
      aPtr->changeSystem(&other);

    liouvillean->swapSystem(*other.liouvillean);

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
      std::vector<shared_ptr<OutputPlugin> >::iterator iPtr1 = outputPlugins.begin(), 
	iPtr2 = other.outputPlugins.begin();
      
      while (iPtr1 != outputPlugins.end())
	{
#ifdef DYNAMO_DEBUG
	  if (typeid(*(*iPtr1)) != typeid(*(*iPtr2)))
	    M_throw() << "Output plugin mismatch while replexing! lists not sorted the same perhaps?";
#endif
	  
	  (*iPtr1)->changeSystem(iPtr2->get());
	  
	  (*iPtr1)->temperatureRescale(scale1 * scale1);
	  (*iPtr2)->temperatureRescale(scale2 * scale2);
	  
	  ++iPtr1; 
	  ++iPtr2;
	}
    }

    //This is swapped last as things need it for calcs
    ensemble->swap(*other.ensemble);
  }

  double
  SimData::calcInternalEnergy() const
  {
    double intECurrent = 0.0;

    BOOST_FOREACH(const shared_ptr<Interaction> & plugptr, interactions)
      intECurrent += plugptr->getInternalEnergy();

    return intECurrent;
  }

  void 
  SimData::setCOMVelocity(const Vector COMVelocity)
  {  
    Vector sumMV(0,0,0);
 
    long double sumMass(0);

    //Determine the discrepancy VECTOR
    BOOST_FOREACH(Particle & Part, particleList)
      {
	Vector  pos(Part.getPosition()), vel(Part.getVelocity());
	BCs->applyBC(pos,vel);
	double mass = species[Part].getMass(Part.getID());
	//Note we sum the negatives!
	sumMV -= vel * mass;
	sumMass += mass;
      }
  
    sumMV /= sumMass;
  
    sumMV += COMVelocity;

    BOOST_FOREACH(Particle & Part, particleList)
      Part.getVelocity() =  Part.getVelocity() + sumMV;
  }

  void 
  SimData::addSystemTicker()
  {
    BOOST_FOREACH(shared_ptr<System>& ptr, systems)
      if (ptr->getName() == "SystemTicker")
	M_throw() << "System Ticker already exists";
  
    systems.push_back(shared_ptr<System>(new SysTicker(this, lastRunMFT, "SystemTicker")));
  }

  double
  SimData::getSimVolume() const
  { 
    double vol = 1.0;
    for (size_t iDim = 0; iDim < NDIM; iDim++)
      vol *= primaryCellSize[iDim];
    return vol;
  }


  double
  SimData::getNumberDensity() const
  {
    return N / getSimVolume();
  }

  double 
  SimData::getPackingFraction() const
  {
    double volume = 0.0;
  
    BOOST_FOREACH(const shared_ptr<Species>& sp, species)
      BOOST_FOREACH(const size_t& ID, *(sp->getRange()))
      volume += sp->getIntPtr()->getExcludedVolume(ID);
  
    return  volume / getSimVolume();
  }

  void
  SimData::SystemOverlapTest()
  {
    liouvillean->updateAllParticles();

    std::vector<Particle>::const_iterator iPtr1, iPtr2;
  
    for (iPtr1 = particleList.begin(); iPtr1 != particleList.end(); ++iPtr1)
      for (iPtr2 = iPtr1 + 1; iPtr2 != particleList.end(); ++iPtr2)    
	getInteraction(*iPtr1, *iPtr2)->checkOverlaps(*iPtr1, *iPtr2);

    BOOST_FOREACH(const Particle& part, particleList)
      BOOST_FOREACH(const shared_ptr<Local>& lcl, locals)
      if (lcl->isInteraction(part))
	lcl->checkOverlaps(part);
  }

}
