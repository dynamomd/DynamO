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

#include <dynamo/simulation.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/outputplugins/tickerproperty/ticker.hpp>
#include <dynamo/BC/include.hpp>
#include <dynamo/systems/sysTicker.hpp>
#include <dynamo/locals/local.hpp>
#include <dynamo/species/species.hpp>
#include <dynamo/topology/topology.hpp>
#include <dynamo/globals/global.hpp>
#include <dynamo/interactions/interaction.hpp>
#include <dynamo/outputplugins/misc.hpp>
#include <dynamo/globals/PBCSentinel.hpp>
#include <boost/filesystem.hpp>
#include <dynamo/BC/BC.hpp>
#include <iomanip>
#include <set>

//! The configuration file version, a version mismatch prevents an XML file load.
static const std::string configFileVersion("1.5.0");

namespace dynamo
{
  Simulation::Simulation():
    Base("Simulation"),
    systemTime(0.0),
    eventCount(0),
    endEventCount(100000),
    eventPrintInterval(50000),
    nextPrintEvent(0),
    _force_unwrapped(false),
    primaryCellSize({1,1,1}),
    ranGenerator(std::random_device()()),
    lastRunMFT(0.0),
    simID(0),
    stateID(0),
    replexExchangeNumber(0),
    status(START)
  {}

  namespace {
    /*! \brief Hidden functor used for sorting containers of
        shared_ptr's holiding OutputPlugin classes.
     */
    struct OutputPluginSort
    {
      bool operator()(const shared_ptr<OutputPlugin>& lhs, const shared_ptr<OutputPlugin>& rhs)
      {
	      return (*lhs) < (*rhs);
      }
    };
  }

  void
  Simulation::reset()
  {
    if (status != INITIALISED)
      M_throw() << "Cannot reinitialise an un-initialised simulation";
    status = START;
    outputPlugins.clear();
    dynamics->updateAllParticles();
    systemTime = 0.0;
    eventCount = 0;
    nextPrintEvent = 0;
    lastRunMFT = 0.0;
  }

  void
  Simulation::initialise()
  {
    if (status != START)
      M_throw() << "Sim initialised at wrong time";
    
    for (shared_ptr<Species>& ptr : species)
      ptr->initialise();

    dout << "Validating Species definitions" << std::endl;
    unsigned int count = 0;
    //Now confirm that every species has only one species type!
    for (const Particle& part : particles)
      {
	for (shared_ptr<Species>& ptr : species)
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
      for (shared_ptr<Species>& ptr : species)
	tot += ptr->getCount();
    
      if (tot < N())
	M_throw() << "The particle count according to the species definition is too low\n"
		  << "discrepancy = " << tot - N()
		  << "\nN = " << N();
    
      if (tot > N())
	M_throw() << "The particle count according to the species definition is too high\n"
		  << "discrepancy = " << tot - N()
		  << "\nN = " << N();
    }

    status = SPECIES_INIT;

    dout << "Validating self-Interaction definitions" << std::endl;
    //Check that each particle has a representative interaction
    for (const Particle& particle : particles) {
      try {
	getInteraction(particle, particle);
      } catch (...) {
	M_throw() << "The particle (ID=" << particle.getID()
		  << ") does not have a self Interaction defined. Self Interactions are not used for the dynamics of the system, but are used to draw/visualise the particles, as well as calculate the excluded volume and other properties. Please add a self-Interaction";
      }
    }

    dout << "Initialising the Dynamics" << std::endl;
    dynamics->initialise();

    dout << "DOF = " << dynamics->getParticleDOF() << std::endl;
    
    status = DYNAMICS_INIT;

    dout << "Initialising Scheduler Neighbourlist" << std::endl;
    ptrScheduler->initialiseNBlist();
    
    dout << "Initialising Interactions" << std::endl;
    {
      size_t ID=0;
      
      for (shared_ptr<Interaction>& ptr : interactions)
	ptr->initialise(ID++);
    }
    
    if (std::dynamic_pointer_cast<BCPeriodic>(BCs))
      {
	double max_interaction_dist = getLongestInteraction();
	//Check that each simulation length is greater than 2x the
	//maximum interaction distance, otherwise particles can
	//interact with two periodic images!
	for (size_t i(0); i < NDIM; ++i)
	  if (primaryCellSize[i] <= (2.0 * max_interaction_dist))
	    M_throw() << "When using periodic boundary conditions, the size of the "
	      "primary image must be at least 2x the maximum interaction "
	      "distance in all dimensions, otherwise one particle can interact "
	      "with multiple periodic images of another particle."
		      << "\nprimaryCellSize[" << i << "] = "  << primaryCellSize[i]
		      << "\nLongest interaction distance = "  << max_interaction_dist;
      }

    status = INTERACTION_INIT;
	
    dout << "Initialising Locals" << std::endl;
    {
      size_t ID=0;
      //Must be initialised before globals. Neighbour lists are
      //implemented as globals and must initialise where locals are and their ID.
      for (shared_ptr<Local>& ptr : locals)
	ptr->initialise(ID++);
    }

    status = LOCAL_INIT;

    dout << "Initialising Globals" << std::endl;
    /* Add the Periodic Boundary Condition sentinel (if required). */
    if (std::dynamic_pointer_cast<BCPeriodic>(BCs))
      globals.push_back(shared_ptr<Global>(new GPBCSentinel(this, "PBCSentinel")));

    {
      size_t ID=0;
      for (shared_ptr<Global>& ptr : globals)
	ptr->initialise(ID++);
    }

    status = GLOBAL_INIT;

    dout << "Initialising Systems" << std::endl;
    //Search to check if a ticker System is needed
    for (shared_ptr<OutputPlugin>& Ptr : outputPlugins)
      if (std::dynamic_pointer_cast<OPTicker>(Ptr))
	{
	  addSystemTicker();
	  break;
	}

    {
      size_t ID=0;  
      for (shared_ptr<System>& ptr : systems)
	ptr->initialise(ID++);
    }

    status = SYSTEM_INIT;

    ensemble->initialise();

    status = ENSEMBLE_INIT;

    if (ptrScheduler == NULL)
      M_throw() << "The scheduler has not been set!";      

    dout << "Initialising Scheduler" << std::endl;
    if (endEventCount) 
      //Only initialise the scheduler if we're simulating
      ptrScheduler->initialise();

    status = SCHEDULER_INIT;

    dout << "Initialising OutputPlugins" << std::endl;
    //This sorting must be done according to the derived plugins sort
    //operators.
    std::sort(outputPlugins.begin(), outputPlugins.end(), OutputPluginSort());

    for (shared_ptr<OutputPlugin> & Ptr : outputPlugins)
      Ptr->initialise();

    status = OUTPUTPLUGIN_INIT;

    _nextPrint = eventCount + eventPrintInterval;
    status = INITIALISED;
  }

  Event 
  Simulation::getEvent(const Particle& p1, const Particle& p2) const
  {
    for (const shared_ptr<Interaction>& ptr : interactions)
      if (ptr->isInteraction(p1,p2))
	return ptr->getEvent(p1,p2);
    
    M_throw() << "Could not find the right interaction to test for";
  }

  void 
  Simulation::stream(const double dt)
  {
    BCs->update(dt);

    dynamics->stream(dt);

    for (shared_ptr<System>& ptr : systems)
      ptr->stream(dt);
  }

  double 
  Simulation::getLongestInteraction() const
  {
    double maxval = 0.0;

    for (const shared_ptr<Interaction>& ptr : interactions)
      if (ptr->maxIntDist() > maxval)
	maxval = ptr->maxIntDist();

    //Should the locals be included here?
    
    return maxval;
  }

  const shared_ptr<Interaction>&
  Simulation::getInteraction(const Particle& p1, const Particle& p2) const 
  {
    for (const shared_ptr<Interaction>& ptr : interactions)
      if (ptr->isInteraction(p1,p2))
	return ptr;
  
    M_throw() << "Could not find an Interaction between particles " << p1.getID() << " and " << p2.getID() << ". All particle pairings must have a corresponding Interaction defined.";
  }

  const shared_ptr<Species>& 
  Simulation::SpeciesContainer::operator[](const Particle& p1) const 
  {
    for (const shared_ptr<Species>& ptr : *this)
      if (ptr->isSpecies(p1)) return ptr;
    
    M_throw() << "Could not find the species corresponding to particle ID=" 
	      << p1.getID(); 
  }

  shared_ptr<Species>& 
  Simulation::SpeciesContainer::operator[](const Particle& p1) 
  {
    for (shared_ptr<Species>& ptr : *this)
      if (ptr->isSpecies(p1)) return ptr;
    
    M_throw() << "Could not find the species corresponding to particle ID=" 
	      << p1.getID(); 
  }
  
  void Simulation::addSpecies(shared_ptr<Species> sp)
  {
    if (status >= INITIALISED)
      M_throw() << "Cannot add species after simulation initialisation";
    
    species.push_back(sp);
  }

  void checkNodeNameAttribute(magnet::xml::Node node)
  {
    //Check for unique names
    std::set<std::string> names;
    for (; node.valid(); ++node)
      {
	std::string currentname = node.getAttribute("Name").as<std::string>();
	if (names.find(currentname) != names.end())
	  M_throw() << node.getName() << " at path :" << node.getPath() << "\n Does not have a unique name (Name=\"" << currentname << "\")";
	names.insert(currentname);
      }
  }

  void
  Simulation::loadXMLfile(std::string fileName)
  {
    if (status != START)
      M_throw() << "Loading config at wrong time, status = " << status;

    using namespace magnet::xml;
    
    dout << "Reading the XML input file, " << fileName << ", into memory" << std::endl;
    if (!boost::filesystem::exists(fileName))
      M_throw() << "Could not find the XML file named " << fileName
		<< "\nPlease check the file exists.";
    dout << "Parsing the XML" << std::endl;

    Document doc(fileName);

    dout << "Loading tags from the XML" << std::endl;

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

    _properties << mainNode;

    //Load the Primary cell's size
    primaryCellSize << simNode.getNode("SimulationSize");
    primaryCellSize /= units.unitLength();

    if (simNode.hasNode("Topology"))
      {
	checkNodeNameAttribute(simNode.getNode("Topology").findNode("Structure"));
	size_t i(0);
	for (magnet::xml::Node node = simNode.getNode("Topology").findNode("Structure"); node.valid(); ++node, ++i)
	  topology.push_back(Topology::getClass(node, this, i));
      }

    {
      checkNodeNameAttribute(simNode.getNode("Genus").findNode("Species"));
      size_t i(0);
      for (magnet::xml::Node node = simNode.getNode("Genus").findNode("Species"); node.valid(); ++node, ++i)
	species.push_back(Species::getClass(node, this, i));
    }
    
    BCs = BoundaryCondition::getClass(simNode.getNode("BC"), this);
    dynamics = Dynamics::getClass(simNode.getNode("Dynamics"), this);
    dynamics->loadParticleXMLData(mainNode);
    
    checkNodeNameAttribute(simNode.getNode("Interactions").findNode("Interaction"));
    for (magnet::xml::Node node = simNode.getNode("Interactions").findNode("Interaction"); node.valid(); ++node)
      interactions.push_back(Interaction::getClass(node, this));

    if (simNode.hasNode("Locals"))
      {
	checkNodeNameAttribute(simNode.getNode("Locals").findNode("Local"));
	for (magnet::xml::Node node = simNode.getNode("Locals").findNode("Local"); node.valid(); ++node)
	  locals.push_back(Local::getClass(node, this));
      }

    if (simNode.hasNode("Globals"))
      {
	checkNodeNameAttribute(simNode.getNode("Globals").findNode("Global"));
	for (magnet::xml::Node node = simNode.getNode("Globals").findNode("Global"); node.valid(); ++node)
	  globals.push_back(Global::getClass(node, this));
      }

    if (simNode.hasNode("SystemEvents"))
      {
	checkNodeNameAttribute(simNode.getNode("SystemEvents").findNode("System"));
	for (magnet::xml::Node node = simNode.getNode("SystemEvents").findNode("System"); node.valid(); ++node)
	  systems.push_back(System::getClass(node, this));
      }

    ptrScheduler = Scheduler::getClass(simNode.getNode("Scheduler"), this);
  
    //Fixes or conversions once system is loaded
    lastRunMFT *= units.unitTime();

    //Scale the loaded properties to the simulation units
    _properties.rescaleUnit(Property::Units::L, units.unitLength());
    _properties.rescaleUnit(Property::Units::T, units.unitTime());
    _properties.rescaleUnit(Property::Units::M, units.unitMass());

    ensemble = dynamo::Ensemble::loadEnsemble(*this);
  }

  void
  Simulation::writeXMLfile(std::string fileName, bool applyBC, bool round)
  {
    //Facilitate forced unwrapping when needed
    applyBC = applyBC && !_force_unwrapped;
    
    namespace xml = magnet::xml;
    xml::XmlStream XML;
    XML.setFormatXML(true);

    dynamics->updateAllParticles();

    //Rescale the properties to the configuration file units
    _properties.rescaleUnit(Property::Units::L, 1.0 / units.unitLength());
    _properties.rescaleUnit(Property::Units::T, 1.0 / units.unitTime());
    _properties.rescaleUnit(Property::Units::M, 1.0 / units.unitMass());
    
    XML << std::setprecision(std::numeric_limits<double>::digits10 + 2 - 4 * round)
	<< xml::prolog()
	<< xml::tag("DynamOconfig")
	<< xml::attr("version") << configFileVersion
	<< xml::tag("Simulation");
    
    //Allow this block to fail if need be
    if (getOutputPlugin<OPMisc>())
      {
	double mft = getOutputPlugin<OPMisc>()->getMFT();
	if (!std::isinf(mft) && !std::isnan(mft))
	  XML << xml::attr("lastMFT") << mft;
	else
	  XML << xml::attr("lastMFT") << lastRunMFT;
      }

    XML << xml::tag("Scheduler")
	<< ptrScheduler
	<< xml::endtag("Scheduler")
	<< xml::tag("SimulationSize")
	<< primaryCellSize / units.unitLength()
	<< xml::endtag("SimulationSize")
      	<< xml::tag("Genus");
  
    for (const shared_ptr<Species>& ptr : species)
      XML << xml::tag("Species") 
	  << ptr
	  << xml::endtag("Species");
  
    XML << xml::endtag("Genus")
      	<< xml::tag("BC")
	<< BCs
	<< xml::endtag("BC")
	<< xml::tag("Topology");
  
    for (const shared_ptr<Topology>& ptr : topology)
      XML << xml::tag("Structure")
	  << ptr
	  << xml::endtag("Structure");
    
    XML << xml::endtag("Topology")
	<< xml::tag("Interactions");
  
    for (const shared_ptr<Interaction>& ptr : interactions)
      XML << xml::tag("Interaction")
	  << ptr
	  << xml::endtag("Interaction");
  
    XML << xml::endtag("Interactions")
	<< xml::tag("Locals");
    
    for (const shared_ptr<Local>& ptr : locals)
      XML << xml::tag("Local")
	  << ptr
	  << xml::endtag("Local");
    
    XML << xml::endtag("Locals")
      	<< xml::tag("Globals");
    
    for (const shared_ptr<Global>& ptr : globals) XML << ptr;
    
    XML << xml::endtag("Globals")
	<< xml::tag("SystemEvents");
    
    for (const shared_ptr<System>& ptr : systems) XML << ptr;
  
    XML << xml::endtag("SystemEvents")
      	<< xml::tag("Dynamics")
	<< dynamics
	<< xml::endtag("Dynamics")
	<< xml::endtag("Simulation")
	<< _properties;

    dynamics->outputParticleXMLData(XML, applyBC);

    XML << xml::endtag("DynamOconfig");

    dout << "Config written to " << fileName << std::endl;

    //Rescale the properties back to the simulation units
    _properties.rescaleUnit(Property::Units::L, units.unitLength());
    _properties.rescaleUnit(Property::Units::T, units.unitTime());
    _properties.rescaleUnit(Property::Units::M, units.unitMass());

    XML.write_file(fileName);
  }
  
  void 
  Simulation::replexerSwap(Simulation& other)
  {
    //Get all particles up to date and zero the pecTimes
    dynamics->updateAllParticles();
    other.dynamics->updateAllParticles();
      
    std::swap(systemTime, other.systemTime);
    std::swap(eventCount, other.eventCount);
    std::swap(stateID, other.stateID);
    
    for (size_t i(0); i < systems.size(); ++i)
      systems[i]->replicaExchange(*other.systems[i]);
    
    dynamics->replicaExchange(*other.dynamics);
    
    //Rescale the velocities
    double scale1(sqrt(other.ensemble->getEnsembleVals()[2] / ensemble->getEnsembleVals()[2]));
    for (Particle& part : particles)
      part.getVelocity() *= scale1;
    //This assumes that scaling the velocities just changes the time
    //unit of the simulation. This is not true for systems with external forces!
    other.ptrScheduler->rescaleTimes(scale1);
    
    double scale2(1.0 / scale1);
    for (Particle& part : other.particles)
      part.getVelocity() *= scale2;
    ptrScheduler->rescaleTimes(scale2);

    ptrScheduler->rebuildSystemEvents();
    other.ptrScheduler->rebuildSystemEvents();    

    //Globals?
#ifdef DYNAMO_DEBUG
    if (outputPlugins.size() != other.outputPlugins.size())
      std::cerr << "Error, could not swap output plugin lists as they are not equal in size";
#endif

    for (size_t i(0); i < outputPlugins.size(); ++i)
      {
#ifdef DYNAMO_DEBUG
	if (typeid(*outputPlugins[i]) != typeid(*other.outputPlugins[i]))
	  M_throw() << "Output plugin mismatch while replexing! lists not sorted the same perhaps?";
#endif
	outputPlugins[i]->replicaExchange(*other.outputPlugins[i]);
	outputPlugins[i]->temperatureRescale(scale1 * scale1);
	other.outputPlugins[i]->temperatureRescale(scale2 * scale2);
      }
    
    //This is swapped last as things need it for calcs
    ensemble->swap(*other.ensemble);

//#ifdef DYNAMO_DEBUG
//    //Here we check that all plugins etc have the correct simulation pointer.
//    for (auto Sim : {this, &other}) {
//      for (const auto& plugin: Sim->outputPlugins)
//	if (plugin->getSimPointer() != Sim)
//	  M_throw() << "Programming error: Mismatch of Sim pointer in output plugins after replica exchange. Please file a bug report.";
//
//      for (const auto& plugin: Sim->interactions)
//	if (plugin->getSimPointer() != Sim)
//	  M_throw() << "Programming error: Mismatch of Sim pointer in Interactions after replica exchange. Please file a bug report.";
//
//      for (const auto& plugin: Sim->locals)
//	if (plugin->getSimPointer() != Sim)
//	  M_throw() << "Programming error: Mismatch of Sim pointer in Local after replica exchange. Please file a bug report.";
//
//      for (const auto& plugin: Sim->globals)
//	if (plugin->getSimPointer() != Sim)
//	  M_throw() << "Programming error: Mismatch of Sim pointer in Global after replica exchange. Please file a bug report.";
//
//      for (const auto& plugin: Sim->systems)
//	if (plugin->getSimPointer() != Sim)
//	  M_throw() << "Programming error: Mismatch of Sim pointer in System after replica exchange. Please file a bug report.";
//
//      for (const auto& plugin: Sim->topology)
//	if (plugin->getSimPointer() != Sim)
//	  M_throw() << "Programming error: Mismatch of Sim pointer in Topology after replica exchange. Please file a bug report.";
//      
//      if (Sim->ptrScheduler->getSimPointer() != Sim)
//	M_throw() << "Programming error: Mismatch of Sim pointer in Scheduler after replica exchange. Please file a bug report.";
//
//      if (Sim->dynamics->getSimPointer() != Sim)
//	M_throw() << "Programming error: Mismatch of Sim pointer in Dynamics after replica exchange. Please file a bug report.";
//
//      if (Sim->BCs->getSimPointer() != Sim)
//	M_throw() << "Programming error: Mismatch of Sim pointer in BCs after replica exchange. Please file a bug report.";
//    }
//#endif
  }

  double
  Simulation::calcInternalEnergy() const
  {
    double intECurrent = 0.0;

    for (const shared_ptr<Interaction> & plugptr : interactions)
      intECurrent += plugptr->getInternalEnergy();

    return intECurrent;
  }

  void 
  Simulation::setCOMVelocity(const Vector COMVelocity)
  {  
    magnet::math::NVector<long double, 3> sumMV({0,0,0});
    long double sumMass(0);

    //Determine the momentum discrepancy vector
    for (Particle & Part : particles)
      {
	const double mass = species[Part]->getMass(Part.getID());
	if (std::isinf(mass)) continue;
	Vector pos(Part.getPosition()), vel(Part.getVelocity());
	BCs->applyBC(pos,vel);
	sumMV += vel * mass;
	sumMass += mass;
      }
  
    sumMV /= sumMass;
  
    Vector change = COMVelocity - sumMV;
    for (Particle & Part : particles)
      {
	double mass = species[Part]->getMass(Part.getID());
	if (std::isinf(mass)) continue;
	Part.getVelocity() =  Part.getVelocity() + change;
      }
  }

  void 
  Simulation::addSystemTicker()
  {
    for (shared_ptr<System>& ptr : systems)
      if (ptr->getName() == "SystemTicker")
	M_throw() << "System Ticker already exists";
  
    systems.push_back(shared_ptr<System>(new SysTicker(this, lastRunMFT, "SystemTicker")));
  }

  double
  Simulation::getSimVolume() const
  { 
    double vol = 1.0;
    for (size_t iDim = 0; iDim < NDIM; iDim++)
      vol *= primaryCellSize[iDim];
    return vol;
  }


  double
  Simulation::getNumberDensity() const
  {
    return N() / getSimVolume();
  }

  double 
  Simulation::getPackingFraction() const
  {
    double volume = 0.0;
  
    for (const Particle& particle : particles)
      volume += getInteraction(particle, particle)->getExcludedVolume(particle.getID());
  
    return volume / getSimVolume();
  }

  size_t
  Simulation::checkSystem()
  {
    dynamics->updateAllParticles();

    size_t errors = 0;
    std::vector<Particle>::const_iterator iPtr1, iPtr2;
  
    for (const shared_ptr<Interaction>& interaction_ptr : interactions)
      {
	dout << "Checking Interaction \"" << interaction_ptr->getName() << "\"" << std::endl;
	errors += interaction_ptr->validateState();
      }

    dout << "Testing all particle pairs for invalid states" << std::endl;
    for (iPtr1 = particles.begin(); iPtr1 != particles.end(); ++iPtr1)
      for (iPtr2 = iPtr1 + 1; iPtr2 != particles.end(); ++iPtr2)
	errors += getInteraction(*iPtr1, *iPtr2)->validateState(*iPtr1, *iPtr2);

    for (const Particle& part : particles)
      for (const shared_ptr<Local>& lcl : locals)
	if (lcl->isInteraction(part))
	  errors += lcl->validateState(part);
    
    return errors;
  }

  void
  Simulation::outputData(std::string filename)
  {
    if (status < INITIALISED)
      M_throw() << "Cannot output data when not initialised!";

    namespace xml = magnet::xml;
    xml::XmlStream XML;
    XML.setFormatXML(true);
    
    XML << std::setprecision(std::numeric_limits<double>::digits10 + 2)
	<< xml::prolog() << xml::tag("OutputData");
  
    //Output the data and delete the outputplugins
    for (shared_ptr<OutputPlugin> & Ptr : outputPlugins)
      Ptr->output(XML);
  
    for (shared_ptr<Interaction> & Ptr : interactions)
      Ptr->outputData(XML);

    for (shared_ptr<Local> & Ptr : locals)
      Ptr->outputData(XML);

    for (shared_ptr<System> & Ptr : systems)
      Ptr->outputData(XML);

    XML << xml::endtag("OutputData");

    dout << "Output written to " << filename << std::endl;

    XML.write_file(filename);
  }

  void 
  Simulation::setTickerPeriod(double nP)
  {
    shared_ptr<SysTicker> ptr = std::dynamic_pointer_cast<SysTicker>(systems["SystemTicker"]);
    if (!ptr)
      M_throw() << "Could not find system ticker (maybe not required?)";

    ptr->setTickerPeriod(nP * units.unitTime());
  }

  void 
  Simulation::scaleTickerPeriod(double nP)
  {
    shared_ptr<SysTicker> ptr = std::dynamic_pointer_cast<SysTicker>(systems["SystemTicker"]);
    if (!ptr)
      M_throw() << "Could not find system ticker (maybe not required?)";

    ptr->setTickerPeriod(nP * ptr->getPeriod());
  }

  void 
  Simulation::addOutputPlugin(std::string Name)
  {
    if (status >= INITIALISED)
      M_throw() << "Cannot add plugins now";
  
    dout << "Loading output plugin string " << Name << std::endl;

    shared_ptr<OutputPlugin> tempPlug(OutputPlugin::getPlugin(Name, this));
    outputPlugins.push_back(tempPlug);
  }

  void 
  Simulation::simShutdown()
  { nextPrintEvent = endEventCount = eventCount; }

  bool
  Simulation::runSimulationStep(bool silentMode)
  {
    if (status < INITIALISED)
      M_throw() << "Bad state for runSimulation()";

    try
      {
	ptrScheduler->runNextEvent();
	
	//Periodic work
	if ((eventCount >= _nextPrint) && !silentMode && outputPlugins.size())
	  {
	    //Print the screen data plugins
	    for (shared_ptr<OutputPlugin> & Ptr : outputPlugins)
	      Ptr->periodicOutput();
	    
	    _nextPrint = eventCount + eventPrintInterval;
	    std::cout << std::endl;
	  }
      }
    catch (std::exception &cep)
      {
	M_throw() << "Exception caught while executing event " 
		  << eventCount << "\n" << cep.what();
      }

    return (eventCount < endEventCount);
  }
}
