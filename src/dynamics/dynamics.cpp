/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "dynamics.hpp"
#include "include.hpp"
#include <boost/foreach.hpp>
#include "../datatypes/vector.hpp"
#include "../datatypes/vector.xml.hpp"
#include "../extcode/xmlwriter.hpp"
#include "../extcode/xmlParser.h"
#include <magnet/exception.hpp>
#include <cmath>
#include "../base/is_simdata.hpp"
#include "NparticleEventData.hpp"
#include "systems/sysTicker.hpp"
#include "../schedulers/scheduler.hpp"

Dynamics::Dynamics(DYNAMO::SimData* tmp): 
  SimBase(tmp,"Dynamics",IC_purple),
  p_BC(NULL), 
  p_liouvillean(NULL),
  p_units(NULL)
{}

Dynamics::Dynamics(const XMLNode& XML, DYNAMO::SimData* tmp): 
  SimBase(tmp, "Dynamics", IC_purple),
  p_BC(NULL), 
  p_units(NULL)
{ operator<<(XML); }

Dynamics::~Dynamics() {}

ClonePtr<Topology>& 
Dynamics::getTopology(std::string name)
{
  BOOST_FOREACH(ClonePtr<Topology>& sysPtr, topology)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  M_throw() << "Could not find the topology " << name;
}

const ClonePtr<Topology>& 
Dynamics::getTopology(std::string name) const
{
  BOOST_FOREACH(const ClonePtr<Topology>& sysPtr, topology)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  M_throw() << "Could not find the topology " << name;
}

const Species& 
Dynamics::getSpecies(const Particle& p1) const 
{
  BOOST_FOREACH(const ClonePtr<Species>& ptr, species)
    if (ptr->isSpecies(p1))
      return *ptr;
  
  M_throw() << "Could not find the requested species"
	    << "\nID = " << p1.getID();
}

xml::XmlStream& operator<<(xml::XmlStream& XML, 
			    const Dynamics& g)
{
  g.outputXML(XML);
  return XML;
}

const Species& 
Dynamics::getSpecies(std::string name) const
{
  BOOST_FOREACH(const ClonePtr<Species>& ptr, species)
    if (ptr->getName() == name)
      return *ptr;
  
  M_throw() << "Could not find the " << name << " species"; 
}

ClonePtr<System>&
Dynamics::getSystem(std::string name)
{
  BOOST_FOREACH(ClonePtr<System>& sysPtr, systems)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  M_throw() << "Could not find system plugin";
}

const ClonePtr<System>&
Dynamics::getSystem(std::string name) const
{
  BOOST_FOREACH(const ClonePtr<System>& sysPtr, systems)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  M_throw() << "Could not find system plugin";
}

ClonePtr<Global>&
Dynamics::getGlobal(std::string name)
{
  BOOST_FOREACH(ClonePtr<Global>& sysPtr, globals)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  M_throw() << "Could not find global plugin";
}

const ClonePtr<Global>&
Dynamics::getGlobal(std::string name) const
{
  BOOST_FOREACH(const ClonePtr<Global>& sysPtr, globals)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  M_throw() << "Could not find global plugin";
}

ClonePtr<Local>&
Dynamics::getLocal(std::string name)
{
  BOOST_FOREACH(ClonePtr<Local>& sysPtr, locals)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  M_throw() << "Could not find local plugin";
}

const ClonePtr<Local>&
Dynamics::getLocal(std::string name) const
{
  BOOST_FOREACH(const ClonePtr<Local>& sysPtr, locals)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  M_throw() << "Could not find local plugin";
}

ClonePtr<Interaction>&
Dynamics::getInteraction(std::string name)
{
  BOOST_FOREACH(ClonePtr<Interaction>& sysPtr, interactions)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  M_throw() << "Could not find interaction plugin";
}

const ClonePtr<Interaction>&
Dynamics::getInteraction(std::string name) const
{
  BOOST_FOREACH(const ClonePtr<Interaction>& sysPtr, interactions)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  M_throw() << "Could not find interaction plugin";
}

void 
Dynamics::addSpecies(const ClonePtr<Species>& CSpe)
{
  if (Sim->status >= INITIALISED)
    M_throw() << "Cannot add species after simulation initialisation";

  species.push_back(CSpe);

  BOOST_FOREACH(ClonePtr<Interaction>& intPtr , interactions)
    {
      if (intPtr->isInteraction(*species.back()))
	{
	  species.back()->setIntPtr(intPtr.get_ptr());
	  break;
	}
    }
}

void 
Dynamics::addGlobal(Global* newGlobal)
{
  if (Sim->status >= INITIALISED)
    M_throw() << "Cannot add global events after simulation initialisation";

  ClonePtr<Global> 
    tempPlug(newGlobal);
  
  globals.push_back(tempPlug);
}

void 
Dynamics::addLocal(Local* newLocal)
{
  if (Sim->status >= INITIALISED)
    M_throw() << "Cannot add local events after simulation initialisation";

  ClonePtr<Local> 
    tempPlug(newLocal);
  
  locals.push_back(tempPlug);
}

void
Dynamics::addSystem(System* newSystem)
{
  if (Sim->status >= INITIALISED)
    M_throw() << "Cannot add system events at this time, system is initialised";
  
  ClonePtr<System> 
    tempPlug(newSystem);
  
  systems.push_back(tempPlug); 
}

void
Dynamics::addStructure(Topology* newSystem)
{ 
  if (Sim->status >= INITIALISED)
    M_throw() << "Cannot add structure after simulation initialisation";

  ClonePtr<Topology> 
    tempPlug(newSystem);
  
  topology.push_back(tempPlug); 
}

void 
Dynamics::addSystemTicker()
{
  if (Sim->status >= INITIALISED)
    M_throw() << "Cannot add the system ticker now";

  BOOST_FOREACH(ClonePtr<System>& ptr, systems)
    if (ptr->getName() == "SystemTicker")
      M_throw() << "System Ticker already exists";
  
    addSystem(new CSTicker(Sim, Sim->lastRunMFT, "SystemTicker"));
}

Interaction* 
Dynamics::addInteraction(Interaction* CInt)
{
  ClonePtr<Interaction> tempPlug(CInt);
  interactions.push_back(tempPlug);
  return interactions.back().get_ptr();
}

void 
Dynamics::initialise()
{
  BOOST_FOREACH(ClonePtr<Species>& ptr, species)
    ptr->initialise();
  
  unsigned int count = 0;
  //Now confirm that every species has only one species type!
  BOOST_FOREACH(const Particle& part, Sim->particleList)
    {
      BOOST_FOREACH(ClonePtr<Species>& ptr, species)
	if (ptr->isSpecies(part)) { count++; break; }
      
      if (count < 1)
	M_throw() << "Particle ID=" << part.getID() << " has no species";

      if (count > 1)
	M_throw() << "Particle ID=" << part.getID() << " has more than one species";
      count = 0;
    }

  //Now confirm that there are not more counts from each species than there are particles
  {
    unsigned long tot = 0;
    BOOST_FOREACH(ClonePtr<Species>& ptr, species)
      tot += ptr->getCount();
    
    if (tot < Sim->N)
      M_throw() << "The particle count according to the species definition is too low\n"
		<< "discrepancy = " << tot - Sim->N
		<< "\nN = " << Sim->N;
    
    if (tot > Sim->N)
      M_throw() << "The particle count according to the species definition is too high\n"
		<< "discrepancy = " << tot - Sim->N
		<< "\nN = " << Sim->N;
 }

  p_liouvillean->initialise();

  size_t ID=0;

  BOOST_FOREACH(ClonePtr<Interaction>& ptr, interactions)
    ptr->initialise(ID++);

  ID=0;

  //Must be initialised before globals. Neighbour lists are
  //implemented as globals and must initialise where locals are and their ID.
  BOOST_FOREACH(ClonePtr<Local>& ptr, locals)
    ptr->initialise(ID++);

  ID=0;

  BOOST_FOREACH(ClonePtr<Global>& ptr, globals)
    ptr->initialise(ID++);

  ID=0;

  BOOST_FOREACH(ClonePtr<System>& ptr, systems)
    ptr->initialise(ID++);
}

const ClonePtr<Interaction>&
Dynamics::getInteraction(const Particle& p1, const Particle& p2) const 
{
  BOOST_FOREACH(const ClonePtr<Interaction>& ptr, interactions)
    if (ptr->isInteraction(p1,p2))
      return ptr;
  
  M_throw() << "Could not find the interaction requested";
}

Dynamics::Dynamics(const Dynamics &dyn):
  SimBase(dyn),
  p_BC(dyn.p_BC), 
  p_units(dyn.p_units)
{}

void 
Dynamics::stream(const Iflt& dt)
{
  p_BC->update(dt);

  p_liouvillean->stream(dt);

  BOOST_FOREACH(ClonePtr<System>& ptr, systems)
    ptr->stream(dt);
}


Iflt
Dynamics::calcInternalEnergy() const
{
  Iflt intECurrent = 0.0;

  BOOST_FOREACH(const ClonePtr<Interaction> & plugptr, 
		Sim->dynamics.getInteractions())
    intECurrent += plugptr->getInternalEnergy();

  return intECurrent;
}

Iflt 
Dynamics::getNumberDensity() const
{
  return Sim->N / Sim->dynamics.units().simVolume();
}

Iflt 
Dynamics::getPackingFraction() const
{
  Iflt volume = 0.0;
  
  BOOST_FOREACH(const ClonePtr<Species>& sp, Sim->dynamics.getSpecies())
    volume += pow(sp->getIntPtr()->hardCoreDiam(), NDIM) * sp->getCount();
  
  return  PI * volume / (6 * (Sim->dynamics.units().simVolume()));
}

void 
Dynamics::setCOMVelocity(const Vector COMVelocity)
{  
  Vector sumMV(0,0,0);
 
  lIflt sumMass(0);

  //Determine the discrepancy VECTOR
  BOOST_FOREACH(Particle & Part, Sim->particleList)
    {
      Vector  pos(Part.getPosition()), vel(Part.getVelocity());
      BCs().applyBC(pos,vel);

      //Note we sum the negatives!
      sumMV -= vel * getSpecies(Part).getMass();
      sumMass += getSpecies(Part).getMass();
    }
  
  sumMV /= sumMass;
  
  sumMV += COMVelocity;

  BOOST_FOREACH(Particle & Part, Sim->particleList)
    Part.getVelocity() =  Part.getVelocity() + sumMV;
}

void
Dynamics::operator<<(const XMLNode& XML)
{
  I_cout() << "Loading dynamics from XML";
  
  XMLNode xDynamics=XML.getChildNode("Dynamics"), xSubNode;
  
  //Load the aspect ratio
  if (xDynamics.hasChild("Aspect_Ratio"))
    {
      xSubNode = xDynamics.getChildNode("Aspect_Ratio");
      Sim->aspectRatio << xSubNode;
    }
  
  xSubNode = xDynamics.getChildNode("Units");
  p_units.set_ptr(Units::loadUnits(xSubNode,Sim));
  
  //Now load the BC part, after the aspect ratio!
  xSubNode = xDynamics.getChildNode("BC");
  p_BC.set_ptr(BoundaryCondition::loadClass(xSubNode, Sim));
  
  if (xDynamics.hasChild("Topology"))
    {
      xSubNode = xDynamics.getChildNode("Topology");  
      for (long i=0; i < xSubNode.nChildNode("Structure"); i++)
	{
	  ClonePtr<Topology> tempPlug(Topology::loadClass(xSubNode.getChildNode("Structure",i),Sim,i));
	  topology.push_back(tempPlug);
	}
    }
  
  xSubNode = xDynamics.getChildNode("Genus");  
  for (long i=0; i < xSubNode.nChildNode("Species"); i++)
    species.push_back(ClonePtr<Species>
		      (Species::getClass(xSubNode.getChildNode("Species",i),Sim,i)));
  
  xSubNode = xDynamics.getChildNode("Liouvillean");
  p_liouvillean.set_ptr(Liouvillean::loadClass(xSubNode,Sim));  
  
  xSubNode = xDynamics.getChildNode("Interactions");
  for (long i=0; i < xSubNode.nChildNode("Interaction"); i++)
    {
      ClonePtr<Interaction> tempPlug(Interaction::getClass
					 (xSubNode.getChildNode("Interaction",
								i),Sim));
      interactions.push_back(tempPlug);
    }  
  
  //Link the species and interactions
  BOOST_FOREACH(ClonePtr<Species>& sp , species)
    BOOST_FOREACH(ClonePtr<Interaction>& intPtr , interactions)
    if (intPtr->isInteraction(*sp))
      {
	sp->setIntPtr(intPtr.get_ptr());
	break;
      }
  
  if (xDynamics.hasChild("Globals"))
    {
      xSubNode = xDynamics.getChildNode("Globals");  
      for (long i = 0; i < xSubNode.nChildNode("Global"); ++i)
	{
	  ClonePtr<Global> 
	    tempPlug(Global::getClass(xSubNode.getChildNode("Global", i), Sim));
	  
	  globals.push_back(tempPlug);
	}
    }

  if (xDynamics.hasChild("Locals"))
    {
      xSubNode = xDynamics.getChildNode("Locals");
      
      for (long i = 0; i < xSubNode.nChildNode("Local"); ++i)
	{
	  ClonePtr<Local> 
	    tempPlug(Local::getClass(xSubNode.getChildNode("Local", i), Sim));
	  
	  locals.push_back(tempPlug);
	}
    }
  
  if (xDynamics.hasChild("SystemEvents"))
    {
      xSubNode = xDynamics.getChildNode("SystemEvents");
      
      for (long i=0; i < xSubNode.nChildNode("System"); i++)
	{
	  ClonePtr<System> 
	    tempPlug(System::getClass(xSubNode.getChildNode("System",i),Sim));
	  
	  systems.push_back(tempPlug);
	}
    }
}

void
Dynamics::outputXML(xml::XmlStream &XML) const
{
  XML << xml::tag("Dynamics")
      << xml::tag("Aspect_Ratio")
      << Sim->aspectRatio
      << xml::endtag("Aspect_Ratio")
      << xml::tag("Units")
      << p_units
      << xml::endtag("Units")
      << xml::tag("BC")
      << p_BC
      << xml::endtag("BC")
      << xml::tag("Genus");
  
  BOOST_FOREACH(const ClonePtr<Species>& ptr, species)
    XML << xml::tag("Species")
	<< ptr
	<< xml::endtag("Species");
  
  XML << xml::endtag("Genus")
      << xml::tag("Topology");
  
  BOOST_FOREACH(const ClonePtr<Topology>& ptr, topology)
    XML << xml::tag("Structure")
	<< ptr
	<< xml::endtag("Structure");
  
  XML << xml::endtag("Topology")
      << xml::tag("SystemEvents");
  
  BOOST_FOREACH(const ClonePtr<System>& ptr, systems)
    XML << ptr;
  
  XML << xml::endtag("SystemEvents")
      << xml::tag("Globals");
  
  BOOST_FOREACH(const ClonePtr<Global>& ptr, globals)
    XML << xml::tag("Global")
	<< ptr
	<< xml::endtag("Global");
  
  XML << xml::endtag("Globals")
      << xml::tag("Locals");
  
  BOOST_FOREACH(const ClonePtr<Local>& ptr, locals)
    XML << xml::tag("Local")
	<< ptr
	<< xml::endtag("Local");
  
  XML << xml::endtag("Locals")
      << xml::tag("Interactions");
  
  BOOST_FOREACH(const ClonePtr<Interaction>& ptr, interactions)
    XML << xml::tag("Interaction")
	<< ptr
	<< xml::endtag("Interaction");
  
  XML << xml::endtag("Interactions")
      << xml::tag("Liouvillean")
      << p_liouvillean
      << xml::endtag("Liouvillean")
      << xml::endtag("Dynamics");
}

Iflt 
Dynamics::getLongestInteraction() const
{
  Iflt maxval = 0.0;

  BOOST_FOREACH(const ClonePtr<Interaction>& ptr, interactions)
    if (ptr->maxIntDist() > maxval)
      maxval = ptr->maxIntDist();

  return maxval;
}

void 
Dynamics::rescaleLengths(Iflt val)
{
  BOOST_FOREACH(ClonePtr<Interaction>& ptr, interactions)
    ptr->rescaleLengths(val);

  p_units->rescaleLength(val);
}

void 
Dynamics::SystemOverlapTest()
{
  p_liouvillean->updateAllParticles();

  std::vector<Particle>::const_iterator iPtr1, iPtr2;
  
  for (iPtr1 = Sim->particleList.begin(); iPtr1 != Sim->particleList.end(); ++iPtr1)
    for (iPtr2 = iPtr1 + 1; iPtr2 != Sim->particleList.end(); ++iPtr2)    
      getInteraction(*iPtr1, *iPtr2)->checkOverlaps(*iPtr1, *iPtr2);

  BOOST_FOREACH(const Particle& part, Sim->particleList)
    BOOST_FOREACH(const ClonePtr<Local>& lcl, locals)
    if (lcl->isInteraction(part))
      lcl->checkOverlaps(part);
    
}

void 
Dynamics::setUnits(Units* Uptr) 
{ p_units.set_ptr(Uptr); }

void 
Dynamics::setLiouvillean(Liouvillean* Uptr) 
{ p_liouvillean.set_ptr(Uptr); }
