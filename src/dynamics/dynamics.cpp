/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
#include "../base/is_exception.hpp"
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

smrtPlugPtr<CTopology>& 
Dynamics::getTopology(std::string name)
{
  BOOST_FOREACH(smrtPlugPtr<CTopology>& sysPtr, topology)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  D_throw() << "Could not find the topology " << name;
}

const smrtPlugPtr<CTopology>& 
Dynamics::getTopology(std::string name) const
{
  BOOST_FOREACH(const smrtPlugPtr<CTopology>& sysPtr, topology)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  D_throw() << "Could not find the topology " << name;
}

const CSpecies& 
Dynamics::getSpecies(const CParticle& p1) const 
{
  BOOST_FOREACH(const smrtPlugPtr<CSpecies>& ptr, species)
    if (ptr->isSpecies(p1))
      return *ptr;
  
  D_throw() << "Could not find the requested species"
	    << "\nID = " << p1.getID();
}

xmlw::XmlStream& operator<<(xmlw::XmlStream& XML, 
			    const Dynamics& g)
{
  g.outputXML(XML);
  return XML;
}

const CSpecies& 
Dynamics::getSpecies(std::string name) const
{
  BOOST_FOREACH(const smrtPlugPtr<CSpecies>& ptr, species)
    if (ptr->getName() == name)
      return *ptr;
  
  D_throw() << "Could not find the " << name << " species"; 
}

smrtPlugPtr<CSystem>&
Dynamics::getSystem(std::string name)
{
  BOOST_FOREACH(smrtPlugPtr<CSystem>& sysPtr, systems)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  D_throw() << "Could not find system plugin";
}

const smrtPlugPtr<CSystem>&
Dynamics::getSystem(std::string name) const
{
  BOOST_FOREACH(const smrtPlugPtr<CSystem>& sysPtr, systems)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  D_throw() << "Could not find system plugin";
}

smrtPlugPtr<CGlobal>&
Dynamics::getGlobal(std::string name)
{
  BOOST_FOREACH(smrtPlugPtr<CGlobal>& sysPtr, globals)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  D_throw() << "Could not find global plugin";
}

const smrtPlugPtr<CGlobal>&
Dynamics::getGlobal(std::string name) const
{
  BOOST_FOREACH(const smrtPlugPtr<CGlobal>& sysPtr, globals)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  D_throw() << "Could not find global plugin";
}

smrtPlugPtr<CLocal>&
Dynamics::getLocal(std::string name)
{
  BOOST_FOREACH(smrtPlugPtr<CLocal>& sysPtr, locals)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  D_throw() << "Could not find local plugin";
}

const smrtPlugPtr<CLocal>&
Dynamics::getLocal(std::string name) const
{
  BOOST_FOREACH(const smrtPlugPtr<CLocal>& sysPtr, locals)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  D_throw() << "Could not find local plugin";
}

smrtPlugPtr<CInteraction>&
Dynamics::getInteraction(std::string name)
{
  BOOST_FOREACH(smrtPlugPtr<CInteraction>& sysPtr, interactions)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  D_throw() << "Could not find interaction plugin";
}

const smrtPlugPtr<CInteraction>&
Dynamics::getInteraction(std::string name) const
{
  BOOST_FOREACH(const smrtPlugPtr<CInteraction>& sysPtr, interactions)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  D_throw() << "Could not find interaction plugin";
}

void 
Dynamics::addSpecies(const smrtPlugPtr<CSpecies>& CSpe)
{
  if (Sim->status >= INITIALISED)
    D_throw() << "Cannot add species after simulation initialisation";

  species.push_back(CSpe);

  BOOST_FOREACH(smrtPlugPtr<CInteraction>& intPtr , interactions)
    {
      if (intPtr->isInteraction(*species.back()))
	{
	  species.back()->setIntPtr(intPtr.get_ptr());
	  break;
	}
    }
}

void 
Dynamics::addGlobal(CGlobal* newGlobal)
{
  if (Sim->status >= INITIALISED)
    D_throw() << "Cannot add global events after simulation initialisation";

  smrtPlugPtr<CGlobal> 
    tempPlug(newGlobal);
  
  globals.push_back(tempPlug);
}

void 
Dynamics::addLocal(CLocal* newLocal)
{
  if (Sim->status >= INITIALISED)
    D_throw() << "Cannot add local events after simulation initialisation";

  smrtPlugPtr<CLocal> 
    tempPlug(newLocal);
  
  locals.push_back(tempPlug);
}

void
Dynamics::addSystem(CSystem* newSystem)
{
  if (Sim->status >= INITIALISED)
    D_throw() << "Cannot add system events at this time, system is initialised";
  
  smrtPlugPtr<CSystem> 
    tempPlug(newSystem);
  
  systems.push_back(tempPlug); 
}

void
Dynamics::addStructure(CTopology* newSystem)
{ 
  if (Sim->status >= INITIALISED)
    D_throw() << "Cannot add structure after simulation initialisation";

  smrtPlugPtr<CTopology> 
    tempPlug(newSystem);
  
  topology.push_back(tempPlug); 
}

void 
Dynamics::addSystemTicker()
{
  if (Sim->status >= INITIALISED)
    D_throw() << "Cannot add the system ticker now";

  BOOST_FOREACH(smrtPlugPtr<CSystem>& ptr, systems)
    if (ptr->getName() == "SystemTicker")
      D_throw() << "System Ticker already exists";
  
    addSystem(new CSTicker(Sim, Sim->lastRunMFT, "SystemTicker"));
}

CInteraction* 
Dynamics::addInteraction(CInteraction* CInt)
{
  smrtPlugPtr<CInteraction> tempPlug(CInt);
  interactions.push_back(tempPlug);
  return interactions.back().get_ptr();
}

void 
Dynamics::initialise()
{
  BOOST_FOREACH(smrtPlugPtr<CSpecies>& ptr, species)
    ptr->initialise();
  
  unsigned int count = 0;
  //Now confirm that every species has only one species type!
  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    {
      BOOST_FOREACH(smrtPlugPtr<CSpecies>& ptr, species)
	if (ptr->isSpecies(part)) { count++; break; }
      
      if (count < 1)
	D_throw() << "Particle ID=" << part.getID() << " has no species";

      if (count > 1)
	D_throw() << "Particle ID=" << part.getID() << " has more than one species";
      count = 0;
    }

  //Now confirm that there are not more counts from each species than there are particles
  {
    unsigned long tot = 0;
    BOOST_FOREACH(smrtPlugPtr<CSpecies>& ptr, species)
      tot += ptr->getCount();
    
    if (tot < Sim->lN)
      D_throw() << "The particle count according to the species definition is too low\n"
		<< "discrepancy = " << tot - Sim->lN
		<< "\nN = " << Sim->lN;
    
    if (tot > Sim->lN)
      D_throw() << "The particle count according to the species definition is too high\n"
		<< "discrepancy = " << tot - Sim->lN
		<< "\nN = " << Sim->lN;
 }

  p_liouvillean->initialise();

  size_t ID=0;

  BOOST_FOREACH(smrtPlugPtr<CInteraction>& ptr, interactions)
    ptr->initialise(ID++);

  ID=0;

  //Must be initialised before globals. Neighbour lists are
  //implemented as globals and must initialise where locals are and their ID.
  BOOST_FOREACH(smrtPlugPtr<CLocal>& ptr, locals)
    ptr->initialise(ID++);

  ID=0;

  BOOST_FOREACH(smrtPlugPtr<CGlobal>& ptr, globals)
    ptr->initialise(ID++);

  ID=0;

  BOOST_FOREACH(smrtPlugPtr<CSystem>& ptr, systems)
    ptr->initialise(ID++);
}

const smrtPlugPtr<CInteraction>&
Dynamics::getInteraction(const CParticle& p1, const CParticle& p2) const 
{
  BOOST_FOREACH(const smrtPlugPtr<CInteraction>& ptr, interactions)
    if (ptr->isInteraction(p1,p2))
      return ptr;
  
  D_throw() << "Could not find the interaction requested";
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

  BOOST_FOREACH(smrtPlugPtr<CSystem>& ptr, systems)
    ptr->stream(dt);
}


Iflt
Dynamics::calcInternalEnergy() const
{
  Iflt intECurrent = 0.0;

  BOOST_FOREACH(const smrtPlugPtr<CInteraction> & plugptr, 
		Sim->dynamics.getInteractions())
    intECurrent += plugptr->getInternalEnergy();

  return intECurrent;
}

Iflt 
Dynamics::getNumberDensity() const
{
  return Sim->lN / Sim->dynamics.units().simVolume();
}

Iflt 
Dynamics::getPackingFraction() const
{
  Iflt volume = 0.0;
  
  BOOST_FOREACH(const smrtPlugPtr<CSpecies>& sp, Sim->dynamics.getSpecies())
    volume += pow(sp->getIntPtr()->hardCoreDiam(), NDIM) * sp->getCount();
  
  return  PI * volume / (6 * (Sim->dynamics.units().simVolume()));
}

void 
Dynamics::zeroMomentum(std::vector<CParticle> &pList)
{  
  Vector sumMV(0,0,0), velvec;
 
  //Determine the discrepancy VECTOR
  BOOST_FOREACH( CParticle & Part, pList)
    {
      Vector  pos(Part.getPosition()), vel(Part.getVelocity());
      BCs().applyBC(pos,vel);

      sumMV += vel * getSpecies(Part).getMass();
    }
  
  sumMV /= pList.size();
  
  BOOST_FOREACH(CParticle & Part, pList)
    Part.getVelocity() =  Part.getVelocity() - (sumMV / getSpecies(Part).getMass());
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
	  smrtPlugPtr<CTopology> tempPlug(CTopology::loadClass(xSubNode.getChildNode("Structure",i),Sim,i));
	  topology.push_back(tempPlug);
	}
    }
  
  xSubNode = xDynamics.getChildNode("Genus");  
  for (long i=0; i < xSubNode.nChildNode("Species"); i++)
    species.push_back(smrtPlugPtr<CSpecies>
		      (CSpecies::getClass(xSubNode.getChildNode("Species",i),Sim,i)));
  
  xSubNode = xDynamics.getChildNode("Liouvillean");
  p_liouvillean.set_ptr(Liouvillean::loadClass(xSubNode,Sim));  
  
  xSubNode = xDynamics.getChildNode("Interactions");
  for (long i=0; i < xSubNode.nChildNode("Interaction"); i++)
    {
      smrtPlugPtr<CInteraction> tempPlug(CInteraction::getClass
					 (xSubNode.getChildNode("Interaction",
								i),Sim));
      interactions.push_back(tempPlug);
    }  
  
  //Link the species and interactions
  BOOST_FOREACH(smrtPlugPtr<CSpecies>& sp , species)
    BOOST_FOREACH(smrtPlugPtr<CInteraction>& intPtr , interactions)
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
	  smrtPlugPtr<CGlobal> 
	    tempPlug(CGlobal::getClass(xSubNode.getChildNode("Global", i), Sim));
	  
	  globals.push_back(tempPlug);
	}
    }

  if (xDynamics.hasChild("Locals"))
    {
      xSubNode = xDynamics.getChildNode("Locals");
      
      for (long i = 0; i < xSubNode.nChildNode("Local"); ++i)
	{
	  smrtPlugPtr<CLocal> 
	    tempPlug(CLocal::getClass(xSubNode.getChildNode("Local", i), Sim));
	  
	  locals.push_back(tempPlug);
	}
    }
  
  if (xDynamics.hasChild("SystemEvents"))
    {
      xSubNode = xDynamics.getChildNode("SystemEvents");
      
      for (long i=0; i < xSubNode.nChildNode("System"); i++)
	{
	  smrtPlugPtr<CSystem> 
	    tempPlug(CSystem::getClass(xSubNode.getChildNode("System",i),Sim));
	  
	  systems.push_back(tempPlug);
	}
    }
}

void
Dynamics::outputXML(xmlw::XmlStream &XML) const
{
  XML << xmlw::tag("Dynamics")
      << xmlw::tag("Aspect_Ratio")
      << Sim->aspectRatio
      << xmlw::endtag("Aspect_Ratio")
      << xmlw::tag("Units")
      << p_units
      << xmlw::endtag("Units")
      << xmlw::tag("BC")
      << p_BC
      << xmlw::endtag("BC")
      << xmlw::tag("Genus");
  
  BOOST_FOREACH(const smrtPlugPtr<CSpecies>& ptr, species)
    XML << xmlw::tag("Species")
	<< ptr
	<< xmlw::endtag("Species");
  
  XML << xmlw::endtag("Genus")
      << xmlw::tag("Topology");
  
  BOOST_FOREACH(const smrtPlugPtr<CTopology>& ptr, topology)
    XML << xmlw::tag("Structure")
	<< ptr
	<< xmlw::endtag("Structure");
  
  XML << xmlw::endtag("Topology")
      << xmlw::tag("SystemEvents");
  
  BOOST_FOREACH(const smrtPlugPtr<CSystem>& ptr, systems)
    XML << ptr;
  
  XML << xmlw::endtag("SystemEvents")
      << xmlw::tag("Globals");
  
  BOOST_FOREACH(const smrtPlugPtr<CGlobal>& ptr, globals)
    XML << xmlw::tag("Global")
	<< ptr
	<< xmlw::endtag("Global");
  
  XML << xmlw::endtag("Globals")
      << xmlw::tag("Locals");
  
  BOOST_FOREACH(const smrtPlugPtr<CLocal>& ptr, locals)
    XML << xmlw::tag("Local")
	<< ptr
	<< xmlw::endtag("Local");
  
  XML << xmlw::endtag("Locals")
      << xmlw::tag("Interactions");
  
  BOOST_FOREACH(const smrtPlugPtr<CInteraction>& ptr, interactions)
    XML << xmlw::tag("Interaction")
	<< ptr
	<< xmlw::endtag("Interaction");
  
  XML << xmlw::endtag("Interactions")
      << xmlw::tag("Liouvillean")
      << p_liouvillean
      << xmlw::endtag("Liouvillean")
      << xmlw::endtag("Dynamics");
}

Iflt 
Dynamics::getLongestInteraction() const
{
  Iflt maxval = 0.0;

  BOOST_FOREACH(const smrtPlugPtr<CInteraction>& ptr, interactions)
    if (ptr->maxIntDist() > maxval)
      maxval = ptr->maxIntDist();

  return maxval;
}

void 
Dynamics::rescaleLengths(Iflt val)
{
  BOOST_FOREACH(smrtPlugPtr<CInteraction>& ptr, interactions)
    ptr->rescaleLengths(val);

  p_units->rescaleLength(val);
}

void 
Dynamics::SystemOverlapTest()
{
  p_liouvillean->updateAllParticles();

  std::vector<CParticle>::const_iterator iPtr1, iPtr2;
  
  for (iPtr1 = Sim->vParticleList.begin(); iPtr1 != Sim->vParticleList.end(); ++iPtr1)
    for (iPtr2 = iPtr1 + 1; iPtr2 != Sim->vParticleList.end(); ++iPtr2)    
      getInteraction(*iPtr1, *iPtr2)->checkOverlaps(*iPtr1, *iPtr2);
}

void 
Dynamics::setUnits(Units* Uptr) 
{ p_units.set_ptr(Uptr); }

void 
Dynamics::setLiouvillean(Liouvillean* Uptr) 
{ p_liouvillean.set_ptr(Uptr); }
