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

CDynamics::CDynamics(DYNAMO::SimData* tmp): 
  SimBase(tmp,"CDynamics",IC_purple),
  p_BC(NULL), 
  p_liouvillean(NULL),
  p_units(NULL)
{}

CDynamics::CDynamics(const XMLNode& XML, DYNAMO::SimData* tmp): 
  SimBase(tmp, "CDynamics", IC_purple),
  p_BC(NULL), 
  p_units(NULL)
{ operator<<(XML); }

CDynamics::~CDynamics() {}

std::vector<smrtPlugPtr<CTopology> >& 
CDynamics::getTopology()
{
  return topology;
}

const std::vector<smrtPlugPtr<CTopology> >& 
CDynamics::getTopology() const
{
  return topology;
}

smrtPlugPtr<CTopology>& 
CDynamics::getTopology(std::string name)
{
  BOOST_FOREACH(smrtPlugPtr<CTopology>& sysPtr, topology)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  D_throw() << "Could not find the topology " << name;
}

const smrtPlugPtr<CTopology>& 
CDynamics::getTopology(std::string name) const
{
  BOOST_FOREACH(const smrtPlugPtr<CTopology>& sysPtr, topology)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  D_throw() << "Could not find the topology " << name;
}

void 
CDynamics::setAspectRatio(const CVector<> &AR)
{ Sim->aspectRatio = AR;}

const std::vector<smrtPlugPtr<CInteraction> >& 
CDynamics::getInteractions() const
{ return interactions; }

std::vector<smrtPlugPtr<CInteraction> >& 
CDynamics::getInteractions()
{ return interactions; }

const std::vector<smrtPlugPtr<CGlobal> >& 
CDynamics::getGlobals() const
{ return globals; }

const std::vector<smrtPlugPtr<CLocal> >& 
CDynamics::getLocals() const
{ return locals; }

const std::vector<CSpecies> & 
CDynamics::getSpecies() const
{ return species; }

std::vector<smrtPlugPtr<CSystem> >& 
CDynamics::getSystemEvents()
{ return systems; }

const std::vector<smrtPlugPtr<CSystem> >& 
CDynamics::getSystemEvents() const
{ return systems; }

const CSpecies& 
CDynamics::getSpecies(const CParticle& p1) const 
{
  BOOST_FOREACH(const CSpecies& ptr, species)
    if (ptr.isSpecies(p1))
      return ptr;
  
  D_throw() << "Could not find the requested species"
	    << "\nID = " << p1.getID();
}

xmlw::XmlStream& operator<<(xmlw::XmlStream& XML, 
			    const CDynamics& g)
{
  g.outputXML(XML);
  return XML;
}

const CSpecies& 
CDynamics::getSpecies(std::string name) const
{
  BOOST_FOREACH(const CSpecies& ptr, species)
    if (ptr.getName() == name)
      return ptr;
  
  D_throw() << "Could not find the " << name << " species"; 
}

smrtPlugPtr<CSystem>&
CDynamics::getSystem(std::string name)
{
  BOOST_FOREACH(smrtPlugPtr<CSystem>& sysPtr, systems)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  D_throw() << "Could not find system plugin";
}

const smrtPlugPtr<CSystem>&
CDynamics::getSystem(std::string name) const
{
  BOOST_FOREACH(const smrtPlugPtr<CSystem>& sysPtr, systems)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  D_throw() << "Could not find system plugin";
}

smrtPlugPtr<CGlobal>&
CDynamics::getGlobal(std::string name)
{
  BOOST_FOREACH(smrtPlugPtr<CGlobal>& sysPtr, globals)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  D_throw() << "Could not find global plugin";
}

const smrtPlugPtr<CGlobal>&
CDynamics::getGlobal(std::string name) const
{
  BOOST_FOREACH(const smrtPlugPtr<CGlobal>& sysPtr, globals)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  D_throw() << "Could not find global plugin";
}

smrtPlugPtr<CLocal>&
CDynamics::getLocal(std::string name)
{
  BOOST_FOREACH(smrtPlugPtr<CLocal>& sysPtr, locals)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  D_throw() << "Could not find local plugin";
}

const smrtPlugPtr<CLocal>&
CDynamics::getLocal(std::string name) const
{
  BOOST_FOREACH(const smrtPlugPtr<CLocal>& sysPtr, locals)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  D_throw() << "Could not find local plugin";
}

smrtPlugPtr<CInteraction>&
CDynamics::getInteraction(std::string name)
{
  BOOST_FOREACH(smrtPlugPtr<CInteraction>& sysPtr, interactions)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  D_throw() << "Could not find interaction plugin";
}

const smrtPlugPtr<CInteraction>&
CDynamics::getInteraction(std::string name) const
{
  BOOST_FOREACH(const smrtPlugPtr<CInteraction>& sysPtr, interactions)
    if (sysPtr->getName() == name)
      return sysPtr;
  
  D_throw() << "Could not find interaction plugin";
}

void 
CDynamics::addSpecies(CSpecies CSpe)
{
  if (Sim->status >= INITIALISED)
    D_throw() << "Cannot add species after simulation initialisation";

  species.push_back(CSpe);

  BOOST_FOREACH(smrtPlugPtr<CInteraction>& intPtr , interactions)
    {
      if (intPtr->isInteraction(species.back()))
	{
	  species.back().setIntPtr(intPtr.get_ptr());
	  break;
	}
    }
}

void 
CDynamics::addGlobal(CGlobal* newGlobal)
{
  if (Sim->status >= INITIALISED)
    D_throw() << "Cannot add global events after simulation initialisation";

  smrtPlugPtr<CGlobal> 
    tempPlug(newGlobal);
  
  globals.push_back(tempPlug);
}

void
CDynamics::addSystem(CSystem* newSystem)
{
  if (Sim->status >= INITIALISED)
    D_throw() << "Cannot add system events at this time,"
      " try using addSystemLate"; 
  
  smrtPlugPtr<CSystem> 
    tempPlug(newSystem);
  
  systems.push_back(tempPlug); 
}

void
CDynamics::addSystemLate(CSystem* newSystem)
{
  if (Sim->status < INITIALISED)
    D_throw() << "Cannot add system events using this when early!"
      " try using addSystem"; 

  smrtPlugPtr<CSystem> 
    tempPlug(newSystem);
  
  systems.push_back(tempPlug); 

  systems.back()->initialise(systems.size()-1);
}

void
CDynamics::addStructure(CTopology* newSystem)
{ 
  if (Sim->status >= INITIALISED)
    D_throw() << "Cannot add structure after simulation initialisation";

  smrtPlugPtr<CTopology> 
    tempPlug(newSystem);
  
  topology.push_back(tempPlug); 
}

void 
CDynamics::addSystemTicker()
{
  BOOST_FOREACH(smrtPlugPtr<CSystem>& ptr, systems)
    if (ptr->getName() == "SystemTicker")
      D_throw() << "System Ticker already exists";

  if (Sim->status >= INITIALISED)
    addSystemLate(new CSTicker(Sim, Sim->lastRunMFT, "SystemTicker"));
  else
    addSystemLate(new CSTicker(Sim, Sim->lastRunMFT, "SystemTicker"));
}

void 
CDynamics::setUnits(CUnits* Uptr)
{
  //If set twice the plugin pointer will delete the old copy
  p_units.set_ptr(Uptr);
}

void 
CDynamics::setLiouvillean(CLiouvillean* Uptr)
{
  //If set twice the plugin pointer will delete the old copy
  p_liouvillean.set_ptr(Uptr);
}


CInteraction* 
CDynamics::addInteraction(CInteraction* CInt)
{
  smrtPlugPtr<CInteraction> tempPlug(CInt);
  interactions.push_back(tempPlug);
  return interactions.back().get_ptr();
}

void 
CDynamics::initialise()
{
  BOOST_FOREACH(CSpecies& ptr, species)
    ptr.initialise();
  
  unsigned int count = 0;
  //Now confirm that every species has only one species type!
  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    {
      BOOST_FOREACH(CSpecies& ptr, species)
	if (ptr.isSpecies(part)) count++;
      
      if (count < 1)
	D_throw() << "Particle ID=" << part.getID() << " has no species";

      if (count > 1)
	D_throw() << "Particle ID=" << part.getID() << " has more than one species";
      count = 0;
    }

  //Now confirm that there are not more counts from each species than there are particles
  {
    unsigned long tot = 0;
    BOOST_FOREACH(CSpecies& ptr, species)
      tot += ptr.getCount();
    
    if (tot < Sim->lN)
      D_throw() << "The particle count according to the species definition is too low\n"
		<< "discrepancy = " << tot - Sim->lN;
    
    if (tot > Sim->lN)
      D_throw() << "The particle count according to the species definition is too high\n"
		<< "discrepancy = " << tot - Sim->lN;
  }

  p_liouvillean->initialise();

  size_t ID=0;

  BOOST_FOREACH(smrtPlugPtr<CInteraction>& ptr, interactions)
    ptr->initialise(ID++);

  ID=0;

  BOOST_FOREACH(smrtPlugPtr<CGlobal>& ptr, globals)
    ptr->initialise(ID++);

  ID=0;

  BOOST_FOREACH(smrtPlugPtr<CLocal>& ptr, locals)
    ptr->initialise(ID++);

  ID=0;

  BOOST_FOREACH(smrtPlugPtr<CSystem>& ptr, systems)
    ptr->initialise(ID++);
}

CGlobEvent 
CDynamics::getEvent(const CParticle &p1) const
{
  CGlobEvent tmp, retval;

  BOOST_FOREACH(const smrtPlugPtr<CGlobal>& ptr, globals)
    if (ptr->isInteraction(p1))
      if (retval > (tmp = ptr->getEvent(p1)) )
	retval = tmp;

  return retval;
}

const smrtPlugPtr<CInteraction>&
CDynamics::getInteraction(const CParticle& p1, const CParticle& p2) const 
{
  BOOST_FOREACH(const smrtPlugPtr<CInteraction>& ptr, interactions)
    if (ptr->isInteraction(p1,p2))
      return ptr;
  
  D_throw() << "Could not find the interaction requested";
}

C2ParticleData
CDynamics::runEvent(const CIntEvent &coll)
{
  return coll.getInteraction().runCollision(coll);
}

CNParticleData
CDynamics::runEvent(const CGlobEvent &coll)
{
  return const_cast<CGlobEvent&>(coll).getGlobal().runEvent(coll);
}

CDynamics::CDynamics(const CDynamics &dyn):
  SimBase(dyn),
  p_BC(dyn.p_BC), 
  p_units(dyn.p_units)
{}

void 
CDynamics::stream(const Iflt& dt)
{
  p_BC->update(dt);

  p_liouvillean->stream(dt);

  BOOST_FOREACH(smrtPlugPtr<CSystem>& ptr, systems)
    ptr->stream(dt);
}

Iflt 
CDynamics::getKineticEnergy() const
{
  Iflt energy = 0.0;
  
  BOOST_FOREACH( const CParticle & part, Sim->vParticleList)
    energy += part.getVelocity().square() * getSpecies(part).getMass();
  
  return 0.5 * energy; 
}

CVector<> 
CDynamics::getVecKineticEnergy() const
{  
  CVector<> energy(0.0);
  
  BOOST_FOREACH( const CParticle & part, Sim->vParticleList)
    energy += part.getVelocity() * part.getVelocity()
    * Sim->Dynamics.getSpecies(part).getMass();
  
  return 0.5 * energy; 
}

Iflt
CDynamics::calcInternalEnergy() const
{
  Iflt intECurrent = 0.0;

  BOOST_FOREACH(const smrtPlugPtr<CInteraction> & plugptr, 
		Sim->Dynamics.getInteractions())
    intECurrent += plugptr->getInternalEnergy();

  return intECurrent;
}

Iflt 
CDynamics::getNumberDensity() const
{
  return Sim->lN / Sim->Dynamics.units().simVolume();
}

Iflt 
CDynamics::getPackingFraction() const
{
  Iflt volume = 0.0;
  
  BOOST_FOREACH(const CSpecies& sp, Sim->Dynamics.getSpecies())
    volume += pow(sp.getIntPtr()->hardCoreDiam(), NDIM) * sp.getCount();
  
  return  PI * volume / (6 * (Sim->Dynamics.units().simVolume()));
}

Iflt 
CDynamics::getParticleEnergy(const CParticle& part) const
{
  return 0.5 * (part.getVelocity().square()) * getSpecies(part).getMass();
}

Iflt 
CDynamics::getkT() const
{
  return 2.0 * getKineticEnergy() / (Sim->lN * static_cast<Iflt>(NDIM));
}
  
CVector<> 
CDynamics::getVeckT() const
{
  return getVecKineticEnergy() * (2.0 / Sim->lN);
}

void 
CDynamics::zeroMomentum(std::vector<CParticle> &pList)
{  
  CVector<> sumMV (0), velvec;
  
  //Determine the discrepancy VECTOR
  BOOST_FOREACH( CParticle & Part, pList)
    sumMV += Part.getVelocity() * getSpecies(Part).getMass();
  
  sumMV /= pList.size();
  
  BOOST_FOREACH(CParticle & Part, pList)
    Part.getVelocity() =  Part.getVelocity() - (sumMV / getSpecies(Part).getMass());
}

void
CDynamics::operator<<(const XMLNode& XML)
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
  p_units.set_ptr(CUnits::loadUnits(xSubNode,Sim));
  
  //Now load the BC part, after the aspect ratio!
  xSubNode = xDynamics.getChildNode("BC");
  p_BC.set_ptr(CBC::loadClass(xSubNode, Sim));
  
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
    species.push_back(CSpecies(xSubNode.getChildNode("Species",i),Sim,i));
  
  xSubNode = xDynamics.getChildNode("Liouvillean");
  p_liouvillean.set_ptr(CLiouvillean::loadClass(xSubNode,Sim));  
  
  xSubNode = xDynamics.getChildNode("Interactions");
  for (long i=0; i < xSubNode.nChildNode("Interaction"); i++)
    {
      smrtPlugPtr<CInteraction> tempPlug(CInteraction::getClass
					 (xSubNode.getChildNode("Interaction",
								i),Sim));
      interactions.push_back(tempPlug);
    }  
  
  //Link the species and interactions
  BOOST_FOREACH(CSpecies& sp , species)
    BOOST_FOREACH(smrtPlugPtr<CInteraction>& intPtr , interactions)
    if (intPtr->isInteraction(sp))
      {
	sp.setIntPtr(intPtr.get_ptr());
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
CDynamics::outputXML(xmlw::XmlStream &XML) const
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
      << xmlw::tag("Liouvillean")
      << p_liouvillean
      << xmlw::endtag("Liouvillean")
      << xmlw::tag("Genus");
  
  BOOST_FOREACH(const CSpecies& ptr, species)
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
      << xmlw::endtag("Dynamics");
}


CVector<>
CDynamics::getLabVelocity(const CParticle &part) const
{
  return part.getVelocity();
}

Iflt 
CDynamics::getLongestInteraction() const
{
  Iflt maxval = 0.0;

  BOOST_FOREACH(const smrtPlugPtr<CInteraction>& ptr, interactions)
    if (ptr->maxIntDist() > maxval)
      maxval = ptr->maxIntDist();

  return maxval;
}

void 
CDynamics::rescaleLengths(Iflt val)
{
  BOOST_FOREACH(smrtPlugPtr<CInteraction>& ptr, interactions)
    ptr->rescaleLengths(val);

  p_units->rescaleLength(val);
}

void 
CDynamics::SystemOverlapTest()
{
  p_liouvillean->updateAllParticles();

  std::vector<CParticle>::const_iterator iPtr1, iPtr2;
  
  for (iPtr1 = Sim->vParticleList.begin(); iPtr1 != Sim->vParticleList.end(); ++iPtr1)
    for (iPtr2 = iPtr1 + 1; iPtr2 != Sim->vParticleList.end(); ++iPtr2)    
      getInteraction(*iPtr1, *iPtr2)->checkOverlaps(*iPtr1, *iPtr2);
}
