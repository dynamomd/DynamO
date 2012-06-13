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

#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/dynamics/include.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/dynamics/NparticleEventData.hpp>
#include <dynamo/dynamics/systems/sysTicker.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/exception.hpp>
#include <magnet/xmlwriter.hpp>
#include <boost/foreach.hpp>
#include <cmath>

namespace dynamo {

  Dynamics::Dynamics(dynamo::SimData* tmp): 
    SimBase(tmp, "Dynamics")
  {}

  Dynamics::Dynamics(const magnet::xml::Node& XML, dynamo::SimData* tmp): 
    SimBase(tmp, "Dynamics")
  { operator<<(XML); }

  Dynamics::~Dynamics() {}

  magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream& XML, 
				     const Dynamics& g)
  {
    g.outputXML(XML);
    return XML;
  }

  shared_ptr<System>&
  Dynamics::getSystem(std::string name)
  {
    BOOST_FOREACH(shared_ptr<System>& sysPtr, systems)
      if (sysPtr->getName() == name) 
	return sysPtr;
  
    M_throw() << "Could not find system plugin called " << name;
  }

  const shared_ptr<System>&
  Dynamics::getSystem(std::string name) const
  {
    BOOST_FOREACH(const shared_ptr<System>& sysPtr, systems)
      if (sysPtr->getName() == name)
	return sysPtr;
  
    M_throw() << "Could not find system plugin called " << name;
  }

  shared_ptr<Global>&
  Dynamics::getGlobal(std::string name)
  {
    BOOST_FOREACH(shared_ptr<Global>& sysPtr, globals)
      if (sysPtr->getName() == name)
	return sysPtr;
  
    M_throw() << "Could not find global plugin";
  }

  const shared_ptr<Global>&
  Dynamics::getGlobal(std::string name) const
  {
    BOOST_FOREACH(const shared_ptr<Global>& sysPtr, globals)
      if (sysPtr->getName() == name)
	return sysPtr;
  
    M_throw() << "Could not find global plugin";
  }

  shared_ptr<Local>&
  Dynamics::getLocal(std::string name)
  {
    BOOST_FOREACH(shared_ptr<Local>& sysPtr, locals)
      if (sysPtr->getName() == name)
	return sysPtr;
  
    M_throw() << "Could not find local plugin";
  }

  const shared_ptr<Local>&
  Dynamics::getLocal(std::string name) const
  {
    BOOST_FOREACH(const shared_ptr<Local>& sysPtr, locals)
      if (sysPtr->getName() == name)
	return sysPtr;
  
    M_throw() << "Could not find local plugin";
  }

  void Dynamics::addGlobal(shared_ptr<Global> ptr)
  {
    if (!ptr) M_throw() << "Cannot add an unset Global";
    if (Sim->status >= INITIALISED)
      M_throw() << "Cannot add global events after simulation initialisation";
    globals.push_back(ptr);
  }
  
  void Dynamics::addLocal(shared_ptr<Local> ptr)
  {
    if (!ptr) M_throw() << "Cannot add an unset Local";
    if (Sim->status >= INITIALISED)
      M_throw() << "Cannot add local events after simulation initialisation";
    locals.push_back(ptr);
  }
  
  void Dynamics::addSystem(shared_ptr<System> ptr)
  {
    if (!ptr) M_throw() << "Cannot add an unset System";
    if (Sim->status >= INITIALISED)
      M_throw() << "Cannot add system events at this time, system is initialised";
    systems.push_back(ptr); 
  }
    
  void 
  Dynamics::addSystemTicker()
  {
    BOOST_FOREACH(shared_ptr<System>& ptr, systems)
      if (ptr->getName() == "SystemTicker")
	M_throw() << "System Ticker already exists";
  
    addSystem(shared_ptr<System>(new SysTicker(Sim, Sim->lastRunMFT, "SystemTicker")));
  }

  void 
  Dynamics::initialise()
  {
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

  void 
  Dynamics::stream(const double& dt)
  {
    Sim->BCs->update(dt);

    Sim->liouvillean->stream(dt);

    BOOST_FOREACH(shared_ptr<System>& ptr, systems)
      ptr->stream(dt);
  }


  double
  Dynamics::calcInternalEnergy() const
  {
    double intECurrent = 0.0;

    BOOST_FOREACH(const shared_ptr<Interaction> & plugptr, 
		  Sim->interactions)
      intECurrent += plugptr->getInternalEnergy();

    return intECurrent;
  }

  double
  Dynamics::getSimVolume() const
  { 
    double vol = 1.0;
    for (size_t iDim = 0; iDim < NDIM; iDim++)
      vol *= Sim->primaryCellSize[iDim];
    return vol;
  }


  double
  Dynamics::getNumberDensity() const
  {
    return Sim->N / getSimVolume();
  }

  double 
  Dynamics::getPackingFraction() const
  {
    double volume = 0.0;
  
    BOOST_FOREACH(const shared_ptr<Species>& sp, Sim->species)
      BOOST_FOREACH(const size_t& ID, *(sp->getRange()))
      volume += sp->getIntPtr()->getExcludedVolume(ID);
  
    return  volume / getSimVolume();
  }

  void 
  Dynamics::setCOMVelocity(const Vector COMVelocity)
  {  
    Vector sumMV(0,0,0);
 
    long double sumMass(0);

    //Determine the discrepancy VECTOR
    BOOST_FOREACH(Particle & Part, Sim->particleList)
      {
	Vector  pos(Part.getPosition()), vel(Part.getVelocity());
	Sim->BCs->applyBC(pos,vel);
	double mass = Sim->species[Part].getMass(Part.getID());
	//Note we sum the negatives!
	sumMV -= vel * mass;
	sumMass += mass;
      }
  
    sumMV /= sumMass;
  
    sumMV += COMVelocity;

    BOOST_FOREACH(Particle & Part, Sim->particleList)
      Part.getVelocity() =  Part.getVelocity() + sumMV;
  }

  void
  Dynamics::operator<<(const magnet::xml::Node& XML)
  {  
    if (XML.hasNode("Globals"))
      for (magnet::xml::Node node = XML.getNode("Globals").fastGetNode("Global"); 
	   node.valid(); ++node)
	globals.push_back(Global::getClass(node, Sim));

    if (XML.hasNode("Locals"))
      for (magnet::xml::Node node = XML.getNode("Locals").fastGetNode("Local"); 
	   node.valid(); ++node)
	locals.push_back(Local::getClass(node, Sim));
  
    if (XML.hasNode("SystemEvents"))
      for (magnet::xml::Node node = XML.getNode("SystemEvents").fastGetNode("System"); 
	   node.valid(); ++node)
	systems.push_back(System::getClass(node, Sim));
  }

  void
  Dynamics::outputXML(magnet::xml::XmlStream &XML) const
  {
    XML << magnet::xml::tag("SystemEvents");
  
    BOOST_FOREACH(const shared_ptr<System>& ptr, systems)
      XML << *ptr;
  
    XML << magnet::xml::endtag("SystemEvents")
	<< magnet::xml::tag("Globals");
  
    BOOST_FOREACH(const shared_ptr<Global>& ptr, globals)
      XML << *ptr;
  
    XML << magnet::xml::endtag("Globals")
	<< magnet::xml::tag("Locals");
  
    BOOST_FOREACH(const shared_ptr<Local>& ptr, locals)
      XML << magnet::xml::tag("Local")
	  << *ptr
	  << magnet::xml::endtag("Local");
  
    XML << magnet::xml::endtag("Locals");
  }

  void 
  Dynamics::SystemOverlapTest()
  {
    Sim->liouvillean->updateAllParticles();

    std::vector<Particle>::const_iterator iPtr1, iPtr2;
  
    for (iPtr1 = Sim->particleList.begin(); iPtr1 != Sim->particleList.end(); ++iPtr1)
      for (iPtr2 = iPtr1 + 1; iPtr2 != Sim->particleList.end(); ++iPtr2)    
	Sim->getInteraction(*iPtr1, *iPtr2)->checkOverlaps(*iPtr1, *iPtr2);

    BOOST_FOREACH(const Particle& part, Sim->particleList)
      BOOST_FOREACH(const shared_ptr<Local>& lcl, locals)
      if (lcl->isInteraction(part))
	lcl->checkOverlaps(part);
    
  }
}
