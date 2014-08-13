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

#include <dynamo/interactions/parallelcubes.hpp>
#include <dynamo/interactions/intEvent.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/2particleEventData.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/ranges/IDRange.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/xmlwriter.hpp>
#include <sstream>
#include <cmath>
#include <iomanip>

namespace dynamo {
  IParallelCubes::IParallelCubes(const magnet::xml::Node& XML, 
					       dynamo::Simulation* tmp):
    Interaction(tmp,NULL)
  { operator<<(XML); }

  std::array<double, 4> IParallelCubes::getGlyphSize(size_t ID) const
  {
    return {{_diameter->getProperty(ID), 0, 0, 0}};
  }

  void 
  IParallelCubes::operator<<(const magnet::xml::Node& XML)
  { 
    Interaction::operator<<(XML);
    _diameter = Sim->_properties.getProperty(XML.getAttribute("Diameter"), Property::Units::Length());
    _e = Sim->_properties.getProperty(XML.getAttribute("Elasticity"), Property::Units::Dimensionless());
  }

  double 
  IParallelCubes::maxIntDist() const 
  { return std::sqrt(double(NDIM)) * _diameter->getMaxValue(); }

  double 
  IParallelCubes::getExcludedVolume(size_t ID) const 
  {
    double diam = _diameter->getProperty(ID);
    return diam * diam * diam; 
  }

  IntEvent 
  IParallelCubes::getEvent(const Particle &p1, const Particle &p2) const 
  { 
#ifdef DYNAMO_DEBUG
    if (!Sim->dynamics->isUpToDate(p1))
      M_throw() << "Particle 1 is not up to date";
  
    if (!Sim->dynamics->isUpToDate(p2))
      M_throw() << "Particle 2 is not up to date";

    if (p1 == p2)
      M_throw() << "You shouldn't pass p1==p2 events to the interactions!";
#endif 

    const double dt = Sim->dynamics->CubeCubeInRoot(p1, p2, _diameter->getProperty(p1, p2));

    if (dt != HUGE_VAL)
      return IntEvent(p1, p2, dt, CORE, *this);
  
    return IntEvent(p1, p2, HUGE_VAL, NONE, *this);
  }

  PairEventData
  IParallelCubes::runEvent(Particle& p1, Particle& p2, const IntEvent& iEvent)
  {
    ++Sim->eventCount;
    return Sim->dynamics->parallelCubeColl(iEvent, _e->getProperty(p1, p2), _diameter->getProperty(p1, p2));
  }
   
  void 
  IParallelCubes::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "ParallelCubes"
	<< magnet::xml::attr("Diameter") << _diameter->getName()
	<< magnet::xml::attr("Elasticity") << _e->getName()
	<< magnet::xml::attr("Name") << intName
	<< range;
  }


  bool
  IParallelCubes::validateState(const Particle& p1, const Particle& p2, bool textoutput) const
  {
    const double d = _diameter->getProperty(p1, p2);

    if (Sim->dynamics->cubeOverlap(p1, p2, d))
      {
	Vector rij = p1.getPosition() - p2.getPosition();
	Sim->BCs->applyBC(rij);
	rij /= Sim->units.unitLength();

	if (textoutput)
	  derr << "Particle " << p1.getID() << " and Particle " << p2.getID()
	       << " have a separation of " << rij.toString() 
	       << " but they are cubes and cannot be closer than " << d / Sim->units.unitLength()
	       << " in any dimension."
	       << std::endl;

	return true;
      }

    return false;
  }
}

