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

#include <dynamo/systems/rotateGravity.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/dynamics/gravity.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/particle.hpp>
#include <dynamo/species/species.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/ranges/include.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/math/quaternion.hpp>
#include <fstream>

namespace dynamo {
  SysRotateGravity::SysRotateGravity(const magnet::xml::Node& XML, dynamo::Simulation* tmp): 
    System(tmp)
  {
    operator<<(XML);
    type = ROTATEGRAVITY;
  }

  SysRotateGravity::SysRotateGravity(dynamo::Simulation* tmp, std::string name, double timestep, double angularvel, Vector axis):
    System(tmp),
    _timestep(timestep),
    _angularvel(angularvel),
    _rotationaxis(axis)
  {
    type = ROTATEGRAVITY;
    sysName = name;
  }

  NEventData
  SysRotateGravity::runEvent()
  {
    NEventData SDat;
    for (const shared_ptr<Species>& species : Sim->species)
      for (const unsigned long& partID : *species->getRange())
      SDat.L1partChanges.push_back(ParticleEventData(Sim->particles[partID], *species, RECALCULATE));
    Sim->dynamics->updateAllParticles();
    
    shared_ptr<DynGravity> dynamics = std::dynamic_pointer_cast<DynGravity>(Sim->dynamics);
    if (!dynamics)
      M_throw() << "The RotateGravity system can only be used with the Gravity type dynamics";
    
    double g = dynamics->getGravityVector().nrm();
    Vector newg = magnet::math::Quaternion::fromAngleAxis(_angularvel * _timestep, _rotationaxis) *  dynamics->getGravityVector();
    dynamics->setGravityVector(newg.normal() * g);

    dt = _timestep;
    return SDat;
  }

  void 
  SysRotateGravity::initialise(size_t nID)
  {
    ID = nID;
    dt = _timestep;
    shared_ptr<DynGravity> dynamics = std::dynamic_pointer_cast<DynGravity>(Sim->dynamics);
    if (!dynamics)
      M_throw() << "The RotateGravity system can only be used with the Gravity type dynamics";
  }

  void 
  SysRotateGravity::operator<<(const magnet::xml::Node& XML)
  {
    _angularvel = XML.getAttribute("AngularVel").as<double>() / Sim->units.unitTime();
    _timestep = XML.getAttribute("TimeStep").as<double>() * Sim->units.unitTime();
    _rotationaxis << XML.getNode("Axis");
    _rotationaxis.normalise();
    sysName = XML.getAttribute("Name");
  }

  void 
  SysRotateGravity::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::tag("System")
	<< magnet::xml::attr("Type") << "RotateGravity"
	<< magnet::xml::attr("Name") << sysName
	<< magnet::xml::attr("AngularVel") << _angularvel * Sim->units.unitTime();
  
    if (_timestep != std::numeric_limits<float>::infinity())
      XML << magnet::xml::attr("TimeStep") << _timestep / Sim->units.unitTime();

    XML << magnet::xml::tag("Axis") 
	<< _rotationaxis
	<< magnet::xml::endtag("Axis");

    XML<< magnet::xml::endtag("System");
  }
}
