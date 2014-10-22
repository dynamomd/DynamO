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

#include <dynamo/locals/lcylinder.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/locals/localEvent.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/math/quaternion.hpp>
#ifdef DYNAMO_visualizer
# include <magnet/GL/objects/primitives/cylinder.hpp>
#endif

namespace dynamo {
  LCylinder::LCylinder(const magnet::xml::Node& XML, dynamo::Simulation* tmp):
    Local(tmp, "LocalCylinder")
  {
    operator<<(XML);
  }

  LocalEvent 
  LCylinder::getEvent(const Particle& part) const
  {
#ifdef ISSS_DEBUG
    if (!Sim->dynamics->isUpToDate(part))
      M_throw() << "Particle is not up to date";
#endif

    double colldist = 0.5 * _diameter->getProperty(part) + _cyl_radius;
    
    return LocalEvent(part, Sim->dynamics->getCylinderWallCollision(part, vPosition, vAxis, colldist), WALL, *this);
  }

  void
  LCylinder::runEvent(Particle& part, const LocalEvent& iEvent) const
  {
    ++Sim->eventCount;

    //Run the collision and catch the data

    NEventData EDat = Sim->dynamics->runCylinderWallCollision(part, vPosition, vAxis, _e->getProperty(part));
    Sim->_sigParticleUpdate(EDat);

    //Now we're past the event update the scheduler and plugins
    Sim->ptrScheduler->fullUpdate(part);
  
    for (shared_ptr<OutputPlugin> & Ptr : Sim->outputPlugins)
      Ptr->eventUpdate(iEvent, EDat);
  }

  void 
  LCylinder::operator<<(const magnet::xml::Node& XML)
  {
    range = shared_ptr<IDRange>(IDRange::getClass(XML.getNode("IDRange"),Sim));
    _diameter = Sim->_properties.getProperty(XML.getAttribute("ParticleDiameter"), Property::Units::Length());
    _cyl_radius = XML.getAttribute("CylinderRadius").as<double>() * Sim->units.unitLength();
    _e = Sim->_properties.getProperty(XML.getAttribute("Elasticity"), Property::Units::Dimensionless());
    localName = XML.getAttribute("Name");

    if (_diameter->getMinValue() + _cyl_radius == 0)
      M_throw() << "Cannot have a wall with a diameter and CylinderRadius of zero!";
      
    magnet::xml::Node xBrowseNode = XML.getNode("Axis");
    vAxis << xBrowseNode;
    if (vAxis.nrm() == 0)
      M_throw() << "The Axis for " << XML.getPath() << " named \"" << getName() << "\" has a length of 0. Cannot load";
    vAxis /= vAxis.nrm();
    xBrowseNode = XML.getNode("Origin");
    vPosition << xBrowseNode;
    vPosition *= Sim->units.unitLength();
  }

  void 
  LCylinder::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "Cylinder" 
	<< magnet::xml::attr("Name") << localName
	<< magnet::xml::attr("Elasticity") << _e->getName()
	<< magnet::xml::attr("ParticleDiameter") << _diameter->getName()
	<< magnet::xml::attr("CylinderRadius") << _cyl_radius / Sim->units.unitLength()
	<< range
	<< magnet::xml::tag("Axis")
	<< vAxis
	<< magnet::xml::endtag("Axis")
	<< magnet::xml::tag("Origin")
	<< vPosition / Sim->units.unitLength()
	<< magnet::xml::endtag("Origin");
  }

  bool 
  LCylinder::validateState(const Particle& part, bool textoutput) const
  {
    Vector pos(part.getPosition() - vPosition);
    Sim->BCs->applyBC(pos);
    
    double diam = 0.5 * _diameter->getProperty(part) + _cyl_radius;
    pos -= (pos | vAxis) * vAxis;
    double r = diam - pos.nrm();
    
    if (r > 0)
      {
	if (textoutput)
	  derr << "Particle " << part.getID() << " is " << r / Sim->units.unitLength() << " far into the cylindrical wall."
	       << "\nWall Position = " << (vPosition / Sim->units.unitLength()).toString()
	       << "\nWall Axis = " << vAxis.toString() << ", d = " << diam / Sim->units.unitLength()
	       << "\nParticle Position = " << (part.getPosition() / Sim->units.unitLength()).toString()
	       << "\nSeparation Vector = " << (pos / Sim->units.unitLength()).toString()
	       << "\nSeparation Distance = " << pos.nrm() / Sim->units.unitLength()
	       << std::endl;
	return true;
      }
    return false;
  }

#ifdef DYNAMO_visualizer
  shared_ptr<coil::RenderObj>
  LCylinder::getCoilRenderObj() const
  {
    if (!_renderObj)
      {
	using namespace magnet::GL::objects::primitives;
	const size_t LOD = 20;
	std::vector<GLfloat> verts = Cylinder::getVertices(LOD);
	magnet::math::Quaternion q = Quaternion::fromToVector(vAxis);
	
	const double axis_length = Sim->primaryCellSize[0] / Sim->units.unitLength();
	const double radius = 2 * _cyl_radius / Sim->units.unitLength();
	for (size_t i(0); i < verts.size() / 3; ++i)
	  {
	    Vector result = q * Vector{radius * verts[3*i+0], radius * verts[3*i+1], axis_length * verts[3*i+2]} + vPosition / Sim->units.unitLength();
	    verts[3*i+0] = result[0];
	    verts[3*i+1] = result[1];
	    verts[3*i+2] = result[2];
	  }

	_renderObj.reset(new coil::RTriangleMesh(getName(), verts, Cylinder::getIndices(LOD)));
      }
  
    return std::static_pointer_cast<coil::RenderObj>(_renderObj);
  }

  void 
  LCylinder::updateRenderData() const
  {}
#endif
}

