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

#include <dynamo/locals/boundary.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>

namespace dynamo {
  LBoundary::LBoundary(const magnet::xml::Node& XML, dynamo::Simulation* tmp):
    Local(tmp, "Boundary")
  {
    operator<<(XML);
  }

  Event 
  LBoundary::getEvent(const Particle& part) const
  {
#ifdef ISSS_DEBUG
    if (!Sim->dynamics->isUpToDate(part))
      M_throw() << "Particle is not up to date";
#endif

    Event event(part, HUGE_VAL, LOCAL, NONE, ID);
    for (size_t objid(0); objid < _objects.size(); ++objid) {
      Event newevent = _objects[objid]->getEvent(part);
      if (newevent <  event) {
	newevent._sourceID = ID;
	newevent._additionalData2 = objid;
	event = newevent;
      }
    }
    return event;
  }

  ParticleEventData
  LBoundary::runEvent(Particle& part, const Event& event) const
  {
    ++Sim->eventCount;
    return _objects[event._additionalData2]->runEvent(part, event);
  }

  void 
  LBoundary::operator<<(const magnet::xml::Node& XML)
  {
    range = shared_ptr<IDRange>(IDRange::getClass(XML.getNode("IDRange"),Sim));
    localName = XML.getAttribute("Name");

    _origin << XML.getNode("Origin");
    _origin *= Sim->units.unitLength();
    _diameter = Sim->_properties.getProperty(XML.getAttribute("Diameter"), Property::Units::Length());
    if (_diameter->getMaxValue() == 0)
      M_throw() << "Cannot have a boundary with a diameter of zero";

    _amplitude = Vector();
    _freq = 0;
    _t_shift = 0;
    
    const size_t data_count = XML.hasNode("Amplitude") + XML.hasAttribute("Frequency") + XML.hasAttribute("Phase");
    if ((data_count < 3) && (data_count > 0))
      M_throw() << "For oscillating walls you must have an Amplitude, Frequency, and Phase specified."
		<< XML.getPath();

    if (data_count == 3) {
      _amplitude << XML.getNode("Amplitude");
      _amplitude *= Sim->units.unitLength();
      _freq = XML.getAttribute("Frequency").as<double>() / Sim->units.unitTime();
      _t_shift = XML.getAttribute("Phase").as<double>() * Sim->units.unitTime();
    }

    _kT = 0;
    if (XML.hasAttribute("kT")) {
      if (data_count == 3)
	M_throw() << "Cannot have both a thermalised wall and a oscillating wall" 
		  << XML.getPath();
      _kT = XML.getAttribute("kT").as<double>() * Sim->units.unitEnergy();
    }

    if (_kT < 0)
      M_throw() << "Temperature is less than zero" << XML.getPath();

    for (magnet::xml::Node node = XML.findNode("Object"); node.valid(); ++node)
      _objects.push_back(boundary::Object::getClass(node, Sim));
  }

  void 
  LBoundary::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "Boundary"
	<< magnet::xml::attr("Name") << localName
	<< magnet::xml::attr("Diameter") << _diameter->getName()
      ;

    if (_kT > 0)//If the kT is >0 the wall is thermalised
      XML << magnet::xml::attr("kT") << _kT / Sim->units.unitEnergy();

    if (_freq != 0) //_freq is non-zero if the system is oscillating
      XML << magnet::xml::attr("Frequency") << _freq * Sim->units.unitTime()
	  << magnet::xml::attr("Phase") << _t_shift / Sim->units.unitTime();
    
    XML << range;
    
    if (_freq != 0)
      XML << magnet::xml::tag("Amplitude")
	  << _amplitude / Sim->units.unitLength()
	  << magnet::xml::endtag("Amplitude")
	;

    XML << magnet::xml::tag("Origin")
	<< _origin / Sim->units.unitLength()
	<< magnet::xml::endtag("Origin");

    for (const shared_ptr<boundary::Object>& obj : _objects) {
      XML << magnet::xml::tag("Object");
      obj->outputXML(XML);
      XML << magnet::xml::endtag("Object");
    }
  }

  bool 
  LBoundary::validateState(const Particle& part, bool textoutput) const
  {
    for (const shared_ptr<boundary::Object>& obj : _objects)
      if (obj->validateState(part, textoutput))
	return true;
    return false;
  }

#ifdef DYNAMO_visualizer
  std::pair<std::vector<float>, std::vector<GLuint> > 
  LBoundary::getTessalatedSurfaces() const {
    std::vector<float> verts;
    std::vector<GLuint> elems;
    return std::make_pair(verts, elems);
  }

  shared_ptr<coil::RenderObj>
  LBoundary::getCoilRenderObj() const
  {
    if (!_renderObj) {
      auto triangles = getTessalatedSurfaces();
      _renderObj.reset(new coil::RTriangleMesh(getName(), triangles.first, triangles.second));
    }

    return std::static_pointer_cast<coil::RenderObj>(_renderObj);
  }

  void
  LBoundary::updateRenderData() const {
    //auto triangles = getTessalatedSurfaces();
    //_renderObj->updateGLData(triangles.first, triangles.second);
  }
#endif
}

