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

#include <dynamo/BC/BC.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/locals/boundary.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/units/units.hpp>

namespace dynamo {
LBoundary::LBoundary(const magnet::xml::Node &XML, dynamo::Simulation *tmp)
    : Local(tmp, "Boundary") {
  operator<<(XML);
}

Event LBoundary::getEvent(const Particle &part) const {
#ifdef ISSS_DEBUG
  if (!Sim->dynamics->isUpToDate(part))
    M_throw() << "Particle is not up to date";
#endif

  Event event(part, std::numeric_limits<float>::infinity(), LOCAL, NONE, ID);

  const double diameter = _diameter->getProperty(part);
  for (size_t objid(0); objid < _objects.size(); ++objid) {
    Event newevent = _objects[objid]->getEvent(part, diameter);
    if (newevent < event) {
      newevent._sourceID = ID;
      newevent._additionalData2 = objid;
      event = newevent;
    }
  }

  return event;
}

ParticleEventData LBoundary::runEvent(Particle &part,
                                      const Event &event) const {
  ++Sim->eventCount;
  const Vector normal =
      _objects[event._additionalData2]->getContactNormal(part, event);
  const double e = 1.0;
  const double diameter = _diameter->getProperty(part);
  return Sim->dynamics->runPlaneEvent(part, normal, e, diameter);
}

void LBoundary::operator<<(const magnet::xml::Node &XML) {
  range = shared_ptr<IDRange>(IDRange::getClass(XML.getNode("IDRange"), Sim));
  localName = XML.getAttribute("Name");

  _oscillationData._origin << XML.getNode("Origin");
  _oscillationData._origin *= Sim->units.unitLength();
  _diameter = Sim->_properties.getProperty(XML.getAttribute("Diameter"),
                                           Property::Units::Length());
  if (_diameter->getMaxValue() == 0)
    M_throw() << "Cannot have a boundary with a diameter of zero";

  _oscillationData._amplitude = Vector();
  _oscillationData._freq = 0;
  _oscillationData._t_shift = 0;

  const size_t data_count = XML.hasNode("Amplitude") +
                            XML.hasAttribute("Frequency") +
                            XML.hasAttribute("Phase");
  if ((data_count != 3) && (data_count != 0))
    M_throw() << "For oscillating walls you must have an Amplitude, Frequency, "
                 "and Phase specified."
              << XML.getPath();

  if (data_count == 3) {
    _oscillationData._freq =
        XML.getAttribute("Frequency").as<double>() / Sim->units.unitTime();
    _oscillationData._t_shift =
        XML.getAttribute("Phase").as<double>() * Sim->units.unitTime();
    _oscillationData._amplitude << XML.getNode("Amplitude");
    _oscillationData._amplitude *= Sim->units.unitLength();
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
    _objects.push_back(boundary::Object::getClass(node, Sim, _oscillationData));

  if (_objects.empty())
    M_throw() << "Boundary Locals must have at least one Object.\n"
              << XML.getPath();
}

void LBoundary::outputXML(magnet::xml::XmlStream &XML) const {
  XML << magnet::xml::attr("Type") << "Boundary" << magnet::xml::attr("Name")
      << localName << magnet::xml::attr("Diameter") << _diameter->getName();

  if (_kT > 0) // If the kT is >0 the wall is thermalised
    XML << magnet::xml::attr("kT") << _kT / Sim->units.unitEnergy();

  if (_oscillationData._freq !=
      0) //_freq is non-zero if the system is oscillating
    XML << magnet::xml::attr("Frequency")
        << _oscillationData._freq * Sim->units.unitTime()
        << magnet::xml::attr("Phase")
        << _oscillationData._t_shift / Sim->units.unitTime();

  XML << range;

  if (_oscillationData._freq != 0)
    XML << magnet::xml::tag("Amplitude")
        << _oscillationData._amplitude / Sim->units.unitLength()
        << magnet::xml::endtag("Amplitude");

  XML << magnet::xml::tag("Origin")
      << _oscillationData._origin / Sim->units.unitLength()
      << magnet::xml::endtag("Origin");

  for (const shared_ptr<boundary::Object> &obj : _objects) {
    XML << magnet::xml::tag("Object");
    obj->outputXML(XML);
    XML << magnet::xml::endtag("Object");
  }
}

bool LBoundary::validateState(const Particle &part, bool textoutput) const {
  for (const shared_ptr<boundary::Object> &obj : _objects)
    if (obj->validateState(part, textoutput))
      return true;
  return false;
}

#ifdef DYNAMO_visualizer
std::pair<std::vector<float>, std::vector<GLuint>>
LBoundary::getTessalatedSurfaces() const {

  typedef std::pair<std::vector<float>, std::vector<GLuint>> ReturnType;
  ReturnType retval;

  for (const shared_ptr<boundary::Object> &obj : _objects) {
    ReturnType val = obj->getTessalatedSurfaces();
    const size_t vertex_offset = retval.first.size() / 3;

    // Merge the vertices
    retval.first.insert(retval.first.end(), val.first.begin(), val.first.end());

    // Merge the indicies
    for (const auto &index : val.second)
      retval.second.push_back(index + vertex_offset);
  }
  return retval;
}

shared_ptr<coil::RenderObj> LBoundary::getCoilRenderObj() const {
  if (!_renderObj) {
    const auto triangles = getTessalatedSurfaces();
    _renderObj.reset(
        new coil::RTriangleMesh(getName(), triangles.first, triangles.second));
  }

  return std::static_pointer_cast<coil::RenderObj>(_renderObj);
}

void LBoundary::updateRenderData() const {
  if (_oscillationData._freq == 0)
    return;
  const auto triangles = getTessalatedSurfaces();
  _renderObj->updateGLData(triangles.first, triangles.second);
}
#endif
} // namespace dynamo
