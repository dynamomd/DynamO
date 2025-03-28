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
#include <dynamo/locals/lwall.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/units/units.hpp>
#include <magnet/overlap/cube_plane.hpp>

namespace dynamo {
LWall::LWall(const magnet::xml::Node &XML, dynamo::Simulation *tmp)
    : Local(tmp, "LocalWall") {
  operator<<(XML);
}

Event LWall::getEvent(const Particle &part) const {
#ifdef ISSS_DEBUG
  if (!Sim->dynamics->isUpToDate(part))
    M_throw() << "Particle is not up to date";
#endif

  double colldist = 0.5 * _diameter->getProperty(part);

  return Event(part,
               Sim->dynamics->getPlaneEvent(part, vPosition, vNorm, colldist),
               LOCAL, WALL, ID);
}

ParticleEventData LWall::runEvent(Particle &part, const Event &iEvent) const {
  ++Sim->eventCount;
  if (_amplitude > 0) {
    const double current_T =
        _sqrtT * _sqrtT +
        _amplitude * std::sin(_frequency * Sim->systemTime + _phase_offset);
    return Sim->dynamics->runAndersenWallCollision(
        part, vNorm, sqrt(current_T), _diameter->getProperty(part), _slip);
  }
  if (_sqrtT > 0)
    return Sim->dynamics->runAndersenWallCollision(
        part, vNorm, _sqrtT, _diameter->getProperty(part), _slip);
  else
    return Sim->dynamics->runPlaneEvent(part, vNorm, _e->getProperty(part),
                                        _diameter->getProperty(part));
}

void LWall::operator<<(const magnet::xml::Node &XML) {
  range = shared_ptr<IDRange>(IDRange::getClass(XML.getNode("IDRange"), Sim));
  _diameter = Sim->_properties.getProperty(XML.getAttribute("Diameter"),
                                           Property::Units::Length());

  if (_diameter->getMaxValue() == 0)
    M_throw() << "Cannot have a wall with a diameter of zero";

  _e = Sim->_properties.getProperty(XML.getAttribute("Elasticity"),
                                    Property::Units::Dimensionless());

  _sqrtT = 0;
  _amplitude = 0;
  _frequency = 0;
  _phase_offset = 0;

  if (XML.hasAttribute("Temperature")) {
    _sqrtT = sqrt(XML.getAttribute("Temperature").as<double>() *
                  Sim->units.unitEnergy());

    if (XML.hasAttribute("Amplitude") || XML.hasAttribute("Frequency") ||
        XML.hasAttribute("Phase_Offset")) {
      _amplitude =
          XML.getAttribute("Amplitude").as<double>() * Sim->units.unitEnergy();
      _frequency = XML.getAttribute("Frequency").as<double>();
      _phase_offset = XML.getAttribute("Phase_Offset").as<double>();
    }
    _slip = 0;
    if (XML.hasAttribute("Slip"))
      _slip = XML.getAttribute("Slip").as<double>();
  }

  if (_sqrtT < 0)
    M_throw() << "Cannot use negative temperatures on a Wall";

  if (_frequency < 0)
    M_throw() << "Cannot use negative frequencies on a Wall";

  if (_amplitude > _sqrtT * _sqrtT)
    M_throw() << "Amplitude of temperature oscillation cannot be bigger than "
                 "main temperature value";

  magnet::xml::Node xBrowseNode = XML.getNode("Norm");
  localName = XML.getAttribute("Name");
  vNorm << xBrowseNode;
  if (vNorm.nrm() == 0)
    M_throw() << "The normal for the Local Wall named \"" << getName()
              << "\" has a length of 0. Cannot load";
  vNorm /= vNorm.nrm();
  xBrowseNode = XML.getNode("Origin");
  vPosition << xBrowseNode;
  vPosition *= Sim->units.unitLength();
}

void LWall::outputXML(magnet::xml::XmlStream &XML) const {
  if (!_e)
    M_throw() << "e is unset!";
  XML << magnet::xml::attr("Type") << "Wall" << magnet::xml::attr("Name")
      << localName << magnet::xml::attr("Elasticity") << _e->getName()
      << magnet::xml::attr("Diameter") << _diameter->getName();

  if (_sqrtT > 0)
    XML << magnet::xml::attr("Temperature")
        << _sqrtT * _sqrtT / Sim->units.unitEnergy()
        << magnet::xml::attr("Slip") << _slip;

  if (_frequency > 0) {
    XML << magnet::xml::attr("Frequency") << _frequency;
    XML << magnet::xml::attr("Amplitude")
        << _amplitude / Sim->units.unitEnergy();
    XML << magnet::xml::attr("Phase_Offset")
        << std::fmod(_frequency * Sim->systemTime, 2 * M_PI);
  }

  XML << range << magnet::xml::tag("Norm") << vNorm
      << magnet::xml::endtag("Norm") << magnet::xml::tag("Origin")
      << vPosition / Sim->units.unitLength() << magnet::xml::endtag("Origin");
}

bool LWall::validateState(const Particle &part, bool textoutput) const {
  Vector pos(part.getPosition() - vPosition);
  Sim->BCs->applyBC(pos);

  double diam = 0.5 * _diameter->getProperty(part);
  double r = diam - std::abs(pos | vNorm);

  if (r > 0) {
    if (textoutput)
      derr << "Particle " << part.getID() << " is "
           << r / Sim->units.unitLength() << " far into the wall."
           << "\nWall Pos = "
           << Vector(vPosition / Sim->units.unitLength()).toString()
           << ", Normal = " << vNorm.toString()
           << ", d = " << diam / Sim->units.unitLength() << std::endl;
    return true;
  }
  return false;
}

#ifdef DYNAMO_visualizer

shared_ptr<coil::RenderObj> LWall::getCoilRenderObj() const {
  if (!_renderObj) {
    // Find out what the directions orthagonal to the norm are
    Vector orth1;
    for (size_t i(0); i < NDIM; ++i) {
      orth1 = Vector{0, 0, 0};
      orth1[i] = 1;
      orth1 = vNorm ^ orth1;
      if (orth1.nrm() != 0) {
        orth1 = orth1 / orth1.nrm();
        break;
      }
    }

    Vector orth2 = vNorm ^ orth1;
    if (orth2.nrm() == 0)
      M_throw() << "Cannot generate orthagonal vectors to plot LWall!";
    orth2 /= orth2.nrm();

    orth1 *= (orth1 | Sim->primaryCellSize);
    orth2 *= (orth2 | Sim->primaryCellSize);

    _renderObj.reset(new coil::RSurface(
        getName(), 10, vPosition - 0.5 * (orth1 + orth2), orth1, orth2, vNorm));
  }

  return std::static_pointer_cast<coil::RenderObj>(_renderObj);
}

void LWall::updateRenderData() const {}
#endif
} // namespace dynamo
