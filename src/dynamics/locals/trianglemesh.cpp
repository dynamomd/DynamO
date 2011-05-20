/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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

#include "trianglemesh.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "localEvent.hpp"
#include "../NparticleEventData.hpp"
#include "../overlapFunc/CubePlane.hpp"
#include "../units/units.hpp"
#include "../../datatypes/vector.xml.hpp"
#include "../../schedulers/scheduler.hpp"


LTriangleMesh::LTriangleMesh(dynamo::SimData* nSim, double e, std::string name, CRange* nRange):
  Local(nRange, nSim, "LocalWall"),
  _e(e)
{
  localName = name;
}

LTriangleMesh::LTriangleMesh(const magnet::xml::Node& XML, dynamo::SimData* tmp):
  Local(tmp, "LocalWall")
{
  operator<<(XML);
}

LocalEvent 
LTriangleMesh::getEvent(const Particle& part) const
{
#ifdef ISSS_DEBUG
  if (!Sim->dynamics.getLiouvillean().isUpToDate(part))
    M_throw() << "Particle is not up to date";
#endif

  M_throw() << "Not implemented";
}

void
LTriangleMesh::runEvent(const Particle& part, const LocalEvent& iEvent) const
{
  M_throw() << "Not implemented";
}

bool 
LTriangleMesh::isInCell(const Vector & Origin, const Vector& CellDim) const
{ return true; }

void 
LTriangleMesh::initialise(size_t nID)
{ ID = nID; }

void 
LTriangleMesh::operator<<(const magnet::xml::Node& XML)
{
  range.set_ptr(CRange::getClass(XML,Sim));
  
  try {
    _e = XML.getAttribute("Elasticity").as<double>();
  } 
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in LTriangleMesh";
    }
}

void 
LTriangleMesh::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Type") << "TriangleMesh" 
      << xml::attr("Name") << localName
      << xml::attr("Elasticity") << _e
      << range;
}

void 
LTriangleMesh::checkOverlaps(const Particle& p1) const
{}
