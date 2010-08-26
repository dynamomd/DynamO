/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "AndersenWall.hpp"
#include "../overlapFunc/CubePlane.hpp"
#include "localEvent.hpp"
#include "../NparticleEventData.hpp"
#include "../../extcode/xmlParser.h"
#include "../../extcode/xmlwriter.hpp"
#include "../../datatypes/vector.xml.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_real.hpp>
#include "../../extcode/include/boost/random/normal_distribution.hpp"
#include "../../base/is_simdata.hpp"
#include "../../simulation/particle.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../units/units.hpp"
#include "../../schedulers/scheduler.hpp"
#include <cmath>

CLAndersenWall::CLAndersenWall(const XMLNode& XML, DYNAMO::SimData* ptrSim):
  Local(ptrSim, "GlobalAndersenWall"),
  sqrtT(1.0)
{
  operator<<(XML);
}

CLAndersenWall::CLAndersenWall(DYNAMO::SimData* nSim, Iflt nsqrtT,
			       Vector  nnorm, Vector norigin, 
			       std::string nname, CRange* nRange):
  Local(nRange, nSim, "AndersenWall"),
  vNorm(nnorm),
  vPosition(norigin),
  sqrtT(nsqrtT)
{
  localName = nname;
}

LocalEvent 
CLAndersenWall::getEvent(const Particle& part) const
{
#ifdef ISSS_DEBUG
  if (!Sim->dynamics.getLiouvillean().isUpToDate(part))
    D_throw() << "Particle is not up to date";
#endif

  return LocalEvent(part, Sim->dynamics.getLiouvillean().getWallCollision(part, vPosition, vNorm), WALL, *this);
}

void
CLAndersenWall::runEvent(const Particle& part, const LocalEvent& iEvent) const
{
  ++Sim->lNColl;
  
  NEventData EDat
    (Sim->dynamics.getLiouvillean().runAndersenWallCollision
     (part, vNorm, sqrtT));
  
  Sim->signalParticleUpdate(EDat);
  
  Sim->ptrScheduler->fullUpdate(part);
  
  BOOST_FOREACH(smrtPlugPtr<OutputPlugin> & Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(iEvent, EDat);
}

bool 
CLAndersenWall::isInCell(const Vector & Origin, 
			 const Vector & CellDim) const
{
  return DYNAMO::OverlapFunctions::CubePlane
    (Origin, CellDim, vPosition, vNorm);
}

void 
CLAndersenWall::initialise(size_t nID)
{
  ID = nID;
}

void 
CLAndersenWall::operator<<(const XMLNode& XML)
{
  range.set_ptr(CRange::loadClass(XML,Sim));
  
  try {
    
    sqrtT = sqrt(boost::lexical_cast<Iflt>(XML.getAttribute("Temperature")) 
		 * Sim->dynamics.units().unitEnergy());

    XMLNode xBrowseNode = XML.getChildNode("Norm");
    localName = XML.getAttribute("Name");
    vNorm << xBrowseNode;
    vNorm /= vNorm.nrm();
    xBrowseNode = XML.getChildNode("Origin");
    vPosition << xBrowseNode;
    vPosition *= Sim->dynamics.units().unitLength();

  } 
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in CLAndersenWall";
    }
}

void
CLAndersenWall::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "AndersenWall"
      << xmlw::attr("Name") << localName
      << xmlw::attr("Temperature") << sqrtT * sqrtT 
    / Sim->dynamics.units().unitEnergy()
      << range
      << xmlw::tag("Norm")
      << vNorm
      << xmlw::endtag("Norm")
      << xmlw::tag("Origin")
      << vPosition / Sim->dynamics.units().unitLength()
      << xmlw::endtag("Origin");
}

void 
CLAndersenWall::write_povray_info(std::ostream& os) const
{
  os << "object {\n plane {\n  <" << vNorm[0] << ", " << vNorm[1] 
     << ", " << vNorm[2] << ">, 0 texture{pigment { color rgb<0.5,0.5,0.5>}}}\n clipped_by{box {\n  <" << -Sim->aspectRatio[0]/2 
     << ", " << -Sim->aspectRatio[1]/2 << ", " << -Sim->aspectRatio[2]/2 
     << ">, <" << Sim->aspectRatio[0]/2 << ", " << Sim->aspectRatio[1]/2 
     << ", " << Sim->aspectRatio[2]/2 << "> }\n}\n translate <" << vPosition[0] << 
    ","<< vPosition[1] << "," << vPosition[2] << ">\n}\n";
}
