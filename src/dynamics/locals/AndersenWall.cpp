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
#include <cmath>

CLAndersenWall::CLAndersenWall(const XMLNode& XML, const DYNAMO::SimData* ptrSim):
  CLocal(ptrSim, "GlobalAndersenWall"),
  sqrtT(1.0)
{
  operator<<(XML);
}

CLAndersenWall::CLAndersenWall(const DYNAMO::SimData* nSim, Iflt nsqrtT,
			       CVector<> nnorm, CVector<> norigin, 
			       std::string nname, CRange* nRange):
  CLocal(nRange, nSim, "AndersenWall"),
  vNorm(nnorm),
  vPosition(norigin),
  sqrtT(nsqrtT)
{
  localName = nname;
}

CLocalEvent 
CLAndersenWall::getEvent(const CParticle& part) const
{
  Sim->Dynamics.Liouvillean().updateParticle(part);

  return CLocalEvent(part, Sim->Dynamics.Liouvillean().getWallCollision(part, vPosition, vNorm), WALL, *this);
}

CNParticleData
CLAndersenWall::runEvent(const CLocalEvent& event) const
{
  return CNParticleData
    (Sim->Dynamics.Liouvillean().runAndersenWallCollision
     (event.getParticle(), vNorm, sqrtT));
}

bool 
CLAndersenWall::isInCell(const CVector<>& Origin, 
			 const CVector<>& CellDim) const
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
		 * Sim->Dynamics.units().unitEnergy());

    XMLNode xBrowseNode = XML.getChildNode("Norm");
    localName = XML.getAttribute("Name");
    vNorm << xBrowseNode;
    vNorm = vNorm.unitVector();
    xBrowseNode = XML.getChildNode("Origin");
    vPosition << xBrowseNode;
    vPosition *= Sim->Dynamics.units().unitLength();

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
    / Sim->Dynamics.units().unitEnergy()
      << range
      << xmlw::tag("Norm")
      << vNorm
      << xmlw::endtag("Norm")
      << xmlw::tag("Origin")
      << vPosition / Sim->Dynamics.units().unitLength()
      << xmlw::endtag("Origin");
}
