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

#include "wall.hpp"
#include "globEvent.hpp"
#include "../NparticleEventData.hpp"
#include "../../extcode/xmlParser.h"
#include "../../extcode/xmlwriter.hpp"
#include "../../datatypes/vector.xml.hpp"
#include <boost/lexical_cast.hpp>
#include "../../base/is_simdata.hpp"
#include "../../simulation/particle.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../units/units.hpp"

CGWall::CGWall(const XMLNode &XML, const DYNAMO::SimData* ptrSim):
  CGlobal(ptrSim),
  e(1.0)
{
  operator<<(XML);
}

CGWall::CGWall(const DYNAMO::SimData* nSim, Iflt ne, CVector<> nnorm, 
	       CVector<> norigin, std::string nname, CRange* nRange):
  CGlobal(nRange, nSim),
  vNorm(nnorm),
  vPosition(norigin),
  e(ne)
{
  globName = nname;
}

CGlobEvent 
CGWall::getEvent(const CParticle& part) const
{
  Sim->Dynamics.Liouvillean().updateParticle(part);

  return CGlobEvent(part, Sim->Dynamics.Liouvillean().getWallCollision(part, vPosition, vNorm), WALL, *this);
}

CNParticleData
CGWall::runEvent(const CGlobEvent& event) const
{
  return CNParticleData(Sim->Dynamics.Liouvillean().runWallCollision(event.getParticle(),vNorm,e));
}

void 
CGWall::initialise(size_t nID)
{
  ID=nID;
}

void 
CGWall::operator<<(const XMLNode& XML)
{
  range.set_ptr(CRange::loadClass(XML,Sim));
  
  try {
    e = boost::lexical_cast<Iflt>(XML.getAttribute("Elasticity"));
    XMLNode xBrowseNode = XML.getChildNode("Norm");
    globName = XML.getAttribute("Name");
    vNorm << xBrowseNode;
    vNorm = vNorm.unitVector();
    xBrowseNode = XML.getChildNode("Origin");
    vPosition << xBrowseNode;
    vPosition *= Sim->Dynamics.units().unitLength();
  } 
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in CGWall";
    }
}

void
CGWall::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "Wall" 
      << xmlw::attr("Name") << globName
      << xmlw::attr("Elasticity") << e
      << range
      << xmlw::tag("Norm")
      << vNorm
      << xmlw::endtag("Norm")
      << xmlw::tag("Origin")
      << vPosition / Sim->Dynamics.units().unitLength()
      << xmlw::endtag("Origin");
}
