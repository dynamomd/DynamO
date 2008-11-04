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

#include "OrientationL.hpp"
#include "../2particleEventData.hpp"
#include <boost/progress.hpp>
#include "../../datatypes/vector.xml.hpp"

void
CLNOrientation::operator<<(const XMLNode& XML)
{
  XMLNode xSubNode(XML.getParentNode().getParentNode()
		   .getChildNode("ParticleData"));
  
  size_t nPart = xSubNode.nChildNode("Pt");
  
  orientationData.resize(nPart);
  
  I_cout() << "Loading orientation data....";

  boost::progress_display prog(nPart);
  
  int xml_iter = 0;

  for (size_t i = 0; i < nPart; ++i)
    {
      XMLNode xBrowseNode(xSubNode.getChildNode("Pt", &xml_iter));
      
      orientationData[i].orientation << xBrowseNode.getChildNode("U");
      orientationData[i].angularMomentum << xBrowseNode.getChildNode("O");
      ++prog;
    }
}

void
CLNOrientation::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") 
      << "NOrientation";
}

void 
CLNOrientation::outputExtraPDatXML(xmlw::XmlStream& XML,
				   const CParticle& part) const
{
  XML << xmlw::tag("O")
      << orientationData[part.getID()].angularMomentum
      << xmlw::endtag("O")
      << xmlw::tag("U")
      << orientationData[part.getID()].orientation
      << xmlw::endtag("U");
}

Iflt 
CLNOrientation::getLineLineCollision() const
{
  D_throw() << "Not implemented";
}

C2ParticleData 
CLNOrientation::runLineLineCollision() const
{
  D_throw() << "Not implemented";
}

void 
CLNOrientation::streamParticle(CParticle& part, const Iflt& dt) const
{
  //First hand over to the newtonian dynamics
  CLNewton::streamParticle(part, dt);
  
  //Now stream the orientation dynamics
}
