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

#include "2RChain.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../../simulation/particle.hpp"

C2RChain::C2RChain(const XMLNode& XML, const DYNAMO::SimData*):
  range1(0),range2(0) 
{ 
  if (strcmp(XML.getAttribute("Range"),"Chain"))
    D_throw() << "Attempting to load a chain from a non chain";
  
  range1 = boost::lexical_cast<unsigned long>(XML.getAttribute("Start"));
  range2 = boost::lexical_cast<unsigned long>(XML.getAttribute("End"));
}

bool 
C2RChain::isInRange(const Particle&p1, const Particle&p2) const
{
  if (p1.getID() > p2.getID())
    {
      if ((p1.getID() - p2.getID()) == 1)
	if ((p2.getID() >= range1) && (p1.getID() <= range2))
	  return true;
    }
  else
    if ((p2.getID() - p1.getID()) == 1)
      if ((p1.getID() >= range1) && (p2.getID() <= range2))
	return true;
  

  return false;
}

void 
C2RChain::operator<<(const XMLNode&)
{
  D_throw() << "Due to problems with CRAll C2RChain operator<< cannot work for this class";
}

void 
C2RChain::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Range") << "Chain" 
      << xmlw::attr("Start")
      << range1
      << xmlw::attr("End")
      << range2;
}

