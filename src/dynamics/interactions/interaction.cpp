/*  DYNAMO:- Event driven molecular dynamics simulator 
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

#include "interaction.hpp"
#include "include.hpp"
#include "../interactions/intEvent.hpp"
#include "../species/species.hpp"
#include "../../base/is_simdata.hpp"
#include <magnet/xmlreader.hpp>
#include <cstring>

Interaction::Interaction(DYNAMO::SimData* tmp, C2Range* nR):
  SimBase(tmp,"Interaction",IC_blue),
  range(nR)
{}

bool 
Interaction::isInteraction(const IntEvent &coll) const
{ 
  return isInteraction(Sim->particleList[coll.getParticle1ID()],
		       Sim->particleList[coll.getParticle2ID()]); 
}

bool 
Interaction::isInteraction(const Species &speci) const
{
  return !(intName.compare(speci.getIntName()));
}

xml::XmlStream& operator<<(xml::XmlStream& XML, 
			    const Interaction& g)
{
  g.outputXML(XML);
  return XML;
}

magnet::ClonePtr<C2Range>& 
Interaction::getRange() 
{ return range; }

const magnet::ClonePtr<C2Range>& 
Interaction::getRange() const
{ return range; }

Interaction*
Interaction::getClass(const magnet::xml::Node& XML, DYNAMO::SimData* Sim)
{
  if (!std::strcmp(XML.getAttribute("Type"),"HardSphere"))
    return new IHardSphere(XML, Sim);
  else if (!std::strcmp(XML.getAttribute("Type"),"RoughHardSphere"))
    return new IRoughHardSphere(XML, Sim);
  else if (!std::strcmp(XML.getAttribute("Type"),"SquareWell"))
    return new ISquareWell(XML, Sim);
  else if (!std::strcmp(XML.getAttribute("Type"),"SquareWellSeq"))
    return new ISWSequence(XML, Sim);
  else if (!std::strcmp(XML.getAttribute("Type"),"SquareBond"))
    return new ISquareBond(XML, Sim);
  else if (!std::strcmp(XML.getAttribute("Type"),"SoftCore"))
    return new ISoftCore(XML, Sim);
  else if (!std::strcmp(XML.getAttribute("Type"),"Null"))
    return new INull(XML, Sim);
  else if (!std::strcmp(XML.getAttribute("Type"),"Lines"))
    return new ILines(XML, Sim);
  else if (!std::strcmp(XML.getAttribute("Type"),"Dumbbells"))
    return new IDumbbells(XML, Sim);
  else if (!std::strcmp(XML.getAttribute("Type"),"RotatedParallelCubes"))
    return new IRotatedParallelCubes(XML, Sim);
  else if (!std::strcmp(XML.getAttribute("Type"),"Stepped"))
    return new IStepped(XML, Sim);
  else 
    M_throw() << XML.getAttribute("Type")
	      << ", Unknown type of interaction encountered";
}
