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

#include "interaction.hpp"
#include "include.hpp"
#include <boost/lexical_cast.hpp>
#include "../interactions/intEvent.hpp"
#include "../../extcode/xmlParser.h"
#include "../../base/is_exception.hpp"
#include "../species/species.hpp"
#include "../../base/is_simdata.hpp"
#include <cstring>

Interaction::Interaction(DYNAMO::SimData* tmp, C2Range* nR):
  SimBase(tmp,"Interaction",IC_blue),
  range(nR)
{}

bool 
Interaction::isInteraction(const IntEvent &coll) const
{ 
  return isInteraction(Sim->vParticleList[coll.getParticle1ID()],
		       Sim->vParticleList[coll.getParticle2ID()]); 
}

bool 
Interaction::isInteraction(const CSpecies &speci) const
{
  return !(intName.compare(speci.getIntName()));
}

xmlw::XmlStream& operator<<(xmlw::XmlStream& XML, 
			    const Interaction& g)
{
  g.outputXML(XML);
  return XML;
}

smrtPlugPtr<C2Range>& 
Interaction::getRange() 
{ return range; }

const smrtPlugPtr<C2Range>& 
Interaction::getRange() const
{ return range; }

Interaction*
Interaction::getClass(const XMLNode& XML, DYNAMO::SimData* Sim)
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
  else if (!std::strcmp(XML.getAttribute("Type"),"ParallelCubes"))
    return new IParallelCubes(XML, Sim);
  else if (!std::strcmp(XML.getAttribute("Type"),"RotatedParallelCubes"))
    return new IRotatedParallelCubes(XML, Sim);
  else if (!std::strcmp(XML.getAttribute("Type"),"Stepped"))
    return new IStepped(XML, Sim);
  else if (!std::strcmp(XML.getAttribute("Type"),"InfiniteMass"))
    return new IInfiniteMass(XML, Sim);
  else 
    D_throw() << XML.getAttribute("Type")
	      << ", Unknown type of interaction encountered";
}
