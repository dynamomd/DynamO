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

#include <dynamo/interactions/interaction.hpp>
#include <dynamo/interactions/include.hpp>
#include <dynamo/interactions/intEvent.hpp>
#include <dynamo/species/species.hpp>
#include <dynamo/simulation.hpp>
#include <magnet/xmlreader.hpp>
#include <cstring>

namespace dynamo {
  Interaction::Interaction(dynamo::Simulation* tmp, C2Range* nR):
    SimBase(tmp, "Interaction"),
    range(nR)
  {}

  void 
  Interaction::operator<<(const magnet::xml::Node& XML)
  { range = shared_ptr<C2Range>(C2Range::getClass(XML,Sim)); }

  bool 
  Interaction::isInteraction(const IntEvent &coll) const
  { 
    return isInteraction(Sim->particles[coll.getParticle1ID()],
			 Sim->particles[coll.getParticle2ID()]); 
  }

  bool 
  Interaction::isInteraction(const Species &speci) const
  {
    return !(intName.compare(speci.getIntName()));
  }

  magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream& XML, 
				     const Interaction& g)
  {
    g.outputXML(XML);
    return XML;
  }

  shared_ptr<C2Range>& 
  Interaction::getRange() 
  { return range; }

  const shared_ptr<C2Range>& 
  Interaction::getRange() const
  { return range; }

  shared_ptr<Interaction>
  Interaction::getClass(const magnet::xml::Node& XML, dynamo::Simulation* Sim)
  {
    if (!std::strcmp(XML.getAttribute("Type"),"HardSphere"))
      return shared_ptr<Interaction>(new IHardSphere(XML, Sim));
    else if (!std::strcmp(XML.getAttribute("Type"),"RoughHardSphere"))
      return shared_ptr<Interaction>(new IRoughHardSphere(XML, Sim));
    else if (!std::strcmp(XML.getAttribute("Type"),"SquareWell"))
      return shared_ptr<Interaction>(new ISquareWell(XML, Sim));
    else if (!std::strcmp(XML.getAttribute("Type"),"PRIME_BB"))
      return shared_ptr<Interaction>(new IPRIME_BB(XML, Sim));
    else if (!std::strcmp(XML.getAttribute("Type"),"ThinThread"))
      return shared_ptr<Interaction>(new IThinThread(XML, Sim));
    else if (!std::strcmp(XML.getAttribute("Type"),"SquareWellSeq"))
      return shared_ptr<Interaction>(new ISWSequence(XML, Sim));
    else if (!std::strcmp(XML.getAttribute("Type"),"SquareBond"))
      return shared_ptr<Interaction>(new ISquareBond(XML, Sim));
    else if (!std::strcmp(XML.getAttribute("Type"),"SoftCore"))
      return shared_ptr<Interaction>(new ISoftCore(XML, Sim));
    else if (!std::strcmp(XML.getAttribute("Type"),"Null"))
      return shared_ptr<Interaction>(new INull(XML, Sim));
    else if (!std::strcmp(XML.getAttribute("Type"),"Lines"))
      return shared_ptr<Interaction>(new ILines(XML, Sim));
    else if (!std::strcmp(XML.getAttribute("Type"),"Dumbbells"))
      return shared_ptr<Interaction>(new IDumbbells(XML, Sim));
    else if (!std::strcmp(XML.getAttribute("Type"),"ParallelCubes"))
      return shared_ptr<Interaction>(new IParallelCubes(XML, Sim));
    else if (!std::strcmp(XML.getAttribute("Type"),"Stepped"))
      return shared_ptr<Interaction>(new IStepped(XML, Sim));
    else 
      M_throw() << XML.getAttribute("Type")
		<< ", Unknown type of interaction encountered";
  }
}
