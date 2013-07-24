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
  Interaction::Interaction(dynamo::Simulation* tmp, IDPairRange* nR):
    SimBase(tmp, "Interaction"),
    range(nR)
  {}

  void 
  Interaction::operator<<(const magnet::xml::Node& XML)
  { range = shared_ptr<IDPairRange>(IDPairRange::getClass(XML.getNode("IDPairRange"), Sim)); }

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

  shared_ptr<IDPairRange>& 
  Interaction::getRange() 
  { return range; }

  const shared_ptr<IDPairRange>& 
  Interaction::getRange() const
  { return range; }

  shared_ptr<Interaction>
  Interaction::getClass(const magnet::xml::Node& XML, dynamo::Simulation* Sim)
  {
    if (!XML.getAttribute("Type").getValue().compare("HardSphere"))
      return shared_ptr<Interaction>(new IHardSphere(XML, Sim));
    else if (!XML.getAttribute("Type").getValue().compare("SquareWell"))
      return shared_ptr<Interaction>(new ISquareWell(XML, Sim));
    else if (!XML.getAttribute("Type").getValue().compare("PRIME_BB"))
      return shared_ptr<Interaction>(new IPRIME_BB(XML, Sim));
    else if (!XML.getAttribute("Type").getValue().compare("ThinThread"))
      return shared_ptr<Interaction>(new IThinThread(XML, Sim));
    else if (!XML.getAttribute("Type").getValue().compare("SquareWellSeq"))
      return shared_ptr<Interaction>(new ISWSequence(XML, Sim));
    else if (!XML.getAttribute("Type").getValue().compare("SquareBond"))
      return shared_ptr<Interaction>(new ISquareBond(XML, Sim));
    else if (!XML.getAttribute("Type").getValue().compare("Null"))
      return shared_ptr<Interaction>(new INull(XML, Sim));
    else if (!XML.getAttribute("Type").getValue().compare("Lines"))
      return shared_ptr<Interaction>(new ILines(XML, Sim));
    else if (!XML.getAttribute("Type").getValue().compare("Dumbbells"))
      return shared_ptr<Interaction>(new IDumbbells(XML, Sim));
    else if (!XML.getAttribute("Type").getValue().compare("ParallelCubes"))
      return shared_ptr<Interaction>(new IParallelCubes(XML, Sim));
    else if (!XML.getAttribute("Type").getValue().compare("Stepped"))
      return shared_ptr<Interaction>(new IStepped(XML, Sim));
    else 
      M_throw() << XML.getAttribute("Type").getValue()
		<< ", Unknown type of interaction encountered";
  }
}
