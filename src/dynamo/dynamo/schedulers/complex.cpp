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

#include <dynamo/schedulers/complex.hpp>
#include <dynamo/dynamics/interactions/intEvent.hpp>
#include <dynamo/simulation/particle.hpp>

#include <dynamo/dynamics/liouvillean/liouvillean.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/dynamics/systems/system.hpp>
#include <dynamo/dynamics/globals/global.hpp>
#include <dynamo/dynamics/globals/globEvent.hpp>
#include <dynamo/dynamics/globals/neighbourList.hpp>
#include <dynamo/dynamics/locals/local.hpp>
#include <dynamo/dynamics/locals/localEvent.hpp>
#include <dynamo/schedulers/complexentries/include.hpp>

#include <magnet/xmlreader.hpp>
#include <boost/bind.hpp>
#include <boost/progress.hpp>
#include <cmath> //for huge val

namespace dynamo {
  void 
  SComplex::operator<<(const magnet::xml::Node& XML)
  {
    Scheduler::operator<<(XML);

    for (magnet::xml::Node node = XML.getNode("Entries").fastGetNode("Entry"); 
	 node.valid(); ++node)
      entries.push_back(shared_ptr<SCEntry>(SCEntry::getClass(node, Sim)));
  }

  void
  SComplex::initialise()
  {
    dout << "Reinitialising on collision " << Sim->eventCount << std::endl;
    std::cout.flush();

    BOOST_FOREACH(shared_ptr<SCEntry>& ent, entries)
      ent->initialise();

    Scheduler::initialise();
  }

  void 
  SComplex::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "Complex"
	<< magnet::xml::tag("Sorter")
	<< *sorter
	<< magnet::xml::endtag("Sorter")
	<< magnet::xml::tag("Entries");
  
    BOOST_FOREACH(const shared_ptr<SCEntry>& ent,  entries)
      XML << magnet::xml::tag("Entry")
	  << ent
	  << magnet::xml::endtag("Entry");

    XML << magnet::xml::endtag("Entries");
  }

  SComplex::SComplex(const magnet::xml::Node& XML, dynamo::SimData* const Sim):
    Scheduler(Sim,"ComplexScheduler", NULL)
  { 
    dout << "Complex Scheduler Algorithmn Loaded" << std::endl;
    operator<<(XML);
  }

  SComplex::SComplex(dynamo::SimData* const Sim, FEL* ns):
    Scheduler(Sim,"ComplexScheduler", ns)
  { dout << "Complex Scheduler Algorithmn Loaded" << std::endl; }

  void 
  SComplex::getParticleNeighbourhood(const Particle& part,
				     const nbHoodFunc& func) const
  {
    BOOST_FOREACH(const shared_ptr<SCEntry>& ent, entries)
      if (ent->isApplicable(part))
	ent->getParticleNeighbourhood(part, func);
  }

  void 
  SComplex::getParticleNeighbourhood(const Vector& vec,
				     const nbHoodFunc2& func) const
  {
    BOOST_FOREACH(const shared_ptr<SCEntry>& ent, entries)
	ent->getParticleNeighbourhood(vec, func);
  }
    
  void 
  SComplex::getLocalNeighbourhood(const Particle& part, 
				  const nbHoodFunc& func) const
  {
    BOOST_FOREACH(const shared_ptr<SCEntry>& ent, entries)
      if (ent->isApplicable(part))
	ent->getLocalNeighbourhood(part, func);
  }
}
