/*  dynamo:- Event driven molecular dynamics simulator 
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

#include <dynamo/schedulers/complexentries/rangeentry.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  SCERange::SCERange(const magnet::xml::Node& XML, dynamo::SimData* const nSim):
    SCEntry(nSim, "ParticleRange")
  { operator<<(XML); }

  void 
  SCERange::operator<<(const magnet::xml::Node& XML)
  {
    range = shared_ptr<Range>(Range::getClass(XML, Sim));    
    _testrange = shared_ptr<Range>(Range::getClass(XML.getNode("OtherParticles"), Sim));
  }

  void 
  SCERange::getParticleNeighbourhood(const Particle& part,
				     const GNeighbourList::nbHoodFunc& func) const
  {
    BOOST_FOREACH(const size_t& ID, *_testrange)
      func(part, ID);
  }

  void 
  SCERange::getParticleNeighbourhood(const Vector& vec, 
				     const GNeighbourList::nbHoodFunc2& func) const
  {
    BOOST_FOREACH(const size_t& ID, *_testrange)
      func(ID);
  }

  void 
  SCERange::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "ParticleRange"
	<< range
	<< magnet::xml::tag("OtherParticles")
	<< _testrange
	<< magnet::xml::endtag("OtherParticles")
      ;
  }
}
