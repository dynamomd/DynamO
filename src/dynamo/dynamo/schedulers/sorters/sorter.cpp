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

#include <dynamo/schedulers/sorters/include.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  EventSorter::EventSorter(const dynamo::SimData* const& SD, const char *aName):
    SimBase_const(SD, aName)
  {}

  shared_ptr<EventSorter>
  EventSorter::getClass(const magnet::xml::Node& XML, const dynamo::SimData* Sim)
  {
    if (std::string(XML.getAttribute("Type")) == ESBoundedPQName<pList>::name())
      return shared_ptr<EventSorter>(new ESBoundedPQ<>(Sim));
    if (std::string(XML.getAttribute("Type")) == ESBoundedPQName<PELSingleEvent>::name())
      return shared_ptr<EventSorter>(new ESBoundedPQ<PELSingleEvent>(Sim));
    if (std::string(XML.getAttribute("Type")) == ESBoundedPQName<MinMaxHeapPList<2> >::name())
      return shared_ptr<EventSorter>(new ESBoundedPQ<MinMaxHeapPList<2> >(Sim));
    if (std::string(XML.getAttribute("Type")) == ESBoundedPQName<MinMaxHeapPList<3> >::name())
      return shared_ptr<EventSorter>(new ESBoundedPQ<MinMaxHeapPList<3> >(Sim));
    if (std::string(XML.getAttribute("Type")) == ESBoundedPQName<MinMaxHeapPList<4> >::name())
      return shared_ptr<EventSorter>(new ESBoundedPQ<MinMaxHeapPList<4> >(Sim));
    if (std::string(XML.getAttribute("Type")) == ESBoundedPQName<MinMaxHeapPList<5> >::name())
      return shared_ptr<EventSorter>(new ESBoundedPQ<MinMaxHeapPList<5> >(Sim));
    if (std::string(XML.getAttribute("Type")) == ESBoundedPQName<MinMaxHeapPList<6> >::name())
      return shared_ptr<EventSorter>(new ESBoundedPQ<MinMaxHeapPList<6> >(Sim));
    if (std::string(XML.getAttribute("Type")) == ESBoundedPQName<MinMaxHeapPList<7> >::name())
      return shared_ptr<EventSorter>(new ESBoundedPQ<MinMaxHeapPList<7> >(Sim));
    if (std::string(XML.getAttribute("Type")) == ESBoundedPQName<MinMaxHeapPList<8> >::name())
      return shared_ptr<EventSorter>(new ESBoundedPQ<MinMaxHeapPList<8> >(Sim));
    else if (std::string(XML.getAttribute("Type")) == std::string("CBT"))
      return shared_ptr<EventSorter>(new ESCBT(Sim));
    else 
      M_throw() << "Unknown type of Sorter encountered";
  }

  magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream& XML, const EventSorter& srtr)
  {
    srtr.outputXML(XML);
    return XML;
  }
}
