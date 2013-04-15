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

#include <dynamo/schedulers/sorters/include.hpp>
#include <dynamo/simulation.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  shared_ptr<FEL>
  FEL::getClass(const magnet::xml::Node& XML)
  {
    if (std::string(XML.getAttribute("Type")) == FELBoundedPQName<PELHeap>::name())
      return shared_ptr<FEL>(new FELBoundedPQ<>());
    if (std::string(XML.getAttribute("Type")) == FELBoundedPQName<PELSingleEvent>::name())
      return shared_ptr<FEL>(new FELBoundedPQ<PELSingleEvent>());
    if (std::string(XML.getAttribute("Type")) == FELBoundedPQName<PELMinMax<2> >::name())
      return shared_ptr<FEL>(new FELBoundedPQ<PELMinMax<2> >());
    if (std::string(XML.getAttribute("Type")) == FELBoundedPQName<PELMinMax<3> >::name())
      return shared_ptr<FEL>(new FELBoundedPQ<PELMinMax<3> >());
    if (std::string(XML.getAttribute("Type")) == FELBoundedPQName<PELMinMax<4> >::name())
      return shared_ptr<FEL>(new FELBoundedPQ<PELMinMax<4> >());
    if (std::string(XML.getAttribute("Type")) == FELBoundedPQName<PELMinMax<5> >::name())
      return shared_ptr<FEL>(new FELBoundedPQ<PELMinMax<5> >());
    if (std::string(XML.getAttribute("Type")) == FELBoundedPQName<PELMinMax<6> >::name())
      return shared_ptr<FEL>(new FELBoundedPQ<PELMinMax<6> >());
    if (std::string(XML.getAttribute("Type")) == FELBoundedPQName<PELMinMax<7> >::name())
      return shared_ptr<FEL>(new FELBoundedPQ<PELMinMax<7> >());
    if (std::string(XML.getAttribute("Type")) == FELBoundedPQName<PELMinMax<8> >::name())
      return shared_ptr<FEL>(new FELBoundedPQ<PELMinMax<8> >());
    else if (std::string(XML.getAttribute("Type")) == std::string("CBT"))
      return shared_ptr<FEL>(new FELCBT());
    else 
      M_throw() << "Unknown type of Sorter encountered";
  }

  magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream& XML, const FEL& srtr)
  {
    srtr.outputXML(XML);
    return XML;
  }
}
