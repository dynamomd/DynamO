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

#include "include.hpp"
#include "../../base/is_simdata.hpp"
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

CSSorter::CSSorter(const DYNAMO::SimData* const& SD, const char *aName):
  SimBase_const(SD, aName, IC_white_brown)
{}

CSSorter* 
CSSorter::getClass(const magnet::xml::Node& XML, const DYNAMO::SimData* Sim)
{
  if (std::string(XML.getAttribute("Type")) == CSSBoundedPQName<pList>::name())
    return new CSSBoundedPQ<>(Sim);
  if (std::string(XML.getAttribute("Type")) == CSSBoundedPQName<MinMaxHeapPList<2> >::name())
    return new CSSBoundedPQ<MinMaxHeapPList<2> >(Sim);
  if (std::string(XML.getAttribute("Type")) == CSSBoundedPQName<MinMaxHeapPList<3> >::name())
    return new CSSBoundedPQ<MinMaxHeapPList<3> >(Sim);
  if (std::string(XML.getAttribute("Type")) == CSSBoundedPQName<MinMaxHeapPList<4> >::name())
    return new CSSBoundedPQ<MinMaxHeapPList<4> >(Sim);
  if (std::string(XML.getAttribute("Type")) == CSSBoundedPQName<MinMaxHeapPList<5> >::name())
    return new CSSBoundedPQ<MinMaxHeapPList<5> >(Sim);
  if (std::string(XML.getAttribute("Type")) == CSSBoundedPQName<MinMaxHeapPList<6> >::name())
    return new CSSBoundedPQ<MinMaxHeapPList<6> >(Sim);
  if (std::string(XML.getAttribute("Type")) == CSSBoundedPQName<MinMaxHeapPList<7> >::name())
    return new CSSBoundedPQ<MinMaxHeapPList<7> >(Sim);
  if (std::string(XML.getAttribute("Type")) == CSSBoundedPQName<MinMaxHeapPList<8> >::name())
    return new CSSBoundedPQ<MinMaxHeapPList<8> >(Sim);
  else if (std::string(XML.getAttribute("Type")) == std::string("CBT"))
    return new CSSCBT(Sim);
  else 
    M_throw() << "Unknown type of Sorter encountered";
}

xml::XmlStream& operator<<(xml::XmlStream& XML, const CSSorter& srtr)
{
  srtr.outputXML(XML);
  return XML;
}
