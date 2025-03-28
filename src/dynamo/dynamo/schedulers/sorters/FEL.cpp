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

#include <dynamo/schedulers/sorters/CBTFEL.hpp>
#include <dynamo/schedulers/sorters/MinMaxPEL.hpp>
#include <dynamo/schedulers/sorters/boundedPQFEL.hpp>
#include <dynamo/schedulers/sorters/heapPEL.hpp>
#include <dynamo/schedulers/sorters/referenceFEL.hpp>
#include <dynamo/simulation.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
shared_ptr<FEL> FEL::getClass(const magnet::xml::Node &XML) {
  if (std::string(XML.getAttribute("Type")) == std::string("BoundedPQHeap"))
    return shared_ptr<FEL>(new BoundedPQFEL<HeapPEL>());
  if (std::string(XML.getAttribute("Type")) == std::string("BoundedPQMinMax2"))
    return shared_ptr<FEL>(new BoundedPQFEL<MinMaxPEL<2>>());
  if (std::string(XML.getAttribute("Type")) == std::string("BoundedPQMinMax3"))
    return shared_ptr<FEL>(new BoundedPQFEL<MinMaxPEL<3>>());
  if (std::string(XML.getAttribute("Type")) == std::string("BoundedPQMinMax4"))
    return shared_ptr<FEL>(new BoundedPQFEL<MinMaxPEL<4>>());
  if (std::string(XML.getAttribute("Type")) == std::string("BoundedPQMinMax5"))
    return shared_ptr<FEL>(new BoundedPQFEL<MinMaxPEL<5>>());
  if (std::string(XML.getAttribute("Type")) == std::string("BoundedPQMinMax6"))
    return shared_ptr<FEL>(new BoundedPQFEL<MinMaxPEL<6>>());
  if (std::string(XML.getAttribute("Type")) == std::string("BoundedPQMinMax7"))
    return shared_ptr<FEL>(new BoundedPQFEL<MinMaxPEL<7>>());
  if (std::string(XML.getAttribute("Type")) == std::string("BoundedPQMinMax8"))
    return shared_ptr<FEL>(new BoundedPQFEL<MinMaxPEL<8>>());
  else if ((std::string(XML.getAttribute("Type")) == std::string("CBT")) ||
           (std::string(XML.getAttribute("Type")) == std::string("CBTHeap")))
    return shared_ptr<FEL>(new CBTFEL<HeapPEL>());
  else
    M_throw() << "Unknown type of Sorter encountered";
}

magnet::xml::XmlStream &operator<<(magnet::xml::XmlStream &XML,
                                   const FEL &srtr) {
  srtr.outputXML(XML);
  return XML;
}
} // namespace dynamo
