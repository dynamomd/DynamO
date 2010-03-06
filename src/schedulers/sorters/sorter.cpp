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

#include "include.hpp"
#include "../../base/is_simdata.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"

CSSorter::CSSorter(const DYNAMO::SimData* const& SD, const char *aName):
  SimBase_const(SD, aName, IC_white_brown)
{

}

CSSorter* 
CSSorter::getClass(const XMLNode& XML, const DYNAMO::SimData* Sim)
{
  if (!strcmp(XML.getAttribute("Type"),"BoundedPQ"))
    return new CSSBoundedPQ<>(Sim);
  if (!strcmp(XML.getAttribute("Type"),"BoundedPQMinMax"))
    {
      MinMaxHeapPList::HeapSize 
	= boost::lexical_cast<size_t>(XML.getAttribute("HeapSize"));

      return new CSSBoundedPQ<MinMaxHeapPList>(Sim);
    }
  else if (!strcmp(XML.getAttribute("Type"),"CBT"))
    return new CSSCBT(Sim);
  else 
    D_throw() << "Unknown type of Sorter encountered";
}

xmlw::XmlStream& operator<<(xmlw::XmlStream& XML, const CSSorter& srtr)
{
  srtr.outputXML(XML);
  return XML;
}
