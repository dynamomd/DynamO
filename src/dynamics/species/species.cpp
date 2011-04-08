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
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

Species* 
Species::getClass(const magnet::xml::Node& XML, DYNAMO::SimData* tmp, size_t nID)
{
  if (!std::strcmp(XML.getAttribute("Type"), "Point"))
    return new SpPoint(XML, tmp, nID);
  else if (!std::strcmp(XML.getAttribute("Type"), "SphericalTop"))
    return new SpSphericalTop(XML, tmp, nID);
  else if (!std::strcmp(XML.getAttribute("Type"), "Lines"))
    return new SpLines(XML, tmp, nID);
  else if (!std::strcmp(XML.getAttribute("Type"), "Dumbbells"))
    return new SpDumbbells(XML, tmp, nID);
  else if (!std::strcmp(XML.getAttribute("Type"), "FixedCollider"))
    return new SpFixedCollider(XML, tmp, nID);
  else 
    M_throw() << XML.getAttribute("Type")
	      << ", Unknown type of species encountered";
}
