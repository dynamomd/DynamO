/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "shear.hpp"
#include "../../extcode/xmlParser.h"
#include "../../extcode/xmlwriter.hpp"
#include "../../base/is_exception.hpp"
#include <cstring>

void 
CUShear::operator<<(const XMLNode& XML)
{
  if (std::strcmp(XML.getAttribute("Type"),"Shear"))
    D_throw() << "Attempting to load CUShear from non shear type";
  
  try {
    UnitOfLength = 1.0 / boost::lexical_cast<Iflt>
      (XML.getAttribute("BoxLength"));
  }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in CUElastic";
    }
}

void 
CUShear::outputXML(xmlw::XmlStream &XML) const
{
  XML << xmlw::attr("Type") << "Shear"
      << xmlw::attr("BoxLength") << 1.0 / UnitOfLength; 
}
