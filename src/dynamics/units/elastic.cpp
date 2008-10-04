/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "elastic.hpp"
#include <boost/lexical_cast.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../../base/is_exception.hpp"

CUElastic::CUElastic(const DYNAMO::SimData* tmp): 
  CUnits(tmp),
  UnitOfLength(1.0)
{
  I_cout() << "Elastic units loaded";
}
  
CUElastic::CUElastic(Iflt diameter, const DYNAMO::SimData* tmp):
  CUnits(tmp),
  UnitOfLength(diameter)
{
  I_cout() << "Elastic units loaded";
}

CUElastic::CUElastic(const XMLNode &XML, const DYNAMO::SimData* tmp):
  CUnits(tmp)
{ 
  operator<<(XML); 
  I_cout() << "Elastic units loaded";
}

CUElastic::~CUElastic() {}

Iflt
CUElastic::unitLength() const
{ return UnitOfLength; }

void 
CUElastic::setUnitLength(Iflt scalar)
{ UnitOfLength = scalar; }

Iflt 
CUElastic::unitTime() const
{ return 1.0; }
  
CUnits* 
CUElastic::Clone() const
{ return new CUElastic(*this); }

void 
CUElastic::rescaleLength(Iflt rs)
{ UnitOfLength += rs * UnitOfLength; }
  
void 
CUElastic::operator<<(const XMLNode &XML)
{
  if (strcmp(XML.getAttribute("Type"),"Elastic"))
    D_throw() << "Attempting to load CUElastic from non elastic type";
  
  try {
    UnitOfLength = 1.0/(boost::lexical_cast<Iflt>(XML.getAttribute("BoxLength")));
  }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in CUElastic";
    }
}

void 
CUElastic::outputXML(xmlw::XmlStream &XML) const
{
  XML << xmlw::attr("Type") << "Elastic"
      << xmlw::attr("BoxLength") << 1.0/UnitOfLength; 
}

