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

#include "elastic.hpp"
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

UHardSphere::UHardSphere(const DYNAMO::SimData* tmp): 
  Units(tmp),
  UnitOfLength(1.0)
{
  I_cout() << "HardSphere units loaded";
}
  
UHardSphere::UHardSphere(double diameter, const DYNAMO::SimData* tmp):
  Units(tmp),
  UnitOfLength(diameter)
{
  I_cout() << "HardSphere units loaded";
}

UHardSphere::UHardSphere(const magnet::xml::Node& XML, const DYNAMO::SimData* tmp):
  Units(tmp)
{ 
  operator<<(XML); 
  I_cout() << "HardSphere sphere units loaded";
}

UHardSphere::~UHardSphere() {}

double
UHardSphere::unitLength() const
{ return UnitOfLength; }

void 
UHardSphere::setUnitLength(double scalar)
{ UnitOfLength = scalar; }

double 
UHardSphere::unitTime() const
{ return 1.0; }
  
Units* 
UHardSphere::Clone() const
{ return new UHardSphere(*this); }

void 
UHardSphere::rescaleLength(double rs)
{ UnitOfLength += rs * UnitOfLength; }
  
void 
UHardSphere::operator<<(const magnet::xml::Node& XML)
{  
  try {
    UnitOfLength = 1.0 / XML.getAttribute("BoxLength").as<double>();
  }
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in UHardSphere";
    }
}

void 
UHardSphere::outputXML(xml::XmlStream &XML) const
{
  XML << xml::attr("Type") << "HardSphere"
      << xml::attr("BoxLength") << 1.0/UnitOfLength; 
}

