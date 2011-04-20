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

#include "SW.hpp"
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <cmath>
#include <cstring>

USquareWell::USquareWell(const DYNAMO::SimData* tmp):
  Units(tmp),
  UnitOfEnergy(0.0),
  UnitOfLength(1.0)
{
  I_cout() << "SW units loaded";
}

USquareWell::USquareWell(double diameter, double energy, const DYNAMO::SimData* tmp):
  Units(tmp),
  UnitOfEnergy(energy),
  UnitOfLength(diameter)
{
  I_cout() << "SW units loaded";
}

USquareWell::USquareWell(const magnet::xml::Node& XML, const DYNAMO::SimData* tmp):
  Units(tmp)
{ 
  operator<<(XML); 
  I_cout() << "SW units loaded";
}

USquareWell::~USquareWell() {}

double 
USquareWell::unitLength() const
{ return UnitOfLength; }

void 
USquareWell::setUnitLength(double scalar)
{ UnitOfLength = scalar; }

void 
USquareWell::rescaleLength(double factor)
{ UnitOfLength *= factor; }

double 
USquareWell::unitTime() const
{
  return sqrt(unitLength()*unitLength()*unitMass()/UnitOfEnergy);
}

Units* 
USquareWell::Clone() const
{ return new USquareWell(*this); }
  
void 
USquareWell::operator<<(const magnet::xml::Node& XML)
{
  if (std::strcmp(XML.getAttribute("Type"),"SW"))
    M_throw() << "Attempting to load USquareWell from non elastic type";
  
  try {
    UnitOfLength = 1.0/XML.getAttribute("BoxLength").as<double>();
    UnitOfEnergy = XML.getAttribute("Energy").as<double>();
  }
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in USquareWell";
    }
}

void 
USquareWell::outputXML(xml::XmlStream &XML) const
{
  XML << xml::attr("Type") << "SW"
      << xml::attr("BoxLength") << 1.0/UnitOfLength 
      << xml::attr("Energy") << UnitOfEnergy; 
}


