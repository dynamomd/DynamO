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

#include "SW.hpp"
#include <boost/lexical_cast.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../../base/is_exception.hpp"
#include <cmath>
#include <cstring>

USquareWell::USquareWell(const DYNAMO::SimData* tmp):
  Units(tmp),
  UnitOfEnergy(0.0),
  UnitOfLength(1.0)
{
  I_cout() << "SW units loaded";
}

USquareWell::USquareWell(Iflt diameter, Iflt energy, const DYNAMO::SimData* tmp):
  Units(tmp),
  UnitOfEnergy(energy),
  UnitOfLength(diameter)
{
  I_cout() << "SW units loaded";
}

USquareWell::USquareWell(const XMLNode &XML, const DYNAMO::SimData* tmp):
  Units(tmp)
{ 
  operator<<(XML); 
  I_cout() << "SW units loaded";
}

USquareWell::~USquareWell() {}

Iflt 
USquareWell::unitLength() const
{ return UnitOfLength; }

void 
USquareWell::setUnitLength(Iflt scalar)
{ UnitOfLength = scalar; }

void 
USquareWell::rescaleLength(Iflt rs)
{ UnitOfLength += rs * UnitOfLength; }

Iflt 
USquareWell::unitTime() const
{
  return sqrt(unitLength()*unitLength()*unitMass()/UnitOfEnergy);
}

Units* 
USquareWell::Clone() const
{ return new USquareWell(*this); }
  
void 
USquareWell::operator<<(const XMLNode &XML)
{
  if (std::strcmp(XML.getAttribute("Type"),"SW"))
    D_throw() << "Attempting to load USquareWell from non elastic type";
  
  try {
    UnitOfLength = 1.0/(boost::lexical_cast<Iflt>(XML.getAttribute("BoxLength")));
    UnitOfEnergy = boost::lexical_cast<Iflt>(XML.getAttribute("Energy"));
  }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in USquareWell";
    }
}

void 
USquareWell::outputXML(xmlw::XmlStream &XML) const
{
  XML << xmlw::attr("Type") << "SW"
      << xmlw::attr("BoxLength") << 1.0/UnitOfLength 
      << xmlw::attr("Energy") << UnitOfEnergy; 
}


