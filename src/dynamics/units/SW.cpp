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

#include "SW.hpp"
#include <boost/lexical_cast.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../../base/is_exception.hpp"
#include <cmath>

CUSW::CUSW(const DYNAMO::SimData* tmp):
  CUnits(tmp),
  UnitOfEnergy(0.0),
  UnitOfLength(1.0)
{
  I_cout() << "SW units loaded";
}

CUSW::CUSW(Iflt diameter, Iflt energy, const DYNAMO::SimData* tmp):
  CUnits(tmp),
  UnitOfEnergy(energy),
  UnitOfLength(diameter)
{
  I_cout() << "SW units loaded";
}

CUSW::CUSW(const XMLNode &XML, const DYNAMO::SimData* tmp):
  CUnits(tmp)
{ 
  operator<<(XML); 
  I_cout() << "SW units loaded";
}

CUSW::~CUSW() {}

Iflt 
CUSW::unitLength() const
{ return UnitOfLength; }

void 
CUSW::setUnitLength(Iflt scalar)
{ UnitOfLength = scalar; }

void 
CUSW::rescaleLength(Iflt rs)
{ UnitOfLength += rs * UnitOfLength; }

Iflt 
CUSW::unitTime() const
{
  return sqrt(unitLength()*unitLength()*unitMass()/UnitOfEnergy);
}

CUnits* 
CUSW::Clone() const
{ return new CUSW(*this); }
  
void 
CUSW::operator<<(const XMLNode &XML)
{
  if (strcmp(XML.getAttribute("Type"),"SW"))
    I_throw() << "Attempting to load CUSW from non elastic type";
  
  try {
    UnitOfLength = 1.0/(boost::lexical_cast<Iflt>(XML.getAttribute("BoxLength")));
    UnitOfEnergy = boost::lexical_cast<Iflt>(XML.getAttribute("Energy"));
  }
  catch (boost::bad_lexical_cast &)
    {
      I_throw() << "Failed a lexical cast in CUSW";
    }
}

void 
CUSW::outputXML(xmlw::XmlStream &XML) const
{
  XML << xmlw::attr("Type") << "SW"
      << xmlw::attr("BoxLength") << 1.0/UnitOfLength 
      << xmlw::attr("Energy") << UnitOfEnergy; 
}


