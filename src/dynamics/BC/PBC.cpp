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

#include "PBC.hpp"
#include "shapes.hpp"
#include "../../extcode/xmlwriter.hpp"

CSPBC::CSPBC(DYNAMO::SimData* tmp):
  Base_Class("SPBC",IC_purple)
{
  Sim = tmp;
}
  
void 
CSPBC::outputXML(xmlw::XmlStream &XML) const
{
  XML << xmlw::attr("Shape") << "Square"
      << xmlw::attr("Boundary") << "PBC";
}

CBC* 
CSPBC::Clone () const 
{ return new CSPBC(*this); }


void 
CSPBC::operator<<(const XMLNode&) 
{}

CRPBC::CRPBC(DYNAMO::SimData* tmp):
  Base_Class("RPBC",IC_purple)
{
  Sim = tmp;
}

void 
CRPBC::outputXML(xmlw::XmlStream &XML) const
{
  XML << xmlw::attr("Shape") << "Rectangular"
      << xmlw::attr("Boundary") << "PBC";
}

void 
CRPBC::operator<<(const XMLNode&) 
{}

CBC* 
CRPBC::Clone () const 
{ return new CRPBC(*this); }


CRNoXPBC::CRNoXPBC(DYNAMO::SimData* tmp):
  Base_Class("RNoXPBC",IC_purple)
{
  Sim = tmp;
}

void 
CRNoXPBC::outputXML(xmlw::XmlStream &XML) const
{
  XML << xmlw::attr("Shape") << "Rectangular"
      << xmlw::attr("Boundary") << "NoXPBC";
}

void 
CRNoXPBC::operator<<(const XMLNode&) 
{}

CBC* 
CRNoXPBC::Clone () const 
{ return new CRNoXPBC(*this); }

