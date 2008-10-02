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

#include "Null.hpp"
#include "../../extcode/xmlwriter.hpp"

CNullBC::CNullBC(const DYNAMO::SimData*):
  Base_Class("NullBC", IC_purple)
{ I_cout() << "No boundary condition loaded"; }

CNullBC::~CNullBC() {}
    
void 
CNullBC::setPBC(CVector<> &)const 
{}

void 
CNullBC::setPBC(CVector<> &, CVector<> &) const 
{}

void 
CNullBC::setPBC(CVector<> &, Iflt) const 
{}

void 
CNullBC::update(const Iflt &) 
{}

void 
CNullBC::outputXML(xmlw::XmlStream &XML) const
{
  XML << xmlw::attr("Boundary") << "Null";
}

void 
CNullBC::operator<<(const XMLNode &)
{}

CBC* 
CNullBC::Clone () const 
{ return new CNullBC(*this); }

void 
CNullBC::rounding(CVector<>&) const 
{}
