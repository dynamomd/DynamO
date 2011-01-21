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

#include "None.hpp"
#include "../../extcode/xmlwriter.hpp"

BCNone::BCNone(const DYNAMO::SimData* Sim):
  BoundaryCondition(Sim, "NullBC", IC_purple)
{ I_cout() << "No boundary condition loaded"; }

BCNone::~BCNone() {}
    
void 
BCNone::applyBC(Vector  &)const 
{}

void 
BCNone::applyBC(Vector  &, Vector &) const 
{}

void 
BCNone::applyBC(Vector  &, const double&) const 
{}

void 
BCNone::update(const double &) 
{}

void 
BCNone::outputXML(xml::XmlStream &XML) const
{
  XML << xml::attr("Boundary") << "None";
}

void 
BCNone::operator<<(const XMLNode &)
{}

BoundaryCondition* 
BCNone::Clone () const 
{ return new BCNone(*this); }

void 
BCNone::rounding(Vector &) const 
{}
