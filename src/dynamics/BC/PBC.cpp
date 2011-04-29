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

#include "PBC.hpp"
#include "../../base/is_simdata.hpp"
#include "../../extcode/mathtemplates.hpp"
#include <magnet/xmlwriter.hpp>

BCPeriodic::BCPeriodic(const DYNAMO::SimData* tmp):
  BoundaryCondition(tmp, "RPBC", IC_purple)
{
  Sim = tmp;
}

void 
BCPeriodic::applyBC(Vector & pos) const
{ 
  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= Sim->primaryCellSize[n] *
      rintfunc (pos[n]/Sim->primaryCellSize[n]);    
}

void 
BCPeriodic::applyBC(Vector & pos, Vector&) const
{
  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= Sim->primaryCellSize[n] *
      rintfunc (pos[n] / Sim->primaryCellSize[n]);    
}

void 
BCPeriodic::applyBC(Vector  &pos, const double&) const 
{
  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= Sim->primaryCellSize[n] *
      rintfunc (pos[n] / Sim->primaryCellSize[n]);    
}

void 
BCPeriodic::outputXML(xml::XmlStream &XML) const
{
  XML << xml::attr("Boundary") << "PBC";
}

void 
BCPeriodic::operator<<(const magnet::xml::Node&) 
{}

BoundaryCondition* 
BCPeriodic::Clone () const 
{ return new BCPeriodic(*this); }


BCPeriodicExceptX::BCPeriodicExceptX(const DYNAMO::SimData* tmp):
  BoundaryCondition(tmp, "NoXPBC",IC_purple)
{ Sim = tmp; }

void 
BCPeriodicExceptX::outputXML(xml::XmlStream &XML) const
{
  XML << xml::attr("Boundary") << "NoXPBC";
}

void 
BCPeriodicExceptX::operator<<(const magnet::xml::Node&) 
{}

BoundaryCondition* 
BCPeriodicExceptX::Clone () const 
{ return new BCPeriodicExceptX(*this); }

void 
BCPeriodicExceptX::applyBC(Vector & pos) const
{ 
  double x = pos[0];

  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= Sim->primaryCellSize[n] *
      rintfunc (pos[n] / Sim->primaryCellSize[n]);    

  pos[0] = x;
}
  
void 
BCPeriodicExceptX::applyBC(Vector & pos, Vector&) const
{ 
  double x = pos[0];

  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= Sim->primaryCellSize[n] *
      rintfunc (pos[n]/Sim->primaryCellSize[n]);    

  pos[0] = x;
}

void 
BCPeriodicExceptX::applyBC(Vector  &pos, const double&) const 
{ 
  double x = pos[0];

  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= Sim->primaryCellSize[n] *
      rintfunc (pos[n] / Sim->primaryCellSize[n]);    
  
  applyBC(pos); 

  pos[0] = x;
}

BCPeriodicXOnly::BCPeriodicXOnly(const DYNAMO::SimData* tmp):
  BoundaryCondition(tmp, "NoXPBC",IC_purple)
{
  Sim = tmp;
}

void 
BCPeriodicXOnly::outputXML(xml::XmlStream &XML) const
{
  XML << xml::attr("Boundary") << "OnlyXPBC";
}

void 
BCPeriodicXOnly::operator<<(const magnet::xml::Node&) 
{}

BoundaryCondition* 
BCPeriodicXOnly::Clone () const 
{ return new BCPeriodicXOnly(*this); }

void 
BCPeriodicXOnly::applyBC(Vector & pos) const
{ 
  pos[0] -= Sim->primaryCellSize[0] 
    * rintfunc (pos[0] / Sim->primaryCellSize[0]);
}
  
void 
BCPeriodicXOnly::applyBC(Vector & pos, Vector&) const
{ 
  applyBC(pos);
}

void 
BCPeriodicXOnly::applyBC(Vector  &pos, const double&) const 
{ 
  applyBC(pos); 
}
