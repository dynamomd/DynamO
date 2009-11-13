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

#include "PBC.hpp"
#include "../../base/is_simdata.hpp"
#include "rintfunc.hpp"
#include "../../extcode/xmlwriter.hpp"

BCSquarePeriodic::BCSquarePeriodic(const DYNAMO::SimData* tmp):
  BoundaryCondition(tmp, "SPBC", IC_purple)
{
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    if (Sim->aspectRatio[iDim] != 1.0)
      D_throw() << "The simulation aspect ratio is not unity for the use of "
	"square PBC's";
  Sim = tmp;
}
  
void 
BCSquarePeriodic::outputXML(xmlw::XmlStream &XML) const
{
  XML << xmlw::attr("Shape") << "Square"
      << xmlw::attr("Boundary") << "PBC";
}

BoundaryCondition* 
BCSquarePeriodic::Clone () const 
{ return new BCSquarePeriodic(*this); }


void 
BCSquarePeriodic::operator<<(const XMLNode&) 
{}

void 
BCSquarePeriodic::applyBC(Vector & pos) const
{ 
  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= rintfunc (pos[n]);
}

void 
BCSquarePeriodic::applyBC(Vector & pos, Vector&) const
{ 
  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= rintfunc (pos[n]);
}

void 
BCSquarePeriodic::applyBC(Vector  &pos, const Iflt&) const 
{ 
  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= rintfunc (pos[n]);
}

BCRectangularPeriodic::BCRectangularPeriodic(const DYNAMO::SimData* tmp):
  BoundaryCondition(tmp, "RPBC", IC_purple)
{
  Sim = tmp;
}

void 
BCRectangularPeriodic::applyBC(Vector & pos) const
{ 
  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= Sim->aspectRatio[n] *
      rintfunc (pos[n]/Sim->aspectRatio[n]);    
}

void 
BCRectangularPeriodic::applyBC(Vector & pos, Vector&) const
{
  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= Sim->aspectRatio[n] *
      rintfunc (pos[n] / Sim->aspectRatio[n]);    
}

void 
BCRectangularPeriodic::applyBC(Vector  &pos, const Iflt&) const 
{
  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= Sim->aspectRatio[n] *
      rintfunc (pos[n] / Sim->aspectRatio[n]);    
}

void 
BCRectangularPeriodic::outputXML(xmlw::XmlStream &XML) const
{
  XML << xmlw::attr("Shape") << "Rectangular"
      << xmlw::attr("Boundary") << "PBC";
}

void 
BCRectangularPeriodic::operator<<(const XMLNode&) 
{}

BoundaryCondition* 
BCRectangularPeriodic::Clone () const 
{ return new BCRectangularPeriodic(*this); }


BCSquarePeriodicExceptX::BCSquarePeriodicExceptX(const DYNAMO::SimData* tmp):
  BoundaryCondition(tmp, "RNoXPBC",IC_purple)
{
  Sim = tmp;
}

void 
BCSquarePeriodicExceptX::outputXML(xmlw::XmlStream &XML) const
{
  XML << xmlw::attr("Shape") << "Rectangular"
      << xmlw::attr("Boundary") << "NoXPBC";
}

void 
BCSquarePeriodicExceptX::operator<<(const XMLNode&) 
{}

BoundaryCondition* 
BCSquarePeriodicExceptX::Clone () const 
{ return new BCSquarePeriodicExceptX(*this); }

void 
BCSquarePeriodicExceptX::applyBC(Vector & pos) const
{ 
  Iflt x = pos[0];

  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= Sim->aspectRatio[n] *
      rintfunc (pos[n] / Sim->aspectRatio[n]);    

  pos[0] = x;
}
  
void 
BCSquarePeriodicExceptX::applyBC(Vector & pos, Vector&) const
{ 
  Iflt x = pos[0];

  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= Sim->aspectRatio[n] *
      rintfunc (pos[n]/Sim->aspectRatio[n]);    

  pos[0] = x;
}

void 
BCSquarePeriodicExceptX::applyBC(Vector  &pos, const Iflt&) const 
{ 
  Iflt x = pos[0];

  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= Sim->aspectRatio[n] *
      rintfunc (pos[n] / Sim->aspectRatio[n]);    
  
  applyBC(pos); 

  pos[0] = x;
}

BCSquarePeriodicXOnly::BCSquarePeriodicXOnly(const DYNAMO::SimData* tmp):
  BoundaryCondition(tmp, "RNoXPBC",IC_purple)
{
  Sim = tmp;
}

void 
BCSquarePeriodicXOnly::outputXML(xmlw::XmlStream &XML) const
{
  XML << xmlw::attr("Shape") << "Rectangular"
      << xmlw::attr("Boundary") << "OnlyXPBC";
}

void 
BCSquarePeriodicXOnly::operator<<(const XMLNode&) 
{}

BoundaryCondition* 
BCSquarePeriodicXOnly::Clone () const 
{ return new BCSquarePeriodicXOnly(*this); }

void 
BCSquarePeriodicXOnly::applyBC(Vector & pos) const
{ 
  pos[0] -= Sim->aspectRatio[0] 
    * rintfunc (pos[0] / Sim->aspectRatio[0]);
}
  
void 
BCSquarePeriodicXOnly::applyBC(Vector & pos, Vector&) const
{ 
  applyBC(pos);
}

void 
BCSquarePeriodicXOnly::applyBC(Vector  &pos, const Iflt&) const 
{ 
  applyBC(pos); 
}
