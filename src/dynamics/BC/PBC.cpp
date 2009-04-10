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
#include "../../base/is_simdata.hpp"
#include "rintfunc.hpp"
#include "../../extcode/xmlwriter.hpp"

CSPBC::CSPBC(const DYNAMO::SimData* tmp):
  CBC(tmp, "SPBC", IC_purple)
{
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    if (Sim->aspectRatio[iDim] != 1.0)
      D_throw() << "The simulation aspect ratio is not unity for the use of "
	"square PBC's";
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

void 
CSPBC::setPBC(Vector & pos) const
{ 
  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= rintfunc (pos[n]);
}

void 
CSPBC::setPBC(Vector & pos, Vector&) const
{ 
  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= rintfunc (pos[n]);
}

void 
CSPBC::setPBC(Vector  &pos, const Iflt&) const 
{ 
  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= rintfunc (pos[n]);
}

CRPBC::CRPBC(const DYNAMO::SimData* tmp):
  CBC(tmp, "RPBC", IC_purple)
{
  Sim = tmp;
}

void 
CRPBC::setPBC(Vector & pos) const
{ 
  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= Sim->aspectRatio[n] *
      rintfunc (pos[n]/Sim->aspectRatio[n]);    
}

void 
CRPBC::setPBC(Vector & pos, Vector&) const
{
  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= Sim->aspectRatio[n] *
      rintfunc (pos[n] / Sim->aspectRatio[n]);    
}

void 
CRPBC::setPBC(Vector  &pos, const Iflt&) const 
{
  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= Sim->aspectRatio[n] *
      rintfunc (pos[n] / Sim->aspectRatio[n]);    
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


CRNoXPBC::CRNoXPBC(const DYNAMO::SimData* tmp):
  CBC(tmp, "RNoXPBC",IC_purple)
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

void 
CRNoXPBC::setPBC(Vector & pos) const
{ 
  Iflt x = pos[0];

  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= Sim->aspectRatio[n] *
      rintfunc (pos[n] / Sim->aspectRatio[n]);    

  pos[0] = x;
}
  
void 
CRNoXPBC::setPBC(Vector & pos, Vector&) const
{ 
  Iflt x = pos[0];

  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= Sim->aspectRatio[n] *
      rintfunc (pos[n]/Sim->aspectRatio[n]);    

  pos[0] = x;
}

void 
CRNoXPBC::setPBC(Vector  &pos, const Iflt&) const 
{ 
  Iflt x = pos[0];

  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= Sim->aspectRatio[n] *
      rintfunc (pos[n] / Sim->aspectRatio[n]);    
  
  setPBC(pos); 

  pos[0] = x;
}
