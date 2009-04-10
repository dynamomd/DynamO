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

#include "LEBC.hpp"
#include <boost/lexical_cast.hpp>
#include <cmath>
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../interactions/intEvent.hpp"
#include "../../base/is_exception.hpp"
#include "../../base/is_simdata.hpp"
#include "rintfunc.hpp"

CRLEBC::CRLEBC(const DYNAMO::SimData* tmp):
  CBC(tmp, "LEBC",IC_purple),
  dxd(0.0) 
{
  Sim = tmp;
  I_cout() << "Rectangular Lee's Edwards BC loaded"; 
}

CRLEBC::CRLEBC(const XMLNode& XML, const DYNAMO::SimData* tmp):
  CBC(tmp, "LEBC",IC_purple),
  dxd(0.0) 
{
  Sim = tmp;
  operator<<(XML);
  I_cout() << "Rectangular Lee's Edwards BC loaded"; 
  I_cout() << "DXD = " << dxd;
}

void 
CRLEBC::outputXML(xmlw::XmlStream &XML) const
{
  XML << xmlw::attr("Shape") << "Rectangular"
      << xmlw::attr("Boundary") << "LE"
      << xmlw::attr("DXD") << dxd;
}

void 
CRLEBC::operator<<(const XMLNode& XML)
{ 
  try 
    {
      dxd = boost::lexical_cast<Iflt>(XML.getAttribute("DXD"));
    }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in LEBC";
    }
}

CBC* 
CRLEBC::Clone () const 
{ return new CRLEBC(*this); }

void 
CRLEBC::setPBC(Vector  &pos) const 
{
  //Shift the x distance due to the Lee's Edwards conditions
  pos[0] -= rint(pos[1] / Sim->aspectRatio[1])*dxd;
  
  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= Sim->aspectRatio[n] *
      rintfunc (pos[n] / Sim->aspectRatio[n]);    
}

void 
CRLEBC::setPBC(Vector  &pos, Vector &vel) const 
{
  //Shift the x distance due to the Lee's Edwards conditions
  pos[0] -= rint(pos[1] / Sim->aspectRatio[1]) * dxd;
  
  //Adjust the velocity due to the box shift
  vel[0] -= rint(pos[1] / Sim->aspectRatio[1]) 
    * ShearRate * Sim->aspectRatio[1];
  
  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= Sim->aspectRatio[n] *
      rintfunc (pos[n] / Sim->aspectRatio[n]);    
}

void 
CRLEBC::setPBC(Vector  &posVec, const Iflt& dt) const 
{ 
  Iflt localdxd = dxd + dt * ShearRate * Sim->aspectRatio[1];
  
  //Shift the x distance due to the Lee's Edwards conditions
  posVec[0] -= rint(posVec[1] / Sim->aspectRatio[1]) * localdxd;
  
  for (size_t n = 0; n < NDIM; ++n)
    posVec[n] -= Sim->aspectRatio[n] *
      rintfunc (posVec[n] / Sim->aspectRatio[n]);    
}

void 
CRLEBC::update(const Iflt& dt) 
{
  //Shift the boundary of the system v_box = \gamma*L
  dxd += dt * ShearRate * Sim->aspectRatio[1];
  
  //PBC for the shift to keep accuracy?
  dxd -= floor(dxd/Sim->aspectRatio[0])*Sim->aspectRatio[0];
}


/////////////////////////////Rectangular/////////////////////////////////////

CSLEBC::CSLEBC(const DYNAMO::SimData* tmp):
  CBC(tmp, "LEBC",IC_purple),
  dxd(0.0) 
{
  Sim = tmp;
  I_cout() << "Square Lee's Edwards BC loaded"; 
}

CSLEBC::CSLEBC(const XMLNode& XML, const DYNAMO::SimData* tmp):
  CBC(tmp, "LEBC",IC_purple),
  dxd(0.0) 
{
  Sim = tmp;
  operator<<(XML);
  I_cout() << "Square Lee's Edwards BC loaded"; 
  I_cout() << "DXD = " << dxd;
}

void 
CSLEBC::outputXML(xmlw::XmlStream &XML) const
{
  XML << xmlw::attr("Shape") << "Square"
      << xmlw::attr("Boundary") << "LE"
      << xmlw::attr("DXD") << dxd;
}

void 
CSLEBC::operator<<(const XMLNode& XML)
{ 
  try 
    {
      dxd = boost::lexical_cast<Iflt>(XML.getAttribute("DXD"));
    }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in LEBC";
    }
}

CBC* 
CSLEBC::Clone () const 
{ return new CSLEBC(*this); }

void 
CSLEBC::setPBC(Vector  &pos) const 
{
  //Shift the x distance due to the Lee's Edwards conditions
  pos[0] -= rint(pos[1]) * dxd;
  
  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= rintfunc (pos[n]);
}

void 
CSLEBC::setPBC(Vector  &pos, Vector &vel) const 
{
  //Shift the x distance due to the Lee's Edwards conditions
  pos[0] -= rint(pos[1]) * dxd;
  
  //Adjust the velocity due to the box shift
  vel[0] -= rint(pos[1]) * ShearRate;
  
  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= rintfunc (pos[n]);
}

void 
CSLEBC::setPBC(Vector  &posVec, const Iflt& dt) const 
{
  Iflt localdxd = dxd + dt * ShearRate;
  
  //Shift the x distance due to the Lee's Edwards conditions
  posVec[0] -= rint(posVec[1]) * localdxd;
  
  for (size_t n = 0; n < NDIM; ++n)
    posVec[n] -= rintfunc (posVec[n]);
}

void 
CSLEBC::update(const Iflt& dt) 
{
  //Shift the boundary of the system v_box = \gamma*L
  dxd += dt * ShearRate;
  
  //PBC for the shift to keep accuracy?
  dxd -= floor(dxd);
}

