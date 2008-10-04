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


CRLEBC::CRLEBC(DYNAMO::SimData* tmp):
  Base_Class("LEBC",IC_purple),
  dxd(0.0) 
{
  Sim = tmp;
  I_cout() << "Rectangular Lee's Edwards BC loaded"; 
}

CRLEBC::CRLEBC(const XMLNode& XML, DYNAMO::SimData* tmp):
  Base_Class("LEBC",IC_purple),
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
CRLEBC::setPBC(CVector<> &pos) const 
{
  //Shift the x distance due to the Lee's Edwards conditions
  pos.data[0] -= rint(pos.data[1]/Sim->aspectRatio[1])*dxd;
  
  rounding(pos);
}

void 
CRLEBC::setPBC(CVector<> &pos, CVector<> &vel) const 
{
  //Shift the x distance due to the Lee's Edwards conditions
  pos.data[0] -= rint(pos.data[1]/Sim->aspectRatio[1]) * dxd;
  
  //Adjust the velocity due to the box shift
  vel.data[0] -= rint(pos.data[1]/Sim->aspectRatio[1]) * ShearRate 
    * Sim->aspectRatio[1];
  
  rounding(pos);
}

void 
CRLEBC::setPBC(CVector<> &posVec, Iflt dt) const 
{ 
  Iflt localdxd = dxd + dt * ShearRate * Sim->aspectRatio[1];
  
  //Shift the x distance due to the Lee's Edwards conditions
  posVec.data[0] -= rint(posVec.data[1]/Sim->aspectRatio[1]) * localdxd;
  
  rounding(posVec);
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

CSLEBC::CSLEBC(DYNAMO::SimData* tmp):
  Base_Class("LEBC",IC_purple),
  dxd(0.0) 
{
  Sim = tmp;
  I_cout() << "Square Lee's Edwards BC loaded"; 
}

CSLEBC::CSLEBC(const XMLNode& XML, DYNAMO::SimData* tmp):
  Base_Class("LEBC",IC_purple),
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
CSLEBC::setPBC(CVector<> &pos) const 
{
  //Shift the x distance due to the Lee's Edwards conditions
  pos.data[0] -= rint(pos.data[1])*dxd;
  
  rounding(pos);
}

void 
CSLEBC::setPBC(CVector<> &pos, CVector<> &vel) const 
{
  //Shift the x distance due to the Lee's Edwards conditions
  pos.data[0] -= rint(pos.data[1]) * dxd;
  
  //Adjust the velocity due to the box shift
  vel.data[0] -= rint(pos.data[1]) * ShearRate;
  
  rounding(pos);
}

void 
CSLEBC::setPBC(CVector<> &posVec, Iflt dt) const 
{
  Iflt localdxd = dxd + dt * ShearRate;
  
  //Shift the x distance due to the Lee's Edwards conditions
  posVec.data[0] -= rint(posVec.data[1]) * localdxd;
  
  rounding(posVec);
}

void 
CSLEBC::update(const Iflt& dt) 
{
  //Shift the boundary of the system v_box = \gamma*L
  dxd += dt * ShearRate;
  
  //PBC for the shift to keep accuracy?
  dxd -= floor(dxd);
}

