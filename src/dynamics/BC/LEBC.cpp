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

#include "LEBC.hpp"
#include "../../extcode/xmlParser.h"
#include "../interactions/intEvent.hpp"
#include "../../base/is_simdata.hpp"
#include "../../extcode/mathtemplates.hpp"
#include "../units/shear.hpp"
#include <boost/lexical_cast.hpp>
#include <magnet/xmlwriter.hpp>
#include <cmath>

BCRectangularLeesEdwards::BCRectangularLeesEdwards(const DYNAMO::SimData* tmp):
  BoundaryCondition(tmp, "LEBC",IC_purple),
  dxd(0.0) 
{
  Sim = tmp;
  I_cout() << "Rectangular Lee's Edwards BC loaded"; 
}

BCRectangularLeesEdwards::BCRectangularLeesEdwards(const XMLNode& XML, const DYNAMO::SimData* tmp):
  BoundaryCondition(tmp, "LEBC",IC_purple),
  dxd(0.0) 
{
  Sim = tmp;
  operator<<(XML);
  I_cout() << "Rectangular Lee's Edwards BC loaded"; 
  I_cout() << "DXD = " << dxd;
}

void 
BCRectangularLeesEdwards::outputXML(xml::XmlStream &XML) const
{
  XML << xml::attr("Shape") << "Rectangular"
      << xml::attr("Boundary") << "LE"
      << xml::attr("DXD") << dxd;
}

void 
BCRectangularLeesEdwards::operator<<(const XMLNode& XML)
{ 
  try 
    {
      if (XML.isAttributeSet("DXD"))
	dxd = boost::lexical_cast<double>(XML.getAttribute("DXD"));
    }
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in LEBC";
    }
}

BoundaryCondition* 
BCRectangularLeesEdwards::Clone () const 
{ return new BCRectangularLeesEdwards(*this); }

void 
BCRectangularLeesEdwards::applyBC(Vector  &pos) const 
{
  //Shift the x distance due to the Lee's Edwards conditions
  pos[0] -= rint(pos[1] / Sim->aspectRatio[1])*dxd;
  
  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= Sim->aspectRatio[n] *
      rintfunc (pos[n] / Sim->aspectRatio[n]);    
}

void 
BCRectangularLeesEdwards::applyBC(Vector  &pos, Vector &vel) const 
{
  //Shift the x distance due to the Lee's Edwards conditions
  pos[0] -= rint(pos[1] / Sim->aspectRatio[1]) * dxd;
  
  //Adjust the velocity due to the box shift
  vel[0] -= rint(pos[1] / Sim->aspectRatio[1]) 
    * UShear::ShearRate() * Sim->aspectRatio[1];
  
  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= Sim->aspectRatio[n] *
      rintfunc (pos[n] / Sim->aspectRatio[n]);    
}

void 
BCRectangularLeesEdwards::applyBC(Vector  &posVec, const double& dt) const 
{ 
  double localdxd = dxd + dt * UShear::ShearRate() * Sim->aspectRatio[1];
  
  //Shift the x distance due to the Lee's Edwards conditions
  posVec[0] -= rint(posVec[1] / Sim->aspectRatio[1]) * localdxd;
  
  for (size_t n = 0; n < NDIM; ++n)
    posVec[n] -= Sim->aspectRatio[n] *
      rintfunc (posVec[n] / Sim->aspectRatio[n]);    
}

void 
BCRectangularLeesEdwards::update(const double& dt) 
{
  //Shift the boundary of the system v_box = \gamma*L
  dxd += dt * UShear::ShearRate() * Sim->aspectRatio[1];
  
  //PBC for the shift to keep accuracy?
  dxd -= floor(dxd/Sim->aspectRatio[0])*Sim->aspectRatio[0];
}


/////////////////////////////Rectangular/////////////////////////////////////

BCSquareLeesEdwards::BCSquareLeesEdwards(const DYNAMO::SimData* tmp):
  BoundaryCondition(tmp, "LEBC",IC_purple),
  dxd(0.0) 
{
  Sim = tmp;
  I_cout() << "Square Lee's Edwards BC loaded"; 
}

BCSquareLeesEdwards::BCSquareLeesEdwards(const XMLNode& XML, const DYNAMO::SimData* tmp):
  BoundaryCondition(tmp, "LEBC",IC_purple),
  dxd(0.0) 
{
  Sim = tmp;
  operator<<(XML);
  I_cout() << "Square Lee's Edwards BC loaded"; 
  I_cout() << "DXD = " << dxd;
}

void 
BCSquareLeesEdwards::outputXML(xml::XmlStream &XML) const
{
  XML << xml::attr("Shape") << "Square"
      << xml::attr("Boundary") << "LE"
      << xml::attr("DXD") << dxd;
}

void 
BCSquareLeesEdwards::operator<<(const XMLNode& XML)
{ 
  try 
    {
      if (XML.isAttributeSet("DXD"))
	dxd = boost::lexical_cast<double>(XML.getAttribute("DXD"));
    }
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in LEBC";
    }
}

BoundaryCondition* 
BCSquareLeesEdwards::Clone () const 
{ return new BCSquareLeesEdwards(*this); }

void 
BCSquareLeesEdwards::applyBC(Vector  &pos) const 
{
  //Shift the x distance due to the Lee's Edwards conditions
  pos[0] -= rint(pos[1]) * dxd;
  
  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= rintfunc (pos[n]);
}

void 
BCSquareLeesEdwards::applyBC(Vector  &pos, Vector &vel) const 
{
  //Shift the x distance due to the Lee's Edwards conditions
  pos[0] -= rint(pos[1]) * dxd;
  
  //Adjust the velocity due to the box shift
  vel[0] -= rint(pos[1]) * UShear::ShearRate();
  
  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= rintfunc (pos[n]);
}

void 
BCSquareLeesEdwards::applyBC(Vector  &posVec, const double& dt) const 
{
  double localdxd = dxd + dt * UShear::ShearRate();
  
  //Shift the x distance due to the Lee's Edwards conditions
  posVec[0] -= rint(posVec[1]) * localdxd;
  
  for (size_t n = 0; n < NDIM; ++n)
    posVec[n] -= rintfunc (posVec[n]);
}

void 
BCSquareLeesEdwards::update(const double& dt) 
{
  //Shift the boundary of the system v_box = \gamma*L
  dxd += dt * UShear::ShearRate();
  
  //PBC for the shift to keep accuracy?
  dxd -= floor(dxd);
}

