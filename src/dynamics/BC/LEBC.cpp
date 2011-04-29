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
#include "../interactions/intEvent.hpp"
#include "../../base/is_simdata.hpp"
#include "../../extcode/mathtemplates.hpp"
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <cmath>

BCLeesEdwards::BCLeesEdwards(const DYNAMO::SimData* tmp):
  BoundaryCondition(tmp, "LEBC",IC_purple),
  dxd(0.0) 
{
  Sim = tmp;
  I_cout() << " Lee's Edwards BC loaded"; 
}

BCLeesEdwards::BCLeesEdwards(const magnet::xml::Node& XML, 
						   const DYNAMO::SimData* tmp):
  BoundaryCondition(tmp, "LEBC",IC_purple),
  dxd(0.0) 
{
  Sim = tmp;
  operator<<(XML);
  I_cout() << " Lee's Edwards BC loaded"; 
  I_cout() << "DXD = " << dxd;
}

void 
BCLeesEdwards::outputXML(xml::XmlStream &XML) const
{
  XML << xml::attr("Type") << "LE"
      << xml::attr("DXD") << dxd;
}

void 
BCLeesEdwards::operator<<(const magnet::xml::Node& XML)
{ 
  try 
    {
      if (XML.getAttribute("DXD").valid())
	dxd = XML.getAttribute("DXD").as<double>();
    }
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in LEBC";
    }
}

BoundaryCondition* 
BCLeesEdwards::Clone () const 
{ return new BCLeesEdwards(*this); }

void 
BCLeesEdwards::applyBC(Vector  &pos) const 
{
  //Shift the x distance due to the Lee's Edwards conditions
  pos[0] -= rint(pos[1] / Sim->primaryCellSize[1])*dxd;
  
  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= Sim->primaryCellSize[n] *
      rintfunc (pos[n] / Sim->primaryCellSize[n]);    
}

void 
BCLeesEdwards::applyBC(Vector  &pos, Vector &vel) const 
{
  //Shift the x distance due to the Lee's Edwards conditions
  pos[0] -= rint(pos[1] / Sim->primaryCellSize[1]) * dxd;
  
  //Adjust the velocity due to the box shift
  vel[0] -= rint(pos[1] / Sim->primaryCellSize[1]) 
    * shearRate() * Sim->primaryCellSize[1];
  
  for (size_t n = 0; n < NDIM; ++n)
    pos[n] -= Sim->primaryCellSize[n] *
      rintfunc (pos[n] / Sim->primaryCellSize[n]);    
}

void 
BCLeesEdwards::applyBC(Vector& posVec, const double& dt) const 
{ 
  double localdxd = dxd + dt * shearRate() * Sim->primaryCellSize[1];
  
  //Shift the x distance due to the Lee's Edwards conditions
  posVec[0] -= rint(posVec[1] / Sim->primaryCellSize[1]) * localdxd;
  
  for (size_t n = 0; n < NDIM; ++n)
    posVec[n] -= Sim->primaryCellSize[n] *
      rintfunc (posVec[n] / Sim->primaryCellSize[n]);    
}

void 
BCLeesEdwards::update(const double& dt) 
{
  //Shift the boundary of the system v_box = \gamma*L
  dxd += dt * shearRate() * Sim->primaryCellSize[1];
  
  //PBC for the shift to keep accuracy?
  dxd -= floor(dxd/Sim->primaryCellSize[0])*Sim->primaryCellSize[0];
}
