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

#include "socells.hpp"
#include "globEvent.hpp"
#include "../NparticleEventData.hpp"
#include "../../extcode/xmlParser.h"
#include "../../extcode/xmlwriter.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../units/units.hpp"
#include "../ranges/1RAll.hpp"
#include "../../schedulers/scheduler.hpp"
#include "../locals/local.hpp"
#include "../BC/LEBC.hpp"
#include <boost/static_assert.hpp>

CGSOCells::CGSOCells(DYNAMO::SimData* nSim, const std::string& name):
  CGlobal(nSim, "SingleOccupancyCells"),
  cellCount(0),
  cellDimension(1,1,1),
  MaxIntDist(0.0)
{
  globName = name;
  I_cout() << "Single occupancy cells loaded";
}

CGSOCells::CGSOCells(const XMLNode &XML, DYNAMO::SimData* ptrSim):
  CGlobal(ptrSim, "SingleOccupancyCells"),
  cellCount(0),
  cellDimension(1,1,1),
  MaxIntDist(0.0)
{
  operator<<(XML);

  I_cout() << "Single occupancy cells loaded";
}

void 
CGSOCells::operator<<(const XMLNode& XML)
{
  try {    
    globName = XML.getAttribute("Name");	
  }
  catch(...)
    {
      D_throw() << "Error loading CGSOCells";
    }
}

CGlobEvent 
CGSOCells::getEvent(const CParticle& part) const
{
#ifdef ISSS_DEBUG
  if (!Sim->Dynamics.Liouvillean().isUpToDate(part))
    D_throw() << "Particle is not up to date";
#endif

  //This 
  //Sim->Dynamics.Liouvillean().updateParticle(part);
  //is not required as we compensate for the delay using 
  //Sim->Dynamics.Liouvillean().getParticleDelay(part)

  D_throw() << "Not implemented";
  
  /*
    return CGlobEvent(part,
		    Sim->Dynamics.Liouvillean().
		    getSquareCellCollision2
		    (part, cells[partCellData[part.getID()].cell].origin, 
		     cellDimension)
		    -Sim->Dynamics.Liouvillean().getParticleDelay(part)
		    ,
		    CELL, *this);
  */
}

void
CGSOCells::runEvent(const CParticle& part) const
{
}

void 
CGSOCells::initialise(size_t nID)
{
  ID=nID;

}

void 
CGSOCells::outputXML(xmlw::XmlStream& XML) const
{
  outputXML(XML, "Cells");
}

void
CGSOCells::outputXML(xmlw::XmlStream& XML, const std::string& name) const
{
  XML << xmlw::attr("Type") << name
      << xmlw::attr("Name") << globName
      << xmlw::attr("CellWidth") 
      << MaxIntDist / Sim->Dynamics.units().unitLength();
}
