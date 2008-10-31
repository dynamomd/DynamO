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

#include "multlist_shear.hpp"
#include "../dynamics/interactions/intEvent.hpp"
#include "../extcode/threadpool.hpp"
#include "../simulation/particle.hpp"
#include "../dynamics/dynamics.hpp"
#include "../dynamics/liouvillean/liouvillean.hpp"
#include "../dynamics/BC/BC.hpp"
#include "../base/is_simdata.hpp"
#include "../base/is_base.hpp"
#include "../dynamics/globals/globEvent.hpp"
#include "../dynamics/systems/system.hpp"
#include <cmath> //for huge val
#include <boost/lexical_cast.hpp>
#include "../extcode/xmlParser.h"
#include "../dynamics/globals/global.hpp"
#include "../dynamics/globals/globEvent.hpp"
#include "../dynamics/ranges/1range.hpp"
#include "../dynamics/ranges/2RSingle.hpp"
#include "../dynamics/BC/LEBC.hpp"

CSMultListShear::CSMultListShear(const XMLNode& XML, const DYNAMO::SimData* Sim):
  CSMultList(Sim, "MultListCellular")
{
  I_cout() << "Multi List Cellular Algorithm with shearing";
  operator<<(XML);
}

CSMultListShear::CSMultListShear(const DYNAMO::SimData* Sim):
  CSMultList(Sim, "MultListCellular")
{ I_cout() << "Multi List Cellular Algorithmn with shearing"; }

void
CSMultListShear::initialise()
{
  if (!Sim->Dynamics.BCTypeTest<CRLEBC>() 
      && !Sim->Dynamics.BCTypeTest<CSLEBC>())
    D_throw() << "This scheduler isn't suitable for non-sheared systems\n" 
      "it would work but you may as well use the normal multlist scheduler";

  reinitialise(Sim->Dynamics.getLongestInteraction()); 
}

void
CSMultListShear::reinitialise(Iflt maxdiam)
{
  I_cout() << "Reinitialising on collision " << Sim->lNColl;
  //Absolutely cannot have overlapping cells  
  lambda = 0.0;
  eventHeap.clear();
  eventHeap.resize(Sim->lN);
  eventCount.clear();
  eventCount.resize(Sim->lN, 0);
  
  //Create the cells
  addCells(maxdiam, false);

  //Link up the LE cells
  link_LE_cells();
  
  //Now initialise the interactions
  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    addNewEvents_init(part);
  
  eventHeap.init();
  
#ifndef CBT
  I_cout() << "BPQ: Number of lists " << eventHeap.NLists();
  I_cout() << "BPQ: Scale Factor " << eventHeap.scaleFactor();
#endif
}

void 
CSMultListShear::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "MultListShear";
  CSCells::outputXML(XML);
}

void 
CSMultListShear::operator<<(const XMLNode& XML)
{
  CSCells::operator<<(XML);
}

void 
CSMultListShear::cellEvent(const CParticle& part) const 
{  
  //Determine the cell transition direction, its saved
  int cellDirection(eventHeap[part.getID()].top().collCounter2);
  //long inPosition;
  int endCell;
  
  Sim->Dynamics.Liouvillean().updateParticle(part);

  //This is required to get the correct sign on the velocity
  CVector<> rpos(part.getPosition() - cells[partCellData[part.getID()].cell].origin);
  CVector<> vel(part.getVelocity());
  Sim->Dynamics.BCs().setPBC(rpos, vel);

  //In this cell Event we must brute force periodic cell transistions
  if (((cells[partCellData[part.getID()].cell].coords[1] == 0)
      && (cellDirection == 1) && std::signbit(vel[1]))
      ||
      ((cells[partCellData[part.getID()].cell].coords[1] == cellCount[1] - 1)
       && (cellDirection == 1) && !std::signbit(vel[1])))
    {
      //Debug section
#ifdef DYNAMO_WallCollDebug
      std::cerr << "\nBoundary transition ";

      if (std::signbit(vel[1])) 
	std::cerr << "Down";
      else
	std::cerr << "Up";
#endif      
      //Bottom heading down
      //Stream it to the boundary //Recheck the dt
      
      double dt = Sim->Dynamics.Liouvillean()
	.getSquareCellCollision(part, cells[partCellData[part.getID()].cell].origin, 
				cellDimension).dt;

      Sim->Dynamics.Liouvillean().advanceUpdateParticle(part, dt);
      
      CVector<> tmpPos = part.getPosition();

      //Add enough of a step to move it into the other cell
      if (std::signbit(vel[1]))
	tmpPos[1] -= 0.5 * cellDimension[1];
      else
	tmpPos[1] += 0.5 * cellDimension[1];
      
      //Now a special predictive setpbc must be used
      Sim->Dynamics.BCs().setPBC(tmpPos, dt);
      //Now use the end coordinates to give you the final cell
      endCell = getID(tmpPos);
    }
  else if (std::signbit(vel[cellDirection])) 
      endCell = cells[partCellData[part.getID()].cell].negCells[cellDirection];
  else
      endCell = cells[partCellData[part.getID()].cell].posCells[cellDirection];

  //Debug section
#ifdef DYNAMO_WallCollDebug
  {      
    CVector<long> tmp = cells[partCellData[part.getID()].cell].coords;
    CVector<long> tmp2 = cells[endCell].coords;
    
    std::cerr << "\nsysdt " 
	      << (eventHeap.next_dt() + Sim->dSysTime)
	      << "  WALL ID "
	      << part.getID()
	      << "  dt " << (eventHeap.next_dt())
	      << "  from <" 
	      << tmp[0] << "," << tmp[1] << "," << tmp[2]
	      << "> to <" 
	      << tmp2[0] << "," << tmp2[1] << "," << tmp2[2] << ">";
  }
#endif  

  //Need to do a full update due to the channels of linked cells
  //at the boundaries being behind the particle when it wraps
  //around, just using cells in the direction of the motion misses
  //these

  removeFromCell(part.getID());
  addToCell(part.getID(), endCell);

  eventHeap[part.getID()].clear();
  addNewEvents(part);
  eventHeap.update(part.getID());
  
  //This is slower! as I haven't used the normal updates for other
  //times, however, I just want this to work!
}
