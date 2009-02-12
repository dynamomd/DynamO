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

#include "gcellsShearing.hpp"
#include "globEvent.hpp"
#include "../NparticleEventData.hpp"
#include "../../extcode/xmlParser.h"
#include "../../extcode/xmlwriter.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../units/units.hpp"
#include "../ranges/1RAll.hpp"
#include "../ranges/1RNone.hpp"
#include "../ranges/2RSingle.hpp"
#include "../../schedulers/scheduler.hpp"
#include "../locals/local.hpp"
#include "../BC/LEBC.hpp"

CGCellsShearing::CGCellsShearing(DYNAMO::SimData* nSim, 
				 const std::string& name):
  CGCells2(nSim, "ShearingCells", NULL)
{
  globName = name;
  I_cout() << "Shearing Cells Loaded";
}

CGCellsShearing::CGCellsShearing(const XMLNode &XML, 
				 DYNAMO::SimData* ptrSim):
  CGCells2(ptrSim, "ShearingCells")
{
  operator<<(XML);

  I_cout() << "Cells in shearing Loaded";
}

void 
CGCellsShearing::initialise(size_t nID)
{
  ID=nID;
  
  if (dynamic_cast<const CLEBC *>(&(Sim->Dynamics.BCs())) == NULL)
    D_throw() << "You cannot use the shearing neighbour list"
	      << " in a system without Lees Edwards BC's";

  reinitialise(Sim->Dynamics.getLongestInteraction());
}

void 
CGCellsShearing::operator<<(const XMLNode& XML)
{
  try {
    if (XML.isAttributeSet("lambda"))
      lambda = boost::lexical_cast<Iflt>
	(XML.getAttribute("lambda"));
    
    globName = XML.getAttribute("Name");	
  }
  catch(...)
    {
      D_throw() << "Error loading CGCellsShearing";
    }
  
  if (lambda < 0.0 || lambda > 1.0)
    D_throw() << "Lambda out of bounds [0,1), lambda = " << lambda;
}

void
CGCellsShearing::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "ShearingCells"
      << xmlw::attr("Lambda") << lambda
      << xmlw::attr("Name") << globName;
}

void 
CGCellsShearing::runEvent(const CParticle& part) const
{

  size_t oldCell(partCellData[part.getID()].cell);

  //Determine the cell transition direction, its saved
  size_t cellDirection(Sim->Dynamics.Liouvillean().
		       getSquareCellCollision3
		       (part, cells[oldCell].origin, 
			cellDimension));

  if ((cellDirection != 1)
      || (cells[oldCell].coords[1] 
	  != (std::signbit(vel[1]) ? 0 : (cellCount[1] - 1))))
    {
      size_t cellpow(1);
      
      for (size_t iDim(0); iDim < cellDirection; ++iDim)
	cellpow *= cellCount[iDim];
      
      int velsign = 2 * (part.getVelocity()[cellDirection] > 0) - 1;
      int offset = (part.getVelocity()[cellDirection] > 0) * (cellCount[cellDirection] - 1);
      int endCell(oldCell + cellpow * velsign), inCell(oldCell + 2 * cellpow * velsign);
      
      {
	int tmpint = velsign * cellpow * cellCount[cellDirection];
	if (cells[oldCell].coords[cellDirection] == offset)
	  {
	    endCell -= tmpint;
	    inCell -= tmpint;
	  }
	else if (cells[oldCell].coords[cellDirection] == offset - velsign)
	  inCell -= tmpint;
      }

      removeFromCell(part.getID());
      addToCell(part.getID(), endCell);
    }
  else
    {
      Sim->Dynamics.Liouvillean().updateParticle(part);
      
      size_t oldCell(partCellData[part.getID()].cell);
      
      //Determine the cell transition direction, its saved
      size_t cellDirection(Sim->Dynamics.Liouvillean().
			   getSquareCellCollision3
			   (part, cells[oldCell].origin, 
			    cellDimension));
      
      /////////////////////////////EDIT MAGINOT LINE
      /*
      //long inPosition;
      size_t endCell;
      
      //This is required to get the correct sign on the velocity
      CVector<> rpos(part.getPosition() 
		 - cells[partCellData[part.getID()].cell].origin);

  CVector<> vel(part.getVelocity());

  Sim->Dynamics.BCs().setPBC(rpos, vel);
  
  //In this cell Event we must brute force periodic cell transistions
  if ((cellDirection == 1)
      && (cells[partCellData[part.getID()].cell].coords[1] 
	  == (std::signbit(vel[1]) ? 0 : (cellCount[1] - 1))))
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
      //Stream it to the boundary 
      //Recheck the dt      
      double dt = Sim->Dynamics.Liouvillean()
	.getSquareCellCollision2(part, 
				 cells[partCellData[part.getID()].cell].origin, 
				 cellDimension);

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
      endCell = getCellID(tmpPos);

      //Need to do a full update due to the channels of linked cells
      //at the boundaries being behind the particle when it wraps
      //around, just using cells in the direction of the motion misses
      //these

      removeFromCell(part.getID());
      addToCell(part.getID(), endCell);
      
      //Get rid of the virtual event that is next, update is delayed till
      //after all events are added
      Sim->ptrScheduler->popNextEvent();

      //Tell about the new locals
      BOOST_FOREACH(const size_t& lID, cells[endCell].locals)
	BOOST_FOREACH(const nbHoodSlot& nbs, sigNewLocalNotify)
	nbs.second(part, lID);

      //Tell about all particles in all linked cells, SLOW BUT SURE,
      BOOST_FOREACH(const int& nb, cells[endCell].neighbours)
	for (int next = cells[nb].list; next != -1; 
	     next = partCellData[next].next)
	  if (part.getID() != static_cast<size_t>(next))
	    BOOST_FOREACH(const nbHoodSlot& nbs,  sigNewNeighbourNotify)
	      nbs.second(part, next);
      
      //Push the next virtual event, this is the reason the scheduler
      //doesn't need a second callback
      Sim->ptrScheduler->pushEvent(part, getEvent(part));
      Sim->ptrScheduler->sort(part);
      
      BOOST_FOREACH(const nbHoodSlot& nbs, sigCellChangeNotify)
	nbs.second(part, oldCell);
  */
    }

  //Get rid of the virtual event that is next, update is delayed till
  //after all events are added
  Sim->ptrScheduler->popNextEvent();
  
  CVector<int> coords(cells[inCell].coords);
  
  //Particle has just arrived into a new cell warn the scheduler about
  //its new neighbours so it can add them to the heap
  //Holds the displacement in each dimension, the unit is cells!
  BOOST_STATIC_ASSERT(NDIM==3);
      
  //These are the two dimensions to walk in
  size_t dim1 = cellDirection + 1 - 3 * (cellDirection > 1),
    dim2 = cellDirection + 2 - 3 * (cellDirection > 0);
      
  size_t dim1pow(1), dim2pow(1);
      
  for (size_t iDim(0); iDim < dim1; ++iDim)
    dim1pow *= cellCount[iDim];
      
  for (size_t iDim(0); iDim < dim2; ++iDim)
    dim2pow *= cellCount[iDim];
      
  if (--coords[dim1] < 0) coords[dim1] = cellCount[dim1] - 1;
  if (--coords[dim2] < 0) coords[dim2] = cellCount[dim2] - 1;
      
  int nb(getCellIDprebounded(coords));
      
  //We now have the lowest cell coord, or corner of the cells to update
  for (int iDim(0); iDim < 3; ++iDim)
    {
      if (coords[dim2] + iDim == cellCount[dim2]) 
	nb -= dim2pow * cellCount[dim2];
	  
      for (int jDim(0); jDim < 3; ++jDim)
	{	  
	  if (coords[dim1] + jDim == cellCount[dim1]) 
	    nb -= dim1pow * cellCount[dim1];
	      
	  for (int next = cells[nb].list; next >= 0; 
	       next = partCellData[next].next)
	    BOOST_FOREACH(const nbHoodSlot& nbs, sigNewNeighbourNotify)
	      nbs.second(part, next);
	      
	  nb += dim1pow;
	}
	  
      if (coords[dim1] + 2 >= cellCount[dim1]) nb += dim1pow * cellCount[dim1];
	  
      nb += dim2pow - 3 * dim1pow;
    }
      
  //Tell about the new locals
  BOOST_FOREACH(const size_t& lID, cells[endCell].locals)
    BOOST_FOREACH(const nbHoodSlot& nbs, sigNewLocalNotify)
    nbs.second(part, lID);
      
  //Push the next virtual event, this is the reason the scheduler
  //doesn't need a second callback
  Sim->ptrScheduler->pushEvent(part, getEvent(part));
  Sim->ptrScheduler->sort(part);

  BOOST_FOREACH(const nbHoodSlot& nbs, sigCellChangeNotify)
    nbs.second(part, oldCell);
  
  //This doesn't stream the system as its a virtual event

  //Debug section
#ifdef DYNAMO_WallCollDebug
  {      
    CVector<int> tmp = cells[partCellData[part.getID()].cell].coords;
    CVector<int> tmp2 = cells[endCell].coords;
    
    std::cerr << "\nCGWall sysdt " 
	      << Sim->dSysTime / Sim->Dynamics.units().unitTime()
	      << "  WALL ID "
	      << part.getID()
	      << "  from <" 
	      << tmp[0] << "," << tmp[1] << "," << tmp[2]
	      << "> to <" 
	      << tmp2[0] << "," << tmp2[1] << "," << tmp2[2] << ">";
  }
#endif

}

void 
CGCellsShearing::getParticleNeighbourhood(const CParticle& part,
				   const nbHoodFunc& func) const
{
  CGCells2::getParticleNeighbourhood(part, func);
  
  CVector<int> coords(cells[partCellData[part.getID()].cell].coords);
  if ((coords[1] == 0) || (coords[1] == cellCount[1] - 1))
    getExtraLEParticleNeighbourhood(part, func);

}

void 
CGCellsShearing::getExtraLEParticleNeighbourhood(const CParticle& part,
					  const nbHoodFunc& func) const
{
  size_t cellID = partCellData[part.getID()].cell;
  CVector<int> coords(cells[cellID].coords);
  
#ifdef DYNAMO_DEBUG
  if ((coords[1] != 0) && (coords[1] != cellCount[1] - 1))
    D_throw() << "Shouldn't call this function unless the particle is at a border in the y dimension";
#endif 

  //Move to the bottom of x
  coords[0] = 0;
  cellID -= coords[0];
  
  //Get the y plane to test in
  if (coords[1])
    {
      //At the top, assuming the above debug statement won't trigger
      coords[1] = 0;
      cellID -= cellCount[0] * (cellCount[1] - 1);
    }
  else
    {
      //At the bottom
      coords[1] = cellCount[1] - 1;
      cellID += cellCount[0] * (cellCount[1] - 1);
    }
  
  ////Move a single cell down in Z
  if (coords[2])
    {
      //Got room to move down
      --coords[2];
      cellID -= cellCount[0] * cellCount[1];
    }
  else
    {
      //Already at the bottom of the Z component, need to wrap
      coords[2] = cellCount[2] - 1;
      cellID += cellCount[0] * cellCount[1] * (cellCount[2] - 1);
    }
  
  
  for (int i(0); i < 3; ++i)
    {
      if (coords[2] + i == cellCount[2]) cellID -= cellCount[0] * cellCount[1] * cellCount[2];

      for (int j(0); j < cellCount[0]; ++j)
	{
	  for (int next(cells[cellID].list);
	       next >= 0; next = partCellData[next].next)
	    if (next != int(part.getID()))
	      func(part, next);

	  ++cellID;
	}
    
      cellID += cellCount[0] * (cellCount[1] - 1);
    }
}

