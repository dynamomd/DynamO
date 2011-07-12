/*  dynamo:- Event driven molecular dynamics simulator 
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

#include "gcellsShearing.hpp"
#include "globEvent.hpp"
#include "../NparticleEventData.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../units/units.hpp"
#include "../ranges/1RAll.hpp"
#include "../ranges/1RNone.hpp"
#include "../ranges/2RSingle.hpp"
#include "../../schedulers/scheduler.hpp"
#include "../locals/local.hpp"
#include "../BC/LEBC.hpp"
#include "../liouvillean/NewtonianGravityL.hpp"
#include <magnet/xmlwriter.hpp>

CGCellsShearing::CGCellsShearing(dynamo::SimData* nSim, 
				 const std::string& name):
  CGCellsMorton(nSim, "ShearingCells", NULL)
{
  globName = name;
  dout << "Shearing Cells Loaded" << std::endl;
}

CGCellsShearing::CGCellsShearing(const magnet::xml::Node& XML, 
				 dynamo::SimData* ptrSim):
  CGCellsMorton(ptrSim, "ShearingCells")
{
  operator<<(XML);

  dout << "Cells in shearing Loaded" << std::endl;
}

void 
CGCellsShearing::initialise(size_t nID)
{
  ID=nID;
   
  if (!(Sim->dynamics.BCTypeTest<BCLeesEdwards>()))
    M_throw() << "You cannot use the shearing neighbour list in a system without Lees Edwards BC's";
  
  if (overlink != 1) M_throw() << "Cannot shear with overlinking yet";
  
  reinitialise(Sim->dynamics.getLongestInteraction());
}

void
CGCellsShearing::outputXML(magnet::xml::XmlStream& XML) const
{
  CGCellsMorton::outputXML(XML, std::string("ShearingCells"));
}

GlobalEvent 
CGCellsShearing::getEvent(const Particle& part) const
{
#ifdef ISSS_DEBUG
  if (!Sim->dynamics.getLiouvillean().isUpToDate(part))
    M_throw() << "Particle is not up to date";
#endif

  //We do not inherit GCells get Event as the calcPosition thing done
  //for infinite systems is breaking it for shearing for some reason.
  return GlobalEvent(part,
		     Sim->dynamics.getLiouvillean().
		     getSquareCellCollision2
		     (part, 
		      calcPosition(partCellData[part.getID()].cell, part),
		      cellDimension)
		     -Sim->dynamics.getLiouvillean().getParticleDelay(part)
		     ,
		     CELL, *this);
}

void 
CGCellsShearing::runEvent(const Particle& part, const double) const
{
  Sim->dynamics.getLiouvillean().updateParticle(part);

  size_t oldCell(partCellData[part.getID()].cell);
   oldCellCoords;
  Vector oldCellOrigin = calcPosition(oldCell, part);

  //Determine the cell transition direction, its saved
  int cellDirectionInt(Sim->dynamics.getLiouvillean().
		       getSquareCellCollision3
		       (part, oldCellOrigin,
			cellDimension));
  

  Vector pos(part.getPosition() - oldCellOrigin), vel(part.getVelocity());
  Sim->dynamics.BCs().applyBC(pos, vel);

  size_t cellDirection = abs(cellDirectionInt) - 1;

  int endCell(-1); //The ID of the cell the particle enters

  if ((cellDirection == 1) &&
      (cells[oldCell].coords[1] == ((cellDirectionInt < 0) ? 0 : (cellCount[1] - 1))))
    {
      //We're wrapping in the y direction, we have to compute
      //which cell its entering
      
      //Calculate the final x value
      //Time till transition, assumes the particle is up to date
      double dt = Sim->dynamics.getLiouvillean().getSquareCellCollision2
	(part, oldCellOrigin, 
	 cellDimension);
     
      //Remove the old x contribution
      endCell = oldCell - cells[oldCell].coords[0];
      
      //Update the y dimension
      if (cellDirectionInt < 0)
	endCell += cellCount[0] * (cellCount[1]-1);
      else
	endCell -= cellCount[0] * (cellCount[1]-1);

      //Predict the position of the particle in the x dimension
      Sim->dynamics.getLiouvillean().advanceUpdateParticle(part, dt);
      Vector  tmpPos = part.getPosition();

      //This just ensures we wrap the image
      if (cellDirectionInt < 0)
	tmpPos[1] -= 0.5 * cellDimension[1];
      else
	tmpPos[1] += 0.5 * cellDimension[1];

      //This rewinds the particle again
      Sim->dynamics.getLiouvillean().updateParticle(part);

      //Determine the x position (in cell coords) of the particle and
      //add it to the endCellID
      Sim->dynamics.BCs().applyBC(tmpPos, dt);
      endCell += int((tmpPos[0] + 0.5 * Sim->primaryCellSize[0]) 
		     / cellLatticeWidth[0]);

      removeFromCell(part.getID());
      addToCell(part.getID(), endCell);
      
      //Get rid of the virtual event that is next, update is delayed till
      //after all events are added
      Sim->ptrScheduler->popNextEvent();

      //Check the entire neighbourhood, could check just the new
      //neighbours and the extra LE neighbourhood strip but its a lot
      //of code
      if (isUsedInScheduler)
	getParticleNeighbourhood(part, magnet::function::MakeDelegate(Sim->ptrScheduler, 
								  &CScheduler::addInteractionEvent));
      
      BOOST_FOREACH(const nbHoodSlot& nbs, sigNewNeighbourNotify)
	getParticleNeighbourhood(part, nbs.second);
      

    }
  else if ((cellDirection == 1) && 
	   (cells[oldCell].coords[1] == ((cellDirectionInt < 0) ? 1 : (cellCount[1] - 2))))
    {
      //We're entering the boundary of the y direction
      //Calculate the end cell, no boundary wrap check required
      endCell = oldCell + cellCount[0] * ((cellDirectionInt < 0) ? -1 : 1);
      
      removeFromCell(part.getID());
      addToCell(part.getID(), endCell);
      
      //Get rid of the virtual event that is next, update is delayed till
      //after all events are added
      Sim->ptrScheduler->popNextEvent();
      
      //Check the extra LE neighbourhood strip
      if (isUsedInScheduler)
	getExtraLEParticleNeighbourhood(part, magnet::function::MakeDelegate(Sim->ptrScheduler,
									 &CScheduler::addInteractionEvent));
      
      BOOST_FOREACH(const nbHoodSlot& nbs, sigNewNeighbourNotify)
	getExtraLEParticleNeighbourhood(part, nbs.second);
    }
  else
    {
      //Here we follow the same procedure (except one more if statement) as the original cell list for new neighbours
      size_t cellpow(1);
      
      for (size_t iDim(0); iDim < cellDirection; ++iDim)
	cellpow *= cellCount[iDim];
      
      int velsign = 2 * (cellDirectionInt > 0) - 1;
      int offset = (cellDirectionInt > 0) * (cellCount[cellDirection] - 1);
      
      endCell = oldCell + cellpow * velsign;
      int inCell = oldCell + 2 * cellpow * velsign;
      
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

      //Get rid of the virtual event that is next, update is delayed till
      //after all events are added
      Sim->ptrScheduler->popNextEvent();
      
      if ((cellDirection == 2) 
	  && ((cells[oldCell].coords[1] == 0) 
	      || (cells[oldCell].coords[1] == cellCount[1] - 1)))
	//We're at the boundary moving in the z direction, we must
	//add the new LE strips as neighbours	
	//We just check the entire Extra LE neighbourhood
	{
	  if (isUsedInScheduler)
	    getExtraLEParticleNeighbourhood(part, magnet::function::MakeDelegate(Sim->ptrScheduler,
									     &CScheduler::addInteractionEvent));

	  BOOST_FOREACH(const nbHoodSlot& nbs, sigNewNeighbourNotify)
	    getExtraLEParticleNeighbourhood(part, nbs.second);
	}

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
      
      //We now have the lowest cell coord, or corner of the cells to
      //check the contents of
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
		{
		  if (isUsedInScheduler)
		    Sim->ptrScheduler->addInteractionEvent(part, next);

		  BOOST_FOREACH(const nbHoodSlot& nbs, sigNewNeighbourNotify)
		    nbs.second(part, next);
		}
	      
	      nb += dim1pow;
	    }
	  
	  if (coords[dim1] + 2 >= cellCount[dim1]) nb += dim1pow * cellCount[dim1];
	  
	  nb += dim2pow - 3 * dim1pow;
	}
      
    }
   
  //Tell about the new locals
  BOOST_FOREACH(const size_t& lID, cells[endCell].locals)
    {
      if (isUsedInScheduler)
	Sim->ptrScheduler->addLocalEvent(part, lID);

      BOOST_FOREACH(const nbHoodSlot& nbs, sigNewLocalNotify)
	nbs.second(part, lID);
    }
      
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
    CVector<int> tmp = cells[oldCell].coords;
    CVector<int> tmp2 = cells[endCell].coords;
    
    std::cerr << "\nCGWall sysdt " 
	      << Sim->dSysTime / Sim->dynamics.units().unitTime()
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
CGCellsShearing::getParticleNeighbourhood(const Particle& part,
				   const nbHoodFunc& func) const
{
  CGCells::getParticleNeighbourhood(part, func);
  
  CVector<int> coords(cells[partCellData[part.getID()].cell].coords);
  if ((coords[1] == 0) || (coords[1] == cellCount[1] - 1))
    getExtraLEParticleNeighbourhood(part, func);

}

void 
CGCellsShearing::getExtraLEParticleNeighbourhood(const Particle& part,
						 const nbHoodFunc& func) const
{
  size_t cellID = partCellData[part.getID()].cell;
  CVector<int> coords(cells[cellID].coords);
  
#ifdef DYNAMO_DEBUG
  if ((coords[1] != 0) && (coords[1] != cellCount[1] - 1))
    M_throw() << "Shouldn't call this function unless the particle is at a border in the y dimension";
#endif 

  //Move to the bottom of x
  cellID -= coords[0];	
  coords[0] = 0;
  
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

