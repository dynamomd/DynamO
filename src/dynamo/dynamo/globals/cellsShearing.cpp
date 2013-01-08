/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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

#include <dynamo/globals/cellsShearing.hpp>
#include <dynamo/globals/globEvent.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/BC/LEBC.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
  GCellsShearing::GCellsShearing(dynamo::Simulation* nSim, 
				 const std::string& globalname):
    GCells(nSim, globalname)
  {
    name = "ShearingCells";
    dout << "Shearing Cells Loaded" << std::endl;
  }

  GCellsShearing::GCellsShearing(const magnet::xml::Node& XML, 
				 dynamo::Simulation* ptrSim):
    GCells(ptrSim, "ShearingCells")
  {
    operator<<(XML);
    dout << "Cells in shearing Loaded" << std::endl;
  }

  void 
  GCellsShearing::initialise(size_t nID)
  {
    ID=nID;

    if (!std::tr1::dynamic_pointer_cast<BCLeesEdwards>(Sim->BCs))
      derr << "You should not use the shearing neighbour list"
	   << " in a system without Lees Edwards BC's" << std::endl;

    if (overlink != 1) M_throw() << "Cannot shear with overlinking yet";

    reinitialise();
  }

  GlobalEvent 
  GCellsShearing::getEvent(const Particle& part) const
  {
#ifdef ISSS_DEBUG
    if (!Sim->dynamics->isUpToDate(part))
      M_throw() << "Particle is not up to date";
#endif

    //We do not inherit GCells get Event as the calcPosition thing done
    //for infinite systems is breaking it for shearing for some reason.
    return GlobalEvent(part,
		       Sim->dynamics->
		       getSquareCellCollision2
		       (part, calcPosition(partCellData[part.getID()]), 
			cellDimension)
		       - Sim->dynamics->getParticleDelay(part),
		       CELL, *this);
  }

  void 
  GCellsShearing::runEvent(Particle& part, const double) const
  {
    Sim->dynamics->updateParticle(part);

    size_t oldCell(partCellData[part.getID()]);
    magnet::math::MortonNumber<3> oldCellCoords(oldCell);
    Vector oldCellPosition(calcPosition(oldCellCoords));

    //Determine the cell transition direction, its saved
    int cellDirectionInt(Sim->dynamics->
			 getSquareCellCollision3(part, oldCellPosition, cellDimension));
  
    size_t cellDirection = abs(cellDirectionInt) - 1;

    magnet::math::MortonNumber<3> endCell = oldCellCoords; //The ID of the cell the particle enters

    if ((cellDirection == 1) &&
	(oldCellCoords[1] == ((cellDirectionInt < 0) ? 0 : (cellCount[1] - 1))))
      {
	//We're wrapping in the y direction, we have to compute
	//which cell its entering
	endCell[1] = (endCell[1].getRealValue() + cellCount[1] 
		      + ((cellDirectionInt < 0) ?  -1 : 1)) % cellCount[1];

	//Remove the old x contribution
	//Calculate the final x value
	//Time till transition, assumes the particle is up to date
	double dt = Sim->dynamics->getSquareCellCollision2(part, oldCellPosition, cellDimension);
     
	//Predict the position of the particle in the x dimension
	Sim->dynamics->advanceUpdateParticle(part, dt);
	Vector tmpPos = part.getPosition();
	//This rewinds the particle again
	Sim->dynamics->updateParticle(part);

	//Adding this extra half cell ensures we get into the next
	//simulation image, to calculate the position of the new cell
	tmpPos[1] += ((cellDirectionInt < 0) ? -0.5 : 0.5) * cellDimension[1];

	//Determine the x position (in cell coords) of the particle and
	//add it to the endCellID
	Sim->BCs->applyBC(tmpPos, dt);
      
	endCell[0] = getCellID(tmpPos)[0];

	removeFromCell(part.getID());
	addToCell(part.getID(), endCell.getMortonNum());
      
	//Get rid of the virtual event that is next, update is delayed till
	//after all events are added
	Sim->ptrScheduler->popNextEvent();

	//Check the entire neighbourhood, could check just the new
	//neighbours and the extra LE neighbourhood strip but its a lot
	//of code
	if (isUsedInScheduler)
	  {
	    BOOST_FOREACH(const size_t& id2, getParticleNeighbours(part))
	      {
		Sim->ptrScheduler->addInteractionEvent(part, id2);

		BOOST_FOREACH(const nbHoodSlot& nbs, sigNewNeighbourNotify)
		  nbs.second(part, id2);
	      }
	  }
      }
    else if ((cellDirection == 1) &&
	     (oldCellCoords[1] == ((cellDirectionInt < 0) ? 1 : (cellCount[1] - 2))))
      {
	//We're entering the boundary of the y direction
	//Calculate the end cell, no boundary wrap check required
	if (cellDirectionInt > 0)
	  endCell[cellDirection] = (endCell[cellDirection].getRealValue() + 1) % cellCount[cellDirection];
	else
	  endCell[cellDirection] = (endCell[cellDirection].getRealValue() 
				    + cellCount[cellDirection] - 1) % cellCount[cellDirection];

	removeFromCell(part.getID());
	addToCell(part.getID(), endCell.getMortonNum());
      
	//Get rid of the virtual event that is next, update is delayed till
	//after all events are added
	Sim->ptrScheduler->popNextEvent();
      
	//Check the extra LE neighbourhood strip
	if (isUsedInScheduler)
	  {
	    BOOST_FOREACH(const size_t& id2, getAdditionalLEParticleNeighbourhood(part))
	      {
		Sim->ptrScheduler->addInteractionEvent(part, id2);
		BOOST_FOREACH(const nbHoodSlot& nbs, sigNewNeighbourNotify)
		  nbs.second(part, id2);
	      }
	  }
      }
    else
      {
	//Here we follow the same procedure (except one more if statement) as the original cell list for new neighbours
	//The coordinates of the new center cell in the neighbourhood of the
	//particle
	magnet::math::MortonNumber<3> newNBCell(oldCell);
	if (cellDirectionInt > 0)
	  {
	    endCell[cellDirection] = (endCell[cellDirection].getRealValue() + 1) % cellCount[cellDirection];
	    newNBCell[cellDirection] = (endCell[cellDirection].getRealValue() + overlink) % cellCount[cellDirection];
	  }
	else
	  {
	    //We use the trick of adding cellCount to convert the
	    //subtraction to an addition, to prevent errors in the modulus
	    //of underflowing unsigned integers.
	    endCell[cellDirection] = (endCell[cellDirection].getRealValue() 
				      + cellCount[cellDirection] - 1) % cellCount[cellDirection];
	    newNBCell[cellDirection] = (endCell[cellDirection].getRealValue() 
					+ cellCount[cellDirection] - overlink) % cellCount[cellDirection];
	  }
    
	removeFromCell(part.getID());
	addToCell(part.getID(), endCell.getMortonNum());

	//Get rid of the virtual event we're running, an updated event is
	//pushed after all other events are added
	Sim->ptrScheduler->popNextEvent();


	if ((cellDirection == 2) &&
	    ((oldCellCoords[1] == 0) || (oldCellCoords[1] == cellCount[1] -1)))
	  {
	    //We're at the boundary moving in the z direction, we must
	    //add the new LE strips as neighbours	
	    //We just check the entire Extra LE neighbourhood
	    BOOST_FOREACH(const size_t& id2, getAdditionalLEParticleNeighbourhood(part))
	      {
		BOOST_FOREACH(const nbHoodSlot& nbs, sigNewNeighbourNotify)
		  nbs.second(part, id2);
	      }
	  }

	//Particle has just arrived into a new cell warn the scheduler about
	//its new neighbours so it can add them to the heap
	//Holds the displacement in each dimension, the unit is cells!

	//These are the two dimensions to walk in
	size_t dim1 = (cellDirection + 1) % 3,
	  dim2 = (cellDirection + 2) % 3;

	newNBCell[dim1] += cellCount[dim1] - overlink;
	newNBCell[dim2] += cellCount[dim2] - overlink;
  
	size_t walkLength = 2 * overlink + 1;

	const magnet::math::DilatedInteger<3> saved_coord(newNBCell[dim1]);

	//We now have the lowest cell coord, or corner of the cells to update
	for (size_t iDim(0); iDim < walkLength; ++iDim)
	  {
	    newNBCell[dim2] %= cellCount[dim2];

	    for (size_t jDim(0); jDim < walkLength; ++jDim)
	      {
		newNBCell[dim1] %= cellCount[dim1];
  
		BOOST_FOREACH(const size_t& next, list[newNBCell.getMortonNum()])
		  BOOST_FOREACH(const nbHoodSlot& nbs, sigNewNeighbourNotify)
		    nbs.second(part, next);
	  
		++newNBCell[dim1];
	      }

	    newNBCell[dim1] = saved_coord; 
	    ++newNBCell[dim2];
	  }
      }
    
    //Push the next virtual event, this is the reason the scheduler
    //doesn't need a second callback
    Sim->ptrScheduler->pushEvent(part, getEvent(part));
    Sim->ptrScheduler->sort(part);

    BOOST_FOREACH(const nbHoodSlot& nbs, sigCellChangeNotify)
      nbs.second(part, oldCell);
  
#ifdef DYNAMO_WallCollDebug
    {
      magnet::math::MortonNumber<3> newNBCellv(oldCell);
      magnet::math::MortonNumber<3> endCellv(endCell);
    
      dout << "CellEvent: sysdt " 
	   << Sim->systemTime / Sim->units.unitTime()
	   << " ID "
	   << part.getID()
	   << "  from <" 
	   << newNBCellv[0].getRealValue() << "," << newNBCellv[1].getRealValue() 
	   << "," << newNBCellv[2].getRealValue()
	   << "> to <" 
	   << endCellv[0].getRealValue() << "," << endCellv[1].getRealValue()
	   << "," << endCellv[2].getRealValue() << ">"
	   << std::endl;
    }
#endif
  }

  IDRangeList
  GCellsShearing::getParticleNeighbours(const Particle& part) const
  {
    return getParticleNeighbours(magnet::math::MortonNumber<3>(partCellData[part.getID()]));
  }

  IDRangeList
  GCellsShearing::getParticleNeighbours(const Vector& vec) const
  {
    return getParticleNeighbours(magnet::math::MortonNumber<3>(getCellID(vec)));
  }

  IDRangeList
  GCellsShearing::getParticleNeighbours(const magnet::math::MortonNumber<3>& cellCoords) const
  {
    IDRangeList retval(GCells::getParticleNeighbours(cellCoords));

    if ((cellCoords[1] == 0) || (cellCoords[1] == dilatedCellMax[1]))
      {
	const std::vector<size_t>& extra(getAdditionalLEParticleNeighbourhood(cellCoords));
	retval.getContainer().insert(retval.getContainer().end(), extra.begin(), extra.end());
      }
    return retval;
  }
  
  std::vector<size_t>
  GCellsShearing::getAdditionalLEParticleNeighbourhood(const Particle& part) const
  {
    return getAdditionalLEParticleNeighbourhood(magnet::math::MortonNumber<3>(partCellData[part.getID()]));
  }

  std::vector<size_t>
  GCellsShearing::getAdditionalLEParticleNeighbourhood(magnet::math::MortonNumber<3> cellCoords) const
  {  
#ifdef DYNAMO_DEBUG
    if ((cellCoords[1] != 0) && (cellCoords[1] != dilatedCellMax[1]))
      M_throw() << "Shouldn't call this function unless the particle is at a border in the y dimension";
#endif 

    //Move to the bottom of x
    cellCoords[0] = 0;
    //Get the correct y-side (its the opposite to the particles current side)
    cellCoords[1] = (cellCoords[1] > 0) ? 0 : dilatedCellMax[1];  
    ////Move te overlink across
    cellCoords[2] = (cellCoords[2].getRealValue() + cellCount[2] - overlink) % cellCount[2];

    std::vector<size_t> retval;
    retval.reserve(32);
    for (size_t i(0); i < 2 * overlink + 1; ++i)
      {
	cellCoords[2] %= cellCount[2];

	for (size_t j(0); j < cellCount[0]; ++j)
	  {
	    const std::vector<size_t>& nbs(list[cellCoords.getMortonNum()]);
	    retval.insert(retval.end(), nbs.begin(), nbs.end());

	    ++cellCoords[0];
	  }
	++cellCoords[2];
	cellCoords[0] = 0;
      }
    return retval;
  }
}
