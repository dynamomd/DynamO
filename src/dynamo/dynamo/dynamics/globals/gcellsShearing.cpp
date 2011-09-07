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

#include <dynamo/dynamics/globals/gcellsShearing.hpp>
#include <dynamo/dynamics/globals/globEvent.hpp>
#include <dynamo/dynamics/NparticleEventData.hpp>
#include <dynamo/dynamics/liouvillean/liouvillean.hpp>
#include <dynamo/dynamics/units/units.hpp>
#include <dynamo/dynamics/ranges/1RAll.hpp>
#include <dynamo/dynamics/ranges/1RNone.hpp>
#include <dynamo/dynamics/ranges/2RSingle.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/dynamics/locals/local.hpp>
#include <dynamo/dynamics/BC/LEBC.hpp>
#include <dynamo/dynamics/liouvillean/NewtonianGravityL.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
  GCellsShearing::GCellsShearing(dynamo::SimData* nSim, 
				 const std::string& globalname):
    GCells(nSim, globalname)
  {
    name = "ShearingCells";
    dout << "Shearing Cells Loaded" << std::endl;
  }

  GCellsShearing::GCellsShearing(const magnet::xml::Node& XML, 
				 dynamo::SimData* ptrSim):
    GCells(ptrSim, "Unknown")
  {
    operator<<(XML);
    name = "ShearingCells";

    dout << "Cells in shearing Loaded" << std::endl;
  }

  void 
  GCellsShearing::initialise(size_t nID)
  {
    ID=nID;
   
    if (!(Sim->dynamics.BCTypeTest<BCLeesEdwards>()))
      derr << "You should not use the shearing neighbour list"
	   << " in a system without Lees Edwards BC's" << std::endl;

    if (overlink != 1) M_throw() << "Cannot shear with overlinking yet";

    reinitialise();
  }

  void
  GCellsShearing::outputXML(magnet::xml::XmlStream& XML) const
  { GCells::outputXML(XML, "ShearingCells"); }

  GlobalEvent 
  GCellsShearing::getEvent(const Particle& part) const
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
		       (part, calcPosition(partCellData[part.getID()].cell), 
			cellDimension)
		       - Sim->dynamics.getLiouvillean().getParticleDelay(part),
		       CELL, *this);
  }

  void 
  GCellsShearing::runEvent(const Particle& part, const double) const
  {
    Sim->dynamics.getLiouvillean().updateParticle(part);

    size_t oldCell(partCellData[part.getID()].cell);
    magnet::math::MortonNumber<3> oldCellCoords(oldCell);
    Vector oldCellPosition(calcPosition(oldCellCoords));

    //Determine the cell transition direction, its saved
    int cellDirectionInt(Sim->dynamics.getLiouvillean().
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
	double dt = Sim->dynamics.getLiouvillean().getSquareCellCollision2(part, oldCellPosition, cellDimension);
     
	//Predict the position of the particle in the x dimension
	Sim->dynamics.getLiouvillean().advanceUpdateParticle(part, dt);
	Vector tmpPos = part.getPosition();
	//This rewinds the particle again
	Sim->dynamics.getLiouvillean().updateParticle(part);

	//Adding this extra half cell ensures we get into the next
	//simulation image, to calculate the position of the new cell
	tmpPos[1] += ((cellDirectionInt < 0) ? -0.5 : 0.5) * cellDimension[1];

	//Determine the x position (in cell coords) of the particle and
	//add it to the endCellID
	Sim->dynamics.BCs().applyBC(tmpPos, dt);
      
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
	  getParticleNeighbourhood(part, magnet::function::MakeDelegate(&(*Sim->ptrScheduler), 
									&Scheduler::addInteractionEvent));
      
	BOOST_FOREACH(const nbHoodSlot& nbs, sigNewNeighbourNotify)
	  getParticleNeighbourhood(part, nbs.second);
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
	  getExtraLEParticleNeighbourhood(part, magnet::function::MakeDelegate(&(*Sim->ptrScheduler),
									       &Scheduler::addInteractionEvent));
      
	BOOST_FOREACH(const nbHoodSlot& nbs, sigNewNeighbourNotify)
	  getExtraLEParticleNeighbourhood(part, nbs.second);
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
	  //We're at the boundary moving in the z direction, we must
	  //add the new LE strips as neighbours	
	  //We just check the entire Extra LE neighbourhood
	  {
	    if (isUsedInScheduler)
	      getExtraLEParticleNeighbourhood(part, magnet::function::MakeDelegate(&(*Sim->ptrScheduler),
										   &Scheduler::addInteractionEvent));
	  
	    BOOST_FOREACH(const nbHoodSlot& nbs, sigNewNeighbourNotify)
	      getExtraLEParticleNeighbourhood(part, nbs.second);
	  }

	//Particle has just arrived into a new cell warn the scheduler about
	//its new neighbours so it can add them to the heap
	//Holds the displacement in each dimension, the unit is cells!

	//These are the two dimensions to walk in
	size_t dim1 = (cellDirection + 1) % 3,
	  dim2 = (cellDirection + 2) % 3;

	newNBCell[dim1] += cellCount[dim1] - overlink;
	newNBCell[dim2] += cellCount[dim1] - overlink;
  
	size_t walkLength = 2 * overlink + 1;

	const magnet::math::DilatedInteger<3> saved_coord(newNBCell[dim1]);

	//We now have the lowest cell coord, or corner of the cells to update
	for (size_t iDim(0); iDim < walkLength; ++iDim)
	  {
	    newNBCell[dim2] %= cellCount[dim2];

	    for (size_t jDim(0); jDim < walkLength; ++jDim)
	      {
		newNBCell[dim1] %= cellCount[dim1];
  
		for (int next = list[newNBCell.getMortonNum()]; next >= 0; 
		     next = partCellData[next].next)
		  {
		    if (isUsedInScheduler)
		      Sim->ptrScheduler->addInteractionEvent(part, next);

		    BOOST_FOREACH(const nbHoodSlot& nbs, sigNewNeighbourNotify)
		      nbs.second(part, next);
		  }
	  
		++newNBCell[dim1];
	      }

	    newNBCell[dim1] = saved_coord; 
	    ++newNBCell[dim2];
	  }
      }
	   
    BOOST_FOREACH(const size_t& lID, cells[endCell.getMortonNum()])
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
  
#ifdef DYNAMO_WallCollDebug
    {
      magnet::math::MortonNumber<3> newNBCellv(oldCell);
      magnet::math::MortonNumber<3> endCellv(endCell);
    
      std::cerr << "\nCGWall sysdt " 
		<< Sim->dSysTime / Sim->dynamics.units().unitTime()
		<< "  WALL ID "
		<< part.getID()
		<< "  from <" 
		<< newNBCellv.data[0].getRealVal() << "," << newNBCellv.data[1].getRealVal() 
		<< "," << newNBCellv.data[2].getRealVal()
		<< "> to <" 
		<< endCellv.data[0].getRealVal() << "," << endCellv.data[1].getRealVal() 
		<< "," << endCellv.data[2].getRealVal();
    }
#endif
  }

  void 
  GCellsShearing::getParticleNeighbourhood(const Particle& part,
					   const nbHoodFunc& func) const
  {
    GCells::getParticleNeighbourhood(part, func);
  
    size_t cell(partCellData[part.getID()].cell);
    magnet::math::MortonNumber<3> cellCoords(cell);

    if ((cellCoords[1] == 0) || (cellCoords[1] == dilatedCellMax[1]))
      getExtraLEParticleNeighbourhood(part, func);
  }

  void 
  GCellsShearing::getParticleNeighbourhood(const Vector& vec,
					   const nbHoodFunc2& func) const
  {
    GCells::getParticleNeighbourhood(vec, func);
  
    magnet::math::MortonNumber<3> cellCoords = getCellID(vec);

    if ((cellCoords[1] == 0) || (cellCoords[1] == dilatedCellMax[1]))
      {
	//Move to the bottom of x
	cellCoords[0] = 0;
	//Get the correct y-side (its the opposite to the particles current side)
	cellCoords[1] = (cellCoords[1] > 0) ? 0 : dilatedCellMax[1];  
	////Move a single cell down in Z
	cellCoords[2] = (cellCoords[2].getRealValue() + cellCount[2] - 1) % cellCount[2];
	
	for (size_t i(0); i < 3; ++i)
	  {
	    cellCoords[2] %= cellCount[2];
	    
	    for (size_t j(0); j < cellCount[0]; ++j)
	      {
		for (int next = list[cellCoords.getMortonNum()]; next >= 0; 
		     next = partCellData[next].next)
		  func(next);
		
		++cellCoords[0];
	      }
	    ++cellCoords[2];
	    cellCoords[0] = 0;
	  }
      }
  }


  void 
  GCellsShearing::getExtraLEParticleNeighbourhood(const Particle& part,
						  const nbHoodFunc& func) const
  {
    size_t cell(partCellData[part.getID()].cell);
    magnet::math::MortonNumber<3> cellCoords(cell);
  
#ifdef DYNAMO_DEBUG
    if ((cellCoords[1] != 0) && (cellCoords[1] != dilatedCellMax[1]))
      M_throw() << "Shouldn't call this function unless the particle is at a border in the y dimension";
#endif 

    //Move to the bottom of x
    cellCoords[0] = 0;
    //Get the correct y-side (its the opposite to the particles current side)
    cellCoords[1] = (cellCoords[1] > 0) ? 0 : dilatedCellMax[1];  
    ////Move a single cell down in Z
    cellCoords[2] = (cellCoords[2].getRealValue() + cellCount[2] - 1) % cellCount[2];

    for (size_t i(0); i < 3; ++i)
      {
	cellCoords[2] %= cellCount[2];

	for (size_t j(0); j < cellCount[0]; ++j)
	  {
	    for (int next = list[cellCoords.getMortonNum()]; next >= 0; 
		 next = partCellData[next].next)
	      func(part, next);

	    ++cellCoords[0];
	  }
	++cellCoords[2];
	cellCoords[0] = 0;
      }
  }
}
