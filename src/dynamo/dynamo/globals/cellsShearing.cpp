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
    setOutputPrefix("ShearingCells");
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
    Global::initialise(nID);

    if (!std::dynamic_pointer_cast<BCLeesEdwards>(Sim->BCs))
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
    return GlobalEvent(part, Sim->dynamics->getSquareCellCollision2(part, calcPosition(_cellData.getCellID(part.getID())), cellDimension) - Sim->dynamics->getParticleDelay(part), CELL, *this);
  }

  void 
  GCellsShearing::runEvent(Particle& part, const double) const
  {
    Sim->dynamics->updateParticle(part);

    //Get rid of the virtual event that is next, update is delayed
    //till after all events are added
    Sim->ptrScheduler->popNextEvent();

    const size_t oldCellIndex(_cellData.getCellID(part.getID()));
    const auto oldCellCoord = _ordering.toCoord(oldCellIndex);
    const Vector oldCellPosition(calcPosition(oldCellCoord));
    const int cellDirectionInt(Sim->dynamics->getSquareCellCollision3(part, oldCellPosition, cellDimension));
    const size_t cellDirection = abs(cellDirectionInt) - 1;

    auto newCellCoord = oldCellCoord;
    newCellCoord[cellDirection] += _ordering.getDimensions()[cellDirection] + ((cellDirectionInt > 0) ? 1 : -1);
    newCellCoord[cellDirection] %= _ordering.getDimensions()[cellDirection];

    if ((cellDirection == 1) && (oldCellCoord[1] == ((cellDirectionInt < 0) ? 0 : (_ordering.getDimensions()[1] - 1))))
      {
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
      
	newCellCoord[0] = getCellCoords(tmpPos)[0];

	_cellData.moveTo(oldCellIndex, _ordering.toIndex(newCellCoord), part.getID());
      
	//Check the entire neighbourhood, could check just the new
	//neighbours and the extra LE neighbourhood strip but its a lot
	//of code
	std::vector<size_t> neighbours;
	GCells::getParticleNeighbours(part, neighbours);
	for (const size_t& id2 : neighbours)
	  _sigNewNeighbour(part, id2);
      }
    else if ((cellDirection == 1) && (oldCellCoord[1] == ((cellDirectionInt < 0) ? 1 : (_ordering.getDimensions()[1] - 2))))
      {
	//We're entering the boundary of the y direction
	//Calculate the end cell, no boundary wrap check required
	_cellData.moveTo(oldCellIndex, _ordering.toIndex(newCellCoord), part.getID());
            
	//Check the extra LE neighbourhood strip
	std::vector<size_t> nbs;
	getAdditionalLEParticleNeighbourhood(part, nbs);
	for (const size_t& id2 : nbs)
	  {
	    Sim->ptrScheduler->addInteractionEvent(part, id2);
	    _sigNewNeighbour(part, id2);
	  }
      }
    else
      {
	_cellData.moveTo(oldCellIndex, _ordering.toIndex(newCellCoord), part.getID());

	auto newNBCellCoord = newCellCoord;
	newNBCellCoord[cellDirection] += _ordering.getDimensions()[cellDirection] + ((cellDirectionInt > 0) ? 1 : -1);
	newNBCellCoord[cellDirection] %= _ordering.getDimensions()[cellDirection];

	if ((cellDirection == 2) && ((oldCellCoord[1] == 0) || (oldCellCoord[1] == _ordering.getDimensions()[1] - 1)))
	  {
	    //We're at the boundary moving in the z direction, we must
	    //add the new LE strips as neighbours	
	    //We just check the entire Extra LE neighbourhood
	    std::vector<size_t> nbs;
	    getAdditionalLEParticleNeighbourhood(part, nbs);
	    for (const size_t& id2 : nbs) _sigNewNeighbour(part, id2);
	  }

	//Particle has just arrived into a new cell warn the scheduler about
	//its new neighbours so it can add them to the heap
	//Holds the displacement in each dimension, the unit is cells!

	//These are the two dimensions to walk in
	std::array<size_t, 3> steps = {{overlink, overlink, overlink}};
	steps[cellDirection] = 0;
	
	for (auto cellIndex : _ordering.getSurroundingIndices(newNBCellCoord, steps))
	  for (const size_t& next : _cellData.getCellContents(cellIndex))
	    _sigNewNeighbour(part, next);
      }
    
    //Push the next virtual event, this is the reason the scheduler
    //doesn't need a second callback
    Sim->ptrScheduler->pushEvent(part, getEvent(part));
    Sim->ptrScheduler->sort(part);

    _sigCellChange(part, oldCellIndex);
  }

  void
  GCellsShearing::getParticleNeighbours(const std::array<size_t, 3>& cellCoords, std::vector<size_t>& retlist) const
  {
    GCells::getParticleNeighbours(cellCoords, retlist);
    if ((cellCoords[1] == 0) || (cellCoords[1] == (_ordering.getDimensions()[1] - 1)))
      getAdditionalLEParticleNeighbourhood(cellCoords, retlist);
  }
  
  void
  GCellsShearing::getAdditionalLEParticleNeighbourhood(const Particle& part, std::vector<size_t>& retlist) const {
    return getAdditionalLEParticleNeighbourhood(_ordering.toCoord(_cellData.getCellID(part.getID())), retlist);
  }

  void
  GCellsShearing::getAdditionalLEParticleNeighbourhood(std::array<size_t, 3> cellCoords, std::vector<size_t>& retlist) const
  {  
#ifdef DYNAMO_DEBUG
    if ((cellCoords[1] != 0) && (cellCoords[1] != (_ordering.getDimensions()[1] - 1)))
      M_throw() << "Shouldn't call this function unless the particle is at a border in the y dimension";
#endif
    std::array<size_t, 3> start = {{0, (cellCoords[1] > 0) ? 0 : _ordering.getDimensions()[1] - 1, cellCoords[2]}};
    std::array<size_t, 3> steps = {{_ordering.getDimensions()[0], 0, overlink}};
    //These are the two dimensions to walk in
    for (auto cellIndex : _ordering.getSurroundingIndices(start, steps))
      {
	const auto neighbours = _cellData.getCellContents(cellIndex);
	retlist.insert(retlist.end(), neighbours.begin(), neighbours.end());
      }
  }
}
