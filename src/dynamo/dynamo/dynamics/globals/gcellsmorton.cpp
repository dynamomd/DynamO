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

#include "gcellsmorton.hpp"
#include "globEvent.hpp"
#include "../NparticleEventData.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../units/units.hpp"
#include "../ranges/1RAll.hpp"
#include "../../schedulers/scheduler.hpp"
#include "../locals/local.hpp"
#include "../BC/LEBC.hpp"
#include "../liouvillean/NewtonianGravityL.hpp"
#include <magnet/math/ctime_pow.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/static_assert.hpp>
#include <cstdio>


CGCellsMorton::CGCellsMorton(dynamo::SimData* nSim, const std::string& name):
  CGNeighbourList(nSim, "MortonCellNeighbourList"),
  cellCount(0),
  cellDimension(1,1,1),
  _oversizeCells(1.0),
  NCells(0),
  overlink(1)
{
  globName = name;
  dout << "Cells Loaded" << std::endl;
}

CGCellsMorton::CGCellsMorton(const magnet::xml::Node& XML, dynamo::SimData* ptrSim):
  CGNeighbourList(ptrSim, "MortonCellNeighbourList"),
  cellCount(0),
  cellDimension(1,1,1),
  _oversizeCells(1.0),
  NCells(0),
  overlink(1)
{
  operator<<(XML);

  dout << "Cells Loaded" << std::endl;
}

CGCellsMorton::CGCellsMorton(dynamo::SimData* ptrSim, const char* nom, void*):
  CGNeighbourList(ptrSim, nom),
  cellCount(0),
  cellDimension(1,1,1),
  _oversizeCells(1.0),
  NCells(0),
  overlink(1)
{}

void 
CGCellsMorton::operator<<(const magnet::xml::Node& XML)
{
  try {
    //If you add anything here then it needs to go in gListAndCells.cpp too
    if (XML.hasAttribute("OverLink"))
      overlink = XML.getAttribute("OverLink").as<size_t>();

    if (XML.hasAttribute("Oversize"))
      _oversizeCells = XML.getAttribute("Oversize").as<double>();

    if (_oversizeCells < 1.0)
      M_throw() << "You must specify an Oversize greater than 1.0, otherwise your cells are too small!";
    
    globName = XML.getAttribute("Name");	
  }
  catch(...)
    {
      M_throw() << "Error loading CGCellsMorton";
    }
}

GlobalEvent 
CGCellsMorton::getEvent(const Particle& part) const
{
#ifdef ISSS_DEBUG
  if (!Sim->dynamics.getLiouvillean().isUpToDate(part))
    M_throw() << "Particle is not up to date";
#endif

  //This 
  //Sim->dynamics.getLiouvillean().updateParticle(part);
  //is not required as we compensate for the delay using 
  //Sim->dynamics.getLiouvillean().getParticleDelay(part)
  
  return GlobalEvent(part,
		    Sim->dynamics.getLiouvillean().
		    getSquareCellCollision2
		     (part, calcPosition(partCellData[part.getID()].cell,
					 part), cellDimension)
		    -Sim->dynamics.getLiouvillean().getParticleDelay(part)
		    ,
		    CELL, *this);
}

void
CGCellsMorton::runEvent(const Particle& part, const double) const
{
  //Despite the system not being streamed this must be done.  This is
  //because the scheduler and all interactions, locals and systems
  //expect the particle to be up to date.
  Sim->dynamics.getLiouvillean().updateParticle(part);

  const size_t oldCell(partCellData[part.getID()].cell);
  size_t endCell;

  //Determine the cell transition direction, its saved
  int cellDirectionInt(Sim->dynamics.getLiouvillean().
		       getSquareCellCollision3
		       (part, calcPosition(oldCell, part), cellDimension));
  
  size_t cellDirection = abs(cellDirectionInt) - 1;

  magnet::math::DilatedVector inCell(oldCell);

  {
    magnet::math::DilatedVector dendCell(inCell);
    
    if (cellDirectionInt > 0)
      {
	++dendCell.data[cellDirection];
	inCell.data[cellDirection] = dendCell.data[cellDirection] + dilatedOverlink;	

	if (dendCell.data[cellDirection] > dilatedCellMax[cellDirection])
	    dendCell.data[cellDirection] = --dendCell.data[cellDirection]
	      - dilatedCellMax[cellDirection];

	if (inCell.data[cellDirection] > dilatedCellMax[cellDirection])
	  inCell.data[cellDirection] = --inCell.data[cellDirection]
	    - dilatedCellMax[cellDirection];
      }
    else
      {
	--dendCell.data[cellDirection];
	inCell.data[cellDirection] = dendCell.data[cellDirection] - dilatedOverlink;

	if (dendCell.data[cellDirection] > dilatedCellMax[cellDirection])
	  dendCell.data[cellDirection] = dendCell.data[cellDirection]
	    - (std::numeric_limits<magnet::math::DilatedInteger>::max() - dilatedCellMax[cellDirection]);

	if (inCell.data[cellDirection] > dilatedCellMax[cellDirection])
	  inCell.data[cellDirection] = inCell.data[cellDirection]
	    - (std::numeric_limits<magnet::math::DilatedInteger>::max() - dilatedCellMax[cellDirection]);
      }

    endCell = dendCell.getMortonNum();
  }
    
  removeFromCell(part.getID());
  addToCell(part.getID(), endCell);

  //Get rid of the virtual event that is next, update is delayed till
  //after all events are added
  Sim->ptrScheduler->popNextEvent();

  //Particle has just arrived into a new cell warn the scheduler about
  //its new neighbours so it can add them to the heap
  //Holds the displacement in each dimension, the unit is cells!
  BOOST_STATIC_ASSERT(NDIM==3);

  //These are the two dimensions to walk in
  size_t dim1 = cellDirection + 1 - 3 * (cellDirection > 1),
    dim2 = cellDirection + 2 - 3 * (cellDirection > 0);

  inCell.data[dim1] = inCell.data[dim1] - dilatedOverlink;
  inCell.data[dim2] = inCell.data[dim2] - dilatedOverlink;
  
  //Test if the data has looped around
  if (inCell.data[dim1] > dilatedCellMax[dim1])
    inCell.data[dim1] = inCell.data[dim1] - (std::numeric_limits<magnet::math::DilatedInteger>::max() - dilatedCellMax[dim1]);

  if (inCell.data[dim2] > dilatedCellMax[dim2]) 
    inCell.data[dim2] = inCell.data[dim2] - (std::numeric_limits<magnet::math::DilatedInteger>::max() - dilatedCellMax[dim2]);

  int walkLength = 2 * overlink + 1;

  const magnet::math::DilatedInteger saved_coord(inCell.data[dim1]);

  //We now have the lowest cell coord, or corner of the cells to update
  for (int iDim(0); iDim < walkLength; ++iDim)
    {
      if (inCell.data[dim2] > dilatedCellMax[dim2])
	inCell.data[dim2].zero();

      for (int jDim(0); jDim < walkLength; ++jDim)
	{
	  if (inCell.data[dim1] > dilatedCellMax[dim1])
	    inCell.data[dim1].zero();
  
	  for (int next = list[inCell.getMortonNum()]; next >= 0; 
	       next = partCellData[next].next)
	    {
	      if (isUsedInScheduler)
		Sim->ptrScheduler->addInteractionEvent(part, next);

	      BOOST_FOREACH(const nbHoodSlot& nbs, sigNewNeighbourNotify)
		nbs.second(part, next);
	    }
	  
	  
	  ++inCell.data[dim1];
	}

      inCell.data[dim1] = saved_coord;
      
      ++inCell.data[dim2];
    }

  //Tell about the new locals
  BOOST_FOREACH(const size_t& lID, cells[endCell])
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
    magnet::math::DilatedVector inCellv(oldCell);
    magnet::math::DilatedVector endCellv(endCell);
    
    std::cerr << "\nCGWall sysdt " 
	      << Sim->dSysTime / Sim->dynamics.units().unitTime()
	      << "  WALL ID "
	      << part.getID()
	      << "  from <" 
	      << inCellv.data[0].getRealVal() << "," << inCellv.data[1].getRealVal() 
	      << "," << inCellv.data[2].getRealVal()
	      << "> to <" 
	      << endCellv.data[0].getRealVal() << "," << endCellv.data[1].getRealVal() 
	      << "," << endCellv.data[2].getRealVal();
  }
#endif
}

void 
CGCellsMorton::initialise(size_t nID)
{
  ID=nID;
  reinitialise(getMaxInteractionLength());
}

void
CGCellsMorton::reinitialise(const double& maxdiam)
{
  dout << "Reinitialising on collision " << Sim->eventCount << std::endl;

  //Create the cells
  addCells(_oversizeCells * (maxdiam * (1.0 + 10 * std::numeric_limits<double>::epsilon())) / overlink);

  addLocalEvents();

  BOOST_FOREACH(const initSlot& nbs, sigReInitNotify)
    nbs.second();

  if (isUsedInScheduler)
    Sim->ptrScheduler->initialise();
}

void
CGCellsMorton::outputXML(magnet::xml::XmlStream& XML) const
{
  //If you add anything here it also needs to go in gListAndCells.cpp too
  XML << magnet::xml::tag("Global")
      << magnet::xml::attr("Type") << "CellsMorton"
      << magnet::xml::attr("Name") << globName;

  if (overlink > 1)   XML << magnet::xml::attr("OverLink") << overlink;
  if (_oversizeCells != 1.0) XML << magnet::xml::attr("Oversize") << _oversizeCells;
  
  XML << magnet::xml::endtag("Global");
}

void
CGCellsMorton::addCells(double maxdiam)
{
  cells.clear();
  partCellData.resize(Sim->N); //Location data for particles

  NCells = 1;
  cellCount = CVector<size_t>(0);

  if ((Sim->primaryCellSize[0] != Sim->primaryCellSize[1]) 
      || (Sim->primaryCellSize[0] != Sim->primaryCellSize[2])) 
    M_throw() << "This cellular neighbor list does not work unless the primary cell is square";

  for (size_t iDim = 0; iDim < NDIM; iDim++)
    {
      cellCount[iDim] = int(Sim->primaryCellSize[iDim] / (maxdiam * (1.0 + 10 * std::numeric_limits<double>::epsilon())));
      
      if (cellCount[iDim] < 3)
	M_throw() << "Not enough cells in " << char('x'+iDim) << " dimension, need 3+";
      
      NCells *= cellCount[iDim];
      
      dilatedCellMax[iDim] = cellCount[iDim] - 1;
      cellLatticeWidth[iDim] = Sim->primaryCellSize[iDim] / cellCount[iDim];
      cellDimension[iDim] = cellLatticeWidth[iDim] + (cellLatticeWidth[iDim] - maxdiam) 
	* lambda;
      cellOffset[iDim] = -(cellLatticeWidth[iDim] - maxdiam) * lambda * 0.5;
    }

  dilatedOverlink = overlink;
  
  dout << "Cells <x,y,z>  " << cellCount[0] << ","
	   << cellCount[1] << "," << cellCount[2] << std::endl;

  dout << "Cell Offset <x,y,z>  "
           << cellOffset[0] / Sim->dynamics.units().unitLength() << ","
	   << cellOffset[1] / Sim->dynamics.units().unitLength() << ","
	   << cellOffset[2] / Sim->dynamics.units().unitLength() << std::endl;

  dout << "Cells Dimension <x,y,z>  " 
	   << cellDimension[0] / Sim->dynamics.units().unitLength()
	   << ","
	   << cellDimension[1] / Sim->dynamics.units().unitLength()
	   << "," 
	   << cellDimension[2] / Sim->dynamics.units().unitLength() << std::endl;

  dout << "Lattice spacing <x,y,z>  " 
	   << cellLatticeWidth[0] / Sim->dynamics.units().unitLength()
	   << ","
	   << cellLatticeWidth[1] / Sim->dynamics.units().unitLength()
	   << "," 
	   << cellLatticeWidth[2] / Sim->dynamics.units().unitLength() << std::endl;

  //Find the required size of the morton array
  size_t sizeReq(1);

  for (int i(0); i < magnet::math::ctime_pow<2,magnet::math::DilatedInteger::digits>::result; ++i)
    {
      sizeReq *= 2*2*2;
      if (sizeReq >= NCells) break;
    }

  cells.resize(sizeReq); //Empty Cells created!
  list.resize(sizeReq); //Empty Cells created!

  dout << "Vector Size <N>  " << sizeReq << std::endl;
  
  for (size_t iDim = 0; iDim < cellCount[0]; ++iDim)
    for (size_t jDim = 0; jDim < cellCount[1]; ++jDim)
      for (size_t kDim = 0; kDim < cellCount[2]; ++kDim)
	{
	  magnet::math::DilatedVector coords(iDim, jDim, kDim);
	  size_t id = coords.getMortonNum();
	  list[id] = -1;
	}

  //Add the particles section
  //Required so particles find the right owning cell
  Sim->dynamics.getLiouvillean().updateAllParticles(); 

  ////initialise the data structures
  BOOST_FOREACH(const Particle& part, Sim->particleList)
    addToCell(part.getID(), getCellID(part.getPosition()).getMortonNum());
}

void 
CGCellsMorton::addLocalEvents()
{
  for (size_t iDim = 0; iDim < cellCount[0]; ++iDim)
    for (size_t jDim = 0; jDim < cellCount[1]; ++jDim)
      for (size_t kDim = 0; kDim < cellCount[2]; ++kDim)
	{
	  magnet::math::DilatedVector coords(iDim, jDim, kDim);
	  size_t id = coords.getMortonNum();
	  cells[id].clear();
	  Vector pos = calcPosition(coords);
	  
	  //We make the box slightly larger to ensure objects on the boundary are included
	  BOOST_FOREACH(const magnet::ClonePtr<Local>& local, Sim->dynamics.getLocals())
	    if (local->isInCell(pos - 0.0001 * cellDimension, 1.0002 * cellDimension))
	      cells[id].push_back(local->getID());
	}
}

magnet::math::DilatedVector
CGCellsMorton::getCellID(const CVector<int>& coordsold) const
{
  //PBC for vectors
  CVector<int> coords(coordsold);

  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    {
      coords[iDim] %= cellCount[iDim];
      if (coords[iDim] < 0) coords[iDim] += cellCount[iDim];
    }
  
  return magnet::math::DilatedVector(coords[0],coords[1],coords[2]);
}

magnet::math::DilatedVector
CGCellsMorton::getCellID(Vector pos) const
{
  Sim->dynamics.BCs().applyBC(pos);
  CVector<int> temp;
  
  for (size_t iDim = 0; iDim < NDIM; iDim++)
    temp[iDim] = std::floor((pos[iDim] + 0.5 * Sim->primaryCellSize[iDim] - cellOffset[iDim])
			    / cellLatticeWidth[iDim]);
  
  return getCellID(temp);
}


void 
CGCellsMorton::getParticleNeighbourhood(const Particle& part,
					const nbHoodFunc& func) const
{
  BOOST_STATIC_ASSERT(NDIM==3);

  const magnet::math::DilatedVector center_coords(partCellData[part.getID()].cell);
  magnet::math::DilatedVector coords(center_coords);

  for (size_t iDim(0); iDim < NDIM; ++iDim)
    {
      coords.data[iDim] -= dilatedOverlink;
      if (coords.data[iDim] > dilatedCellMax[iDim])
	coords.data[iDim] -= 
	  std::numeric_limits<magnet::math::DilatedInteger>::max() - dilatedCellMax[iDim];
    }

  const magnet::math::DilatedVector zero_coords(coords);

  coords = center_coords;
  for (size_t iDim(0); iDim < NDIM; ++iDim)
    {
      coords.data[iDim] += dilatedOverlink + 1;
      if (coords.data[iDim] > dilatedCellMax[iDim])
	coords.data[iDim] -= dilatedCellMax[iDim] + 1;
    }

  const magnet::math::DilatedVector max_coords(coords);

  coords = zero_coords;
  while (coords.data[2] != max_coords.data[2])
    {
      for (int next(list[coords.getMortonNum()]);
	   next >= 0; next = partCellData[next].next)
	if (next != int(part.getID()))
	  func(part, next);
      
      ++coords.data[0];
      if (coords.data[0] > dilatedCellMax[0]) coords.data[0].zero();
      if (coords.data[0] != max_coords.data[0]) continue;
      
      ++coords.data[1];
      coords.data[0] = zero_coords.data[0];
      if (coords.data[1] > dilatedCellMax[1]) coords.data[1].zero();
      if (coords.data[1] != max_coords.data[1]) continue;
      
      ++coords.data[2];
      coords.data[1] = zero_coords.data[1];
      if (coords.data[2] > dilatedCellMax[2]) coords.data[2].zero();
    }
}

void 
CGCellsMorton::getParticleLocalNeighbourhood(const Particle& part, 
				       const nbHoodFunc& func) const
{
  BOOST_FOREACH(const size_t& id, 
		cells[partCellData[part.getID()].cell])
    func(part, id);
}

double 
CGCellsMorton::getMaxSupportedInteractionLength() const
{
  size_t minDiam = 0;

  //As the lambda or overlap is relative to the cellDimension we just
  //find the minimum cell width

  for (size_t i = 1; i < NDIM; ++i)
    if (cellDimension[i] < cellDimension[minDiam])
      minDiam = i;

  return cellLatticeWidth[minDiam] 
    + lambda * (cellLatticeWidth[minDiam] - cellDimension[minDiam]);
}

double 
CGCellsMorton::getMaxInteractionLength() const
{
  return Sim->dynamics.getLongestInteraction();
}

Vector 
CGCellsMorton::calcPosition(const magnet::math::DilatedVector& coords, const Particle& part) const
{
  //We always return the cell that is periodically nearest to the particle
  Vector primaryCell;
  
  for (size_t i(0); i < NDIM; ++i)
    primaryCell[i] = coords.data[i].getRealVal() * cellLatticeWidth[i]
      - 0.5 * Sim->primaryCellSize[i] + cellOffset[i];

  
  Vector imageCell;
  
  for (size_t i = 0; i < NDIM; ++i)
    imageCell[i] = primaryCell[i]
      - Sim->primaryCellSize[i] * lrint((primaryCell[i] - part.getPosition()[i]) 
					/ Sim->primaryCellSize[i]);

  return imageCell;
}

Vector 
CGCellsMorton::calcPosition(const magnet::math::DilatedVector& coords) const
{
  //We always return the cell that is periodically nearest to the particle
  Vector primaryCell;
  
  for (size_t i(0); i < NDIM; ++i)
    primaryCell[i] = coords.data[i].getRealVal() * cellLatticeWidth[i] 
      - 0.5 * Sim->primaryCellSize[i] + cellOffset[i];
  
  return primaryCell;
}
