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

#include "gcellsmorton.hpp"
#include "globEvent.hpp"
#include "../NparticleEventData.hpp"
#include "../../extcode/xmlParser.h"
#include "../liouvillean/liouvillean.hpp"
#include "../units/units.hpp"
#include "../ranges/1RAll.hpp"
#include "../../schedulers/scheduler.hpp"
#include "../locals/local.hpp"
#include "../BC/LEBC.hpp"
#include "../liouvillean/NewtonianGravityL.hpp"
#include <boost/static_assert.hpp>
#include <magnet/math/ctime_pow.hpp>
#include <magnet/xmlwriter.hpp>
#include <cstdio>


CGCellsMorton::CGCellsMorton(DYNAMO::SimData* nSim, const std::string& name):
  CGNeighbourList(nSim, "MortonCellNeighbourList"),
  cellCount(0),
  cellDimension(1),
  _oversizeCells(1.0),
  NCells(0),
  overlink(1)
{
  globName = name;
  I_cout() << "Cells Loaded";
}

CGCellsMorton::CGCellsMorton(const XMLNode &XML, DYNAMO::SimData* ptrSim):
  CGNeighbourList(ptrSim, "MortonCellNeighbourList"),
  cellCount(0),
  cellDimension(1),
  _oversizeCells(1.0),
  NCells(0),
  overlink(1)
{
  operator<<(XML);

  I_cout() << "Cells Loaded";
}

CGCellsMorton::CGCellsMorton(DYNAMO::SimData* ptrSim, const char* nom, void*):
  CGNeighbourList(ptrSim, nom),
  cellCount(0),
  cellDimension(1),
  _oversizeCells(1.0),
  NCells(0),
  overlink(1)
{}

void 
CGCellsMorton::operator<<(const XMLNode& XML)
{
  try {
    //If you add anything here then it needs to go in gListAndCells.cpp too
    if (XML.isAttributeSet("OverLink"))
      overlink = boost::lexical_cast<size_t>
	(XML.getAttribute("OverLink"));

    if (XML.isAttributeSet("Oversize"))
      _oversizeCells = boost::lexical_cast<double>
	(XML.getAttribute("Oversize"));

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
					 part), 
		     Vector(cellDimension,cellDimension,cellDimension))
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
		       (part, calcPosition(oldCell, part), 
			Vector(cellDimension,cellDimension,cellDimension)));
  
  size_t cellDirection = abs(cellDirectionInt) - 1;

  magnet::math::DilatedVector inCell(oldCell);

  {
    magnet::math::DilatedVector dendCell(inCell);
    
    if (cellDirectionInt > 0)
      {
	++dendCell.data[cellDirection];
	inCell.data[cellDirection] = dendCell.data[cellDirection] + dilatedOverlink;	

	if (dendCell.data[cellDirection] > dilatedCellMax)
	    dendCell.data[cellDirection] = --dendCell.data[cellDirection]
	      - dilatedCellMax;

	if (inCell.data[cellDirection] > dilatedCellMax)
	  inCell.data[cellDirection] = --inCell.data[cellDirection]
	    - dilatedCellMax;
      }
    else
      {
	--dendCell.data[cellDirection];
	inCell.data[cellDirection] = dendCell.data[cellDirection] - dilatedOverlink;

	if (dendCell.data[cellDirection] > dilatedCellMax)
	  dendCell.data[cellDirection] = dendCell.data[cellDirection]
	    - (std::numeric_limits<magnet::math::DilatedInteger>::max() - dilatedCellMax);

	if (inCell.data[cellDirection] > dilatedCellMax)
	  inCell.data[cellDirection] = inCell.data[cellDirection]
	    - (std::numeric_limits<magnet::math::DilatedInteger>::max() - dilatedCellMax);
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
  if (inCell.data[dim1] > dilatedCellMax)
    inCell.data[dim1] = inCell.data[dim1] - (std::numeric_limits<magnet::math::DilatedInteger>::max() - dilatedCellMax);

  if (inCell.data[dim2] > dilatedCellMax) 
    inCell.data[dim2] = inCell.data[dim2] - (std::numeric_limits<magnet::math::DilatedInteger>::max() - dilatedCellMax);

  int walkLength = 2 * overlink + 1;

  const magnet::math::DilatedInteger saved_coord(inCell.data[dim1]);

  //We now have the lowest cell coord, or corner of the cells to update
  for (int iDim(0); iDim < walkLength; ++iDim)
    {
      if (inCell.data[dim2] > dilatedCellMax)
	inCell.data[dim2].zero();

      for (int jDim(0); jDim < walkLength; ++jDim)
	{
	  if (inCell.data[dim1] > dilatedCellMax)
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
  
  if (Sim->dynamics.liouvilleanTypeTest<LNewtonianGravity>())
    I_cout() << "Warning, in order for cellular NB lists to work in gravity\n"
	     << "You must add the ParabolaSentinel Global event.";

  reinitialise(getMaxInteractionLength());
}

void
CGCellsMorton::reinitialise(const double& maxdiam)
{
  I_cout() << "Reinitialising on collision " << Sim->eventCount;

  //Create the cells
  addCells(_oversizeCells * (maxdiam * (1.0 + 10 * std::numeric_limits<double>::epsilon())) / overlink);

  addLocalEvents();

  BOOST_FOREACH(const initSlot& nbs, sigReInitNotify)
    nbs.second();

  if (isUsedInScheduler)
    Sim->ptrScheduler->initialise();
}

void
CGCellsMorton::outputXML(xml::XmlStream& XML) const
{
  //If you add anything here it also needs to go in gListAndCells.cpp too
  XML << xml::attr("Type") << "CellsMorton"
      << xml::attr("Name") << globName;

  if (overlink > 1)   XML << xml::attr("OverLink") << overlink;
  if (_oversizeCells != 1.0) XML << xml::attr("Oversize") << _oversizeCells;

}

void
CGCellsMorton::addCells(double maxdiam)
{
  cells.clear();
  partCellData.resize(Sim->N); //Location data for particles

  NCells = 1;
  cellCount = 0;

#ifdef DYNAMO_DEBUG
  for (size_t iDim = 0; iDim < NDIM; iDim++)
    if (Sim->aspectRatio[iDim] != 1.0) 
      I_cerr() << "Warning! A non-square system is certainly not optimal for memory usage when using Morton Ordered Lists";
#endif

  cellCount = int(1 / maxdiam);
  
  if (cellCount < 3)
    M_throw() << "Not enough cells, sim too small, need 3+";

  if (cellCount > std::numeric_limits<unsigned char>::max())
    {
      I_cout() << "Cell count was " << cellCount
	       << "\n Restricting to " << std::numeric_limits<unsigned char>::max() 
	       << " to stop huge amounts of memory being allocated";

      cellCount = std::numeric_limits<unsigned char>::max();
    }

  dilatedCellMax = cellCount - 1;

  dilatedOverlink = overlink;

  NCells = cellCount * cellCount * cellCount;

  cellLatticeWidth = 1.0 / cellCount;
  
  cellDimension = cellLatticeWidth + (cellLatticeWidth - maxdiam) 
      * lambda;
  
  cellOffset = -(cellLatticeWidth - maxdiam) * lambda * 0.5;


  I_cout() << "Cells <N>  " << NCells;

  I_cout() << "Cells dimension <x>  " 
	   << cellDimension / Sim->dynamics.units().unitLength();

  I_cout() << "Lattice spacing <x,y,z>  " 
	   << cellLatticeWidth / Sim->dynamics.units().unitLength();

  fflush(stdout);
  
  //Find the required size of the morton array
  size_t sizeReq(1);

  for (int i(0); i < magnet::math::ctime_pow<2,magnet::math::DilatedInteger::digits>::result; ++i)
    {
      sizeReq *= 2*2*2;
      if (sizeReq >= NCells) break;
    }

  cells.resize(sizeReq); //Empty Cells created!
  list.resize(sizeReq); //Empty Cells created!

  I_cout() << "Vector Size <N>  " << sizeReq;
  
  for (size_t iDim = 0; iDim < cellCount; ++iDim)
    for (size_t jDim = 0; jDim < cellCount; ++jDim)
      for (size_t kDim = 0; kDim < cellCount; ++kDim)
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
  Vector cellDimensionsVector(cellDimension,cellDimension,cellDimension);
  for (size_t iDim = 0; iDim < cellCount; ++iDim)
    for (size_t jDim = 0; jDim < cellCount; ++jDim)
      for (size_t kDim = 0; kDim < cellCount; ++kDim)
	{
	  magnet::math::DilatedVector coords(iDim, jDim, kDim);
	  size_t id = coords.getMortonNum();
	  cells[id].clear();
	  Vector pos = calcPosition(coords);
	  
	  //We make the box slightly larger to ensure objects on the boundary are included
	  BOOST_FOREACH(const magnet::ClonePtr<Local>& local, Sim->dynamics.getLocals())
	    if (local->isInCell(pos - 0.0001 * cellDimensionsVector, 1.0002 * cellDimensionsVector))
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
      coords[iDim] %= cellCount;
      if (coords[iDim] < 0) coords[iDim] += cellCount;
    }
  
  return magnet::math::DilatedVector(coords[0],coords[1],coords[2]);
}

magnet::math::DilatedVector
CGCellsMorton::getCellID(Vector pos) const
{
  Sim->dynamics.BCs().applyBC(pos);
  CVector<int> temp;
  
  for (size_t iDim = 0; iDim < NDIM; iDim++)
    temp[iDim] = std::floor((pos[iDim] + 0.5 * Sim->aspectRatio[iDim] - cellOffset)
			    / cellLatticeWidth);
  
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
      if (coords.data[iDim] > dilatedCellMax)
	coords.data[iDim] -= 
	  std::numeric_limits<magnet::math::DilatedInteger>::max() - dilatedCellMax;
    }

  const magnet::math::DilatedVector zero_coords(coords);

  coords = center_coords;
  for (size_t iDim(0); iDim < NDIM; ++iDim)
    {
      coords.data[iDim] += dilatedOverlink + 1;
      if (coords.data[iDim] > dilatedCellMax)
	coords.data[iDim] -= dilatedCellMax + 1;
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
      if (coords.data[0] > dilatedCellMax) coords.data[0].zero();
      if (coords.data[0] != max_coords.data[0]) continue;
      
      ++coords.data[1];
      coords.data[0] = zero_coords.data[0];
      if (coords.data[1] > dilatedCellMax) coords.data[1].zero();
      if (coords.data[1] != max_coords.data[1]) continue;
      
      ++coords.data[2];
      coords.data[1] = zero_coords.data[1];
      if (coords.data[2] > dilatedCellMax) coords.data[2].zero();
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
  return cellLatticeWidth + lambda * (cellLatticeWidth - cellDimension);
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
    primaryCell[i] = coords.data[i].getRealVal() * cellLatticeWidth - 0.5 * Sim->aspectRatio[i] + cellOffset;

  
  Vector imageCell;
  
  for (size_t i = 0; i < NDIM; ++i)
    imageCell[i] = primaryCell[i]
      - Sim->aspectRatio[i] * rintfunc((primaryCell[i] - part.getPosition()[i]) / Sim->aspectRatio[i]);

  return imageCell;
}

Vector 
CGCellsMorton::calcPosition(const magnet::math::DilatedVector& coords) const
{
  //We always return the cell that is periodically nearest to the particle
  Vector primaryCell;
  
  for (size_t i(0); i < NDIM; ++i)
    primaryCell[i] = coords.data[i].getRealVal() * cellLatticeWidth 
      - 0.5 * Sim->aspectRatio[i] + cellOffset;
  
  return primaryCell;
}
