/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
#include "../../extcode/xmlwriter.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../units/units.hpp"
#include "../ranges/1RAll.hpp"
#include "../../schedulers/scheduler.hpp"
#include "../locals/local.hpp"
#include "../BC/LEBC.hpp"
#include <boost/static_assert.hpp>
#include <cstdio>
#include "../liouvillean/NewtonianGravityL.hpp"

CGCellsMorton::CGCellsMorton(DYNAMO::SimData* nSim, const std::string& name):
  CGNeighbourList(nSim, "MortonCellNeighbourList"),
  cellCount(0),
  cellDimension(1),
  lambda(0.9), //Default to higher overlap
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
  lambda(0.9), //Default to higher overlap
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
  lambda(0.9), //Default to higher overlap
  NCells(0),
  overlink(1)
{}

void 
CGCellsMorton::operator<<(const XMLNode& XML)
{
  try {
    //If you add anything here then it needs to go in gListAndCells.cpp too
    if (XML.isAttributeSet("Lambda"))
      lambda = boost::lexical_cast<Iflt>
	(XML.getAttribute("Lambda"));

    if (XML.isAttributeSet("OverLink"))
      overlink = boost::lexical_cast<size_t>
	(XML.getAttribute("OverLink"));
    
    globName = XML.getAttribute("Name");	
  }
  catch(...)
    {
      M_throw() << "Error loading CGCellsMorton";
    }
  
  if (lambda < 0.0 || lambda > 1.0)
    M_throw() << "Lambda out of bounds [0,1), lambda = " << lambda;
}

void 
CGCellsMorton::setLambda(const Iflt& nL)
{
  lambda = nL;
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
		    (part, calcPosition(partCellData[part.getID()].cell), 
		     Vector(cellDimension,cellDimension,cellDimension))
		    -Sim->dynamics.getLiouvillean().getParticleDelay(part)
		    ,
		    CELL, *this);
}

void
CGCellsMorton::runEvent(const Particle& part) const
{

  //Despite the system not being streamed this must be done.  This is
  //because the scheduler and all interactions, locals and systems
  //expect the particle to be up to date.
  Sim->dynamics.getLiouvillean().updateParticle(part);

  size_t endCell, oldCell(partCellData[part.getID()].cell);
  
  //Determine the cell transition direction, its saved
  int cellDirectionInt(Sim->dynamics.getLiouvillean().
			  getSquareCellCollision3
			  (part, calcPosition(oldCell), 
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
	    - (magnet::math::DilatedInteger(magnet::math::DilatedInteger::DilatedMaxVal,0) - dilatedCellMax);

	if (inCell.data[cellDirection] > dilatedCellMax)
	  inCell.data[cellDirection] = inCell.data[cellDirection]
	    - (magnet::math::DilatedInteger(magnet::math::DilatedInteger::DilatedMaxVal,0) - dilatedCellMax);
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
    inCell.data[dim1] = inCell.data[dim1] - (magnet::math::DilatedInteger(magnet::math::DilatedInteger::DilatedMaxVal,0) - dilatedCellMax);

  if (inCell.data[dim2] > dilatedCellMax) 
    inCell.data[dim2] = inCell.data[dim2] - (magnet::math::DilatedInteger(magnet::math::DilatedInteger::DilatedMaxVal,0) - dilatedCellMax);

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
    CVector<int> tmp = cells[partCellData[part.getID()].cell].coords;
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
CGCellsMorton::initialise(size_t nID)
{
  ID=nID;
  
  if (Sim->dynamics.liouvilleanTypeTest<LNewtonianGravity>())
    I_cout() << "Warning, in order for cellular NB lists to work in gravity\n"
	     << "You must add the ParabolaSentinel Global event.";

  reinitialise(getMaxInteractionLength());
}

void
CGCellsMorton::reinitialise(const Iflt& maxdiam)
{
  I_cout() << "Reinitialising on collision " << Sim->eventCount;

  //Create the cells
  addCells(maxdiam / overlink);

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
      << xml::attr("Lambda") << lambda
      << xml::attr("Name") << globName;

  if (overlink > 1)   XML << xml::attr("OverLink") << overlink;
}

void
CGCellsMorton::addCells(Iflt maxdiam)
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
  
  I_cout() << "Cells <N>  " << NCells;

  I_cout() << "Cells dimension <x>  " 
	   << cellDimension / Sim->dynamics.units().unitLength();

  I_cout() << "Lattice spacing <x,y,z>  " 
	   << cellLatticeWidth / Sim->dynamics.units().unitLength();

  fflush(stdout);
  
  //Find the required size of the morton array
  size_t sizeReq(1);

  for (int i(0); i < ctime_pow<2,magnet::math::DilatedInteger::S>::result; ++i)
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
  for (size_t iDim = 0; iDim < cellCount; ++iDim)
    for (size_t jDim = 0; jDim < cellCount; ++jDim)
      for (size_t kDim = 0; kDim < cellCount; ++kDim)
	{
	  magnet::math::DilatedVector coords(iDim, jDim, kDim);
	  size_t id = coords.getMortonNum();
	  cells[id].clear();
	  Vector pos = calcPosition(coords);
	  
	  BOOST_FOREACH(const ClonePtr<Local>& local, Sim->dynamics.getLocals())
	    if (local->isInCell(pos, Vector(cellDimension,cellDimension,cellDimension)))
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
CGCellsMorton::getCellID(Vector  pos) const
{
  Sim->dynamics.BCs().applyBC(pos);
  CVector<int> temp;
  
  for (size_t iDim = 0; iDim < NDIM; iDim++)
    temp[iDim] = int((pos[iDim] + 0.5 * Sim->aspectRatio[iDim]) 
		     / cellLatticeWidth);
  
  return getCellID(temp);
}


void 
CGCellsMorton::getParticleNeighbourhood(const Particle& part,
					const nbHoodFunc& func) const
{
  magnet::math::DilatedVector coords(partCellData[part.getID()].cell);

  for (size_t iDim(0); iDim < NDIM; ++iDim)
    {
      coords.data[iDim] = coords.data[iDim] - dilatedOverlink;

      if (coords.data[iDim] > dilatedCellMax) 
	coords.data[iDim] = coords.data[iDim]
	  - (magnet::math::DilatedInteger(magnet::math::DilatedInteger::DilatedMaxVal,0) - dilatedCellMax);
    }
  
  //This loop iterates through each neighbour position
  BOOST_STATIC_ASSERT(NDIM==3);

  int walkLength(2*overlink+1);

  const magnet::math::DilatedVector stored_coords(coords);

  for (int iDim(0); iDim < walkLength; ++iDim)
    {
      if (coords.data[2] > dilatedCellMax) 
	coords.data[2].zero();

      for (int jDim(0); jDim < walkLength; ++jDim)
	{	
	  if (coords.data[1] > dilatedCellMax) 
	    coords.data[1].zero();
  
	  for (int kDim(0); kDim < walkLength; ++kDim)
	    {	      
	      if (coords.data[0] > dilatedCellMax)
		coords.data[0].zero();

	      for (int next(list[coords.getMortonNum()]);
		   next >= 0; next = partCellData[next].next)
		if (next != int(part.getID()))
		  func(part, next);

	      ++coords.data[0];
	    }

	  coords.data[0] = stored_coords.data[0];

	  ++coords.data[1];
	}
      
      coords.data[1] = stored_coords.data[1];

      ++coords.data[2];
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

Iflt 
CGCellsMorton::getMaxSupportedInteractionLength() const
{
  return cellLatticeWidth + lambda * (cellLatticeWidth - cellDimension);
}

Iflt 
CGCellsMorton::getMaxInteractionLength() const
{
  return Sim->dynamics.getLongestInteraction();
}

Vector 
CGCellsMorton::calcPosition(const magnet::math::DilatedVector& coords) const
{
  return Vector(coords.data[0].getRealVal() * cellLatticeWidth - 0.5 * Sim->aspectRatio[0],
		coords.data[1].getRealVal() * cellLatticeWidth - 0.5 * Sim->aspectRatio[1],
		coords.data[2].getRealVal() * cellLatticeWidth - 0.5 * Sim->aspectRatio[2]);
}
