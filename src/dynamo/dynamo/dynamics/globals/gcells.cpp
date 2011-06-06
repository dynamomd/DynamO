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

#include "gcells.hpp"
#include "globEvent.hpp"
#include "../NparticleEventData.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../units/units.hpp"
#include "../ranges/1RAll.hpp"
#include "../../schedulers/scheduler.hpp"
#include "../locals/local.hpp"
#include "../BC/LEBC.hpp"
#include "../liouvillean/NewtonianGravityL.hpp"
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/static_assert.hpp>
#include <cstdio>

Vector 
CGCells::calcPosition(const Vector& primaryCell,
		      const Particle& part) const
{
  //We always return the cell that is periodically nearest to the particle
  Vector imageCell;
  
  for (size_t i = 0; i < NDIM; ++i)
    imageCell[i] = primaryCell[i] - Sim->primaryCellSize[i]
      * rintfunc((primaryCell[i] - part.getPosition()[i]) / Sim->primaryCellSize[i]);
  
  return imageCell;
}


CGCells::CGCells(dynamo::SimData* nSim, const std::string& name, 
		 const size_t& overlink):
  CGNeighbourList(nSim, "CellNeighbourList"),
  cellCount(0),
  cellDimension(1,1,1),
  _oversizeCells(1.0),
  NCells(0),
  overlink(overlink),
  interaction(""),
  MaxIntDist(0.0)
{
  globName = name;
  I_cout() << "Cells Loaded, Overlinking set to " << overlink;
}

CGCells::CGCells(const magnet::xml::Node &XML, dynamo::SimData* ptrSim):
  CGNeighbourList(ptrSim, "CellNeighbourList"),
  cellCount(0),
  cellDimension(1,1,1),
  _oversizeCells(1.0),
  NCells(0),
  overlink(1),
  interaction(""),
  MaxIntDist(0.0)
{
  operator<<(XML);

  I_cout() << "Cells Loaded";
}

CGCells::CGCells(dynamo::SimData* ptrSim, const char* nom, void*):
  CGNeighbourList(ptrSim, nom),
  cellCount(0),
  cellDimension(1,1,1),
  _oversizeCells(1.0),
  NCells(0),
  overlink(1),
  interaction(""),
  MaxIntDist(0.0)
{}

void 
CGCells::operator<<(const magnet::xml::Node& XML)
{
  try {
    if (XML.hasAttribute("OverLink"))
      overlink = XML.getAttribute("OverLink").as<size_t>();

    if (XML.hasAttribute("Oversize"))
      _oversizeCells = XML.getAttribute("Oversize").as<double>();

    if (_oversizeCells < 1.0)
      M_throw() << "You must specify an Oversize greater than 1.0, otherwise your cells are too small!";

    if (XML.hasAttribute("Interaction"))
      interaction = XML.getAttribute("Interaction");

    if (XML.hasAttribute("CellWidth"))
      MaxIntDist = XML.getAttribute("CellWidth").as<double>() * Sim->dynamics.units().unitLength();
    
    globName = XML.getAttribute("Name");	
  }
  catch(...)
    {
      M_throw() << "Error loading CGCells";
    }
}


GlobalEvent 
CGCells::getEvent(const Particle& part) const
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
		     (part, calcPosition(cells[partCellData[part.getID()].cell].origin, part), 
		      cellDimension)
		     -Sim->dynamics.getLiouvillean().getParticleDelay(part)
		     ,
		     CELL, *this);
}

void
CGCells::runEvent(const Particle& part, const double) const
{
  //Despite the system not being streamed this must be done.  This is
  //because the scheduler and all interactions, locals and systems
  //expect the particle to be up to date.
  Sim->dynamics.getLiouvillean().updateParticle(part);
  
  size_t oldCell(partCellData[part.getID()].cell);

  //Determine the cell transition direction
  int cellDirectionInt(Sim->dynamics.getLiouvillean().
		       getSquareCellCollision3
		       (part, calcPosition(cells[oldCell].origin, part), 
			cellDimension));
  
  size_t cellDirection = abs(cellDirectionInt) - 1;

  int endCell(oldCell), inCell(oldCell);

  {
    size_t cellpow(1);
    for (size_t iDim(0); iDim < cellDirection; ++iDim)
      cellpow *= cellCount[iDim];

    int mag = cellpow * cellCount[cellDirection];
  
    if (cellDirectionInt > 0)
      {
	endCell += cellpow;
	inCell += (1 + overlink) * cellpow;
	if (cells[oldCell].coords[cellDirection] == cellCount[cellDirection] - 1)
	  {
	    endCell -= mag;
	    inCell  -= mag;
	  }
	else if (cells[oldCell].coords[cellDirection] >= 
		 int(cellCount[cellDirection] - 1 - overlink))
	  inCell  -= mag;
      }
    else
      {
	endCell -= cellpow;
	inCell -= (1+overlink) * cellpow;
	
	if (cells[oldCell].coords[cellDirection] == 0)
	  {
	    endCell += mag;
	    inCell  += mag;
	  }
	else if (cells[oldCell].coords[cellDirection] <= int(overlink))
	  inCell  += mag;
      }
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
  
  CVector<int> coords(cells[inCell].coords);
  coords[dim1] -= overlink;
  coords[dim2] -= overlink;

  if (coords[dim1] < 0) coords[dim1] += cellCount[dim1];
  if (coords[dim2] < 0) coords[dim2] += cellCount[dim2];

  int nb(getCellIDprebounded(coords));

  size_t dim1pow(1), dim2pow(1);
  
  for (size_t iDim(0); iDim < dim1; ++iDim)
    dim1pow *= cellCount[iDim];

  for (size_t iDim(0); iDim < dim2; ++iDim)
    dim2pow *= cellCount[iDim];

  int walkLength = 2 * overlink+1;

  //We now have the lowest cell coord, or corner of the cells to update
  for (int iDim(0); iDim < walkLength; ++iDim)
    {
      if (coords[dim2] + iDim == cellCount[dim2]) 
	nb -= dim2pow * cellCount[dim2];

      for (int jDim(0); jDim < walkLength; ++jDim)
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

      if (coords[dim1] + walkLength - 1 >= cellCount[dim1]) nb += dim1pow * cellCount[dim1];
      
      nb += dim2pow - walkLength * dim1pow;
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
    CVector<int> tmp2 = cells[partCellData[part.getID()].cell].coords;
    CVector<int> tmp = cells[oldCell].coords;
    
    std::cerr << "\nCGCells sysdt " 
	      << Sim->dSysTime / Sim->dynamics.units().unitTime()
	      << "  Global ID "
	      << part.getID()
	      << "  from <" 
	      << tmp[0] << "," << tmp[1] << "," << tmp[2]
	      << "> to <" 
	      << tmp2[0] << "," << tmp2[1] << "," << tmp2[2] << ">";
  }
#endif
}

void 
CGCells::initialise(size_t nID)
{
  lambda=0.99;
  ID=nID;
  
  reinitialise(getMaxInteractionLength());

  if (Sim->dynamics.liouvilleanTypeTest<LNewtonianGravity>())
    I_cout() << "Warning, in order for cellular NB lists to work in gravity\n"
	     << "You must add the ParabolaSentinel Global event.";
}

void
CGCells::reinitialise(const double& maxdiam)
{
  I_cout() << "Reinitialising on collision " << Sim->eventCount;

  //Create the cells
  addCells(_oversizeCells * maxdiam / overlink);

  addLocalEvents();

  BOOST_FOREACH(const initSlot& nbs, sigReInitNotify)
    nbs.second();
  
  if (isUsedInScheduler)
    Sim->ptrScheduler->initialise();
}

void 
CGCells::outputXML(xml::XmlStream& XML) const
{
  outputXML(XML, "Cells");
}

void
CGCells::outputXML(xml::XmlStream& XML, const std::string& name) const
{
  //If you add anything here it also needs to go in gListAndCells.cpp too
  XML << xml::attr("Type") << name
      << xml::attr("Name") << globName;

  if (MaxIntDist != 0.0)
    XML << xml::attr("CellWidth") << MaxIntDist / Sim->dynamics.units().unitLength();
  else if (!interaction.empty())
    XML << xml::attr("Interaction") << interaction;
      
  if (overlink > 1)   XML << xml::attr("OverLink") << overlink;
  if (_oversizeCells != 1.0) XML << xml::attr("Oversize") << _oversizeCells;
}

void
CGCells::addCells(double maxdiam)
{
  cells.clear();
  partCellData.resize(Sim->N); //Location data for particles

  NCells = 1;
  cellCount = CVector<int>(0);

  for (size_t iDim = 0; iDim < NDIM; iDim++)
    {
      cellCount[iDim] = int(Sim->primaryCellSize[iDim] / (maxdiam * (1.0 + 10 * std::numeric_limits<double>::epsilon())));
      
      if (cellCount[iDim] < 3)
	M_throw() << "Not enough cells in " << char('x'+iDim) << " dimension, need 3+";

//      if (cellCount[iDim] > 400)
//	{
//	  I_cout() << "Cell count was " << cellCount[iDim]
//		   << "\n Restricting to 400 to stop this sim failing";
//	  cellCount[iDim] = 400;
//	}

      NCells *= cellCount[iDim];
    }

  for (size_t iDim = 0; iDim < NDIM; iDim++)
    cellLatticeWidth[iDim] = Sim->primaryCellSize[iDim] / cellCount[iDim];
  
  for (size_t iDim = 0; iDim < NDIM; iDim++)
    cellDimension[iDim] = cellLatticeWidth[iDim] 
      + (cellLatticeWidth[iDim] - maxdiam) 
      * lambda;
  
  //This is just to center the grid of cells about the origin (0,0,0) of the system
  for (size_t iDim = 0; iDim < NDIM; iDim++)
    cellOffset[iDim] = -(cellLatticeWidth[iDim] - maxdiam) * lambda / 2;

  I_cout() << "Cells <x,y,z>  " << cellCount[0] << ","
	   << cellCount[1] << "," << cellCount[2];

  I_cout() << "Cell Offset <x,y,z>  "
           << cellOffset[0] / Sim->dynamics.units().unitLength() << ","
	   << cellOffset[1] / Sim->dynamics.units().unitLength() << ","
	   << cellOffset[2] / Sim->dynamics.units().unitLength();

  I_cout() << "Cells Dimension <x,y,z>  " 
	   << cellDimension[0] / Sim->dynamics.units().unitLength()
	   << ","
	   << cellDimension[1] / Sim->dynamics.units().unitLength()
	   << "," 
	   << cellDimension[2] / Sim->dynamics.units().unitLength();

  I_cout() << "Lattice spacing <x,y,z>  " 
	   << cellLatticeWidth[0] / Sim->dynamics.units().unitLength()
	   << ","
	   << cellLatticeWidth[1] / Sim->dynamics.units().unitLength()
	   << "," 
	   << cellLatticeWidth[2] / Sim->dynamics.units().unitLength();

  try {
    cells.resize(NCells); //Empty Cells created!
  }
  catch(std::bad_alloc& er)
    {
      M_throw() << "The number of cells is causing a bad alloc\n"
		<< "Number of cells (" << NCells << ") and thus the system size could be too large\n"
		<< "Max is " << cells.max_size() << ", aborting"
		<< er.what();
    }

  for (size_t id = 0; id < NCells; id++)
    {
      cells[id].coords = getCoordsFromID(id);

      for (size_t iDim = 0; iDim < NDIM; iDim++)
	cells[id].origin[iDim] = cells[id].coords[iDim] 
	  * cellLatticeWidth[iDim] 
	  - 0.5 * Sim->primaryCellSize[iDim]
	  + cellOffset[iDim];
    }

  //Add the particles section
  //Required so particles find the right owning cell
  Sim->dynamics.getLiouvillean().updateAllParticles(); 

#ifdef DYNAMO_WallCollDebug
  BOOST_FOREACH(const Particle& part, Sim->particleList)
    {
      CVector<int> cellcoords = getCoordsFromID(getCellID(part.getPosition()));
      std::cerr << "\n Adding ID=" << part.getID() 
		<< " pos="  << part.getPosition()[0]
		<< ","  << part.getPosition()[1]
		<< ","  << part.getPosition()[2]
		<< " cellID=" << getCellID(part.getPosition())
		<< " cellCoords=" << cellcoords[0]
		<< "," << cellcoords[1]
		<< "," << cellcoords[2]
	;
    }
#endif

  ////initialise the data structures
  BOOST_FOREACH(const Particle& part, Sim->particleList)    
    addToCell(part.getID(), getCellID(part.getPosition()));
}

void 
CGCells::addLocalEvents()
{
  BOOST_FOREACH(cellStruct& cell, cells)
    {
      cell.locals.clear();

      //We make the box slightly larger to ensure objects on the boundary are included
      BOOST_FOREACH(const magnet::ClonePtr<Local>& local, Sim->dynamics.getLocals())
	if (local->isInCell(cell.origin - 0.0001 * cellDimension, 1.0002 * cellDimension))
	  cell.locals.push_back(local->getID());
    }
}

size_t
CGCells::getCellID(const CVector<int>& coordsold) const
{
  //PBC for vectors
  CVector<int> coords(coordsold);

  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    {
      coords[iDim] %= cellCount[iDim];
      if (coords[iDim] < 0) coords[iDim] += cellCount[iDim];
    }
  
  return getCellIDprebounded(coords);
}

size_t
CGCells::getCellIDprebounded(const CVector<int>& coords) const
{
  int nb(coords[0]);

  size_t pow(cellCount[0]);

  for (size_t iDim(1); iDim < NDIM-1; ++iDim)
    {
      nb += coords[iDim] * pow;
      pow *= cellCount[iDim];
    }
    
  return nb + coords[NDIM-1] * pow;
}

CVector<int> 
CGCells::getCoordsFromID(size_t i) const
{
  CVector<int> tmp;
  i = i % NCells; //PBC's for ID //NOT NEEDED
  
  tmp[0] = i % cellCount[0];
  i /= cellCount[0];
  tmp[1] = i % cellCount[1];
  i /= cellCount[1];
  tmp[2] = i % cellCount[2];
  return tmp;
}

size_t
CGCells::getCellID(Vector  pos) const
{
  Sim->dynamics.BCs().applyBC(pos);
  CVector<int> temp;
  
  for (size_t iDim = 0; iDim < NDIM; iDim++)
    temp[iDim] = std::floor((pos[iDim] + 0.5 * Sim->primaryCellSize[iDim] - cellOffset[iDim])
		     / cellLatticeWidth[iDim]);
  
  return getCellID(temp);
}


void 
CGCells::getParticleNeighbourhood(const Particle& part,
				  const nbHoodFunc& func) const
{
  CVector<int> coords(cells[partCellData[part.getID()].cell].coords);


  for (size_t iDim(0); iDim < NDIM; ++iDim)
    {
      coords[iDim] -= overlink;
      if (coords[iDim] < 0)
	coords[iDim] += cellCount[iDim];
    }
  
  int nb(getCellIDprebounded(coords));

  //This loop iterates through each neighbour position
  BOOST_STATIC_ASSERT(NDIM==3);

  int walkLength(2*overlink+1);

  for (int iDim(0); iDim < walkLength; ++iDim)
    {
      if (coords[2] + iDim == cellCount[2])
	nb -= NCells;

      for (int jDim(0); jDim < walkLength; ++jDim)
	{
	  if (coords[1] + jDim == cellCount[1])
	    nb -=  cellCount[1] * cellCount[0];
	  
	  for (int kDim(0); kDim < walkLength; ++kDim)
	    {
 	      if (coords[0] + kDim == cellCount[0])
		nb -= cellCount[0];
	      
	      for (int next(cells[nb++].list);
		   next >= 0; next = partCellData[next].next)
		if (next != int(part.getID()))
		  func(part, next);
	    }

	  nb += (1 + (coords[0] + walkLength - 1 >= cellCount[0])) 
	    * cellCount[0] - walkLength;
	}

      nb += ((1 + (coords[1] + walkLength - 1 >= cellCount[1])) 
	     * cellCount[1] - walkLength) * cellCount[0];
    }
}

void 
CGCells::getParticleLocalNeighbourhood(const Particle& part, 
				       const nbHoodFunc& func) const
{
  BOOST_FOREACH(const size_t& id, 
		cells[partCellData[part.getID()].cell].locals)
    func(part, id);
}

double 
CGCells::getMaxSupportedInteractionLength() const
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
CGCells::getMaxInteractionLength() const
{
  if (MaxIntDist != 0.0)
    return MaxIntDist;
  else if (!interaction.empty())
    return Sim->dynamics.getInteraction(interaction)->maxIntDist();
  else 
    return Sim->dynamics.getLongestInteraction();
}
