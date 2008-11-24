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

#include "gcells.hpp"
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

CGCells::CGCells(DYNAMO::SimData* nSim, const std::string& name):
  CGNeighbourList(nSim, "GlobalCellularEvent"),
  cellCount(0),
  cellDimension(1.0),
  lambda(0.9), //Default to higher overlap
  NCells(0)

{
  globName = name;
  I_cout() << "Cells Loaded";
}

CGCells::CGCells(const XMLNode &XML, DYNAMO::SimData* ptrSim):
  CGNeighbourList(ptrSim, "GlobalCellularEvent"),
  cellCount(0),
  cellDimension(1.0),
  lambda(0.9), //Default to higher overlap
  NCells(0)
{
  operator<<(XML);

  I_cout() << "Cells Loaded";
}

CGCells::CGCells(DYNAMO::SimData* ptrSim, const char* nom, void*):
  CGNeighbourList(ptrSim, nom),
  cellCount(0),
  cellDimension(1.0),
  lambda(0.9), //Default to higher overlap
  NCells(0)
{}

void 
CGCells::operator<<(const XMLNode& XML)
{
  try {
    //If you add anything here then it needs to go in gListAndCells.cpp too
    if (XML.isAttributeSet("lambda"))
      lambda = boost::lexical_cast<Iflt>
	(XML.getAttribute("Lambda"));
    
    globName = XML.getAttribute("Name");	
  }
  catch(...)
    {
      D_throw() << "Error loading CGCells";
    }
  
  if (lambda < 0.0 || lambda > 1.0)
    D_throw() << "Lambda out of bounds [0,1), lambda = " << lambda;
}

void 
CGCells::setLambda(const Iflt& nL)
{
  lambda = nL;
}

CGlobEvent 
CGCells::getEvent(const CParticle& part) const
{

  //This 
  //Sim->Dynamics.Liouvillean().updateParticle(part);
  //is not required as we compensate for the delay using 
  //Sim->Dynamics.Liouvillean().getParticleDelay(part)
  
  return CGlobEvent(part,
		    Sim->Dynamics.Liouvillean().
		    getSquareCellCollision2
		    (part, cells[partCellData[part.getID()].cell].origin, 
		     cellDimension)
		    -Sim->Dynamics.Liouvillean().getParticleDelay(part)
		    ,
		    VIRTUAL, *this);
}

void
CGCells::runEvent(const CParticle& part) const
{
  Sim->Dynamics.Liouvillean().updateParticle(part);
  
  size_t oldCell(partCellData[part.getID()].cell);

  //Determine the cell transition direction, its saved
  size_t cellDirection(Sim->Dynamics.Liouvillean().
		       getSquareCellCollision3
		       (part, cells[oldCell].origin, 
			cellDimension));
  size_t endCell;
  size_t inPosition;

  if (std::signbit(part.getVelocity()[cellDirection])) 
    {
      endCell = cells[partCellData[part.getID()].cell].negCells[cellDirection];
      
      inPosition = cells[cells[endCell].negCells[cellDirection]]
	.coords[cellDirection];
    }
  else
    {
      endCell = cells[partCellData[part.getID()].cell].posCells[cellDirection];
      
      inPosition = cells[cells[endCell].posCells[cellDirection]]
	.coords[cellDirection];
    }
  
  //Debug section
#ifdef DYNAMO_WallCollDebug
  {      
    CVector<long> tmp = cells[partCellData[part.getID()].cell].coords;
    CVector<long> tmp2 = cells[endCell].coords;
    
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

  removeFromCell(part.getID());
  addToCell(part.getID(), endCell);

  //Get rid of the virtual event that is next, update is delayed till
  //after all events are added
  Sim->ptrScheduler->popNextEvent();

  //Particle has just arrived into a new cell warn the scheduler about
  //its new neighbours so it can add them to the heap
  BOOST_FOREACH(const int& nb, cells[endCell].neighbours)
    if (static_cast<size_t>(cells[nb].coords[cellDirection]) == inPosition)
      for (int next = cells[nb].list; next != -1; 
	   next = partCellData[next].next)
	sigNewNeighbourNotify(part, next);

  //Tell about the new locals
  BOOST_FOREACH(const size_t& lID, cells[endCell].locals)
    sigNewLocalNotify(part, lID);
  
  //Push the next virtual event, this is the reason the scheduler
  //doesn't need a second callback
  Sim->ptrScheduler->pushEvent(part, getEvent(part));
  Sim->ptrScheduler->sort(part);

  sigCellChangeNotify(part, oldCell);
  
  //This doesn't stream the system as its a virtual event
}

void 
CGCells::initialise(size_t nID)
{
  ID=nID;
  
  reinitialise(Sim->Dynamics.getLongestInteraction());
}

void
CGCells::reinitialise(const Iflt& maxdiam)
{
  I_cout() << "Reinitialising on collision " << Sim->lNColl;

  //Create the cells
  addCells(maxdiam, false);

  addLocalEvents();

  ReInitNotify();
}

void
CGCells::outputXML(xmlw::XmlStream& XML) const
{
  //If you add anything here it also needs to go in gListAndCells.cpp too
  XML << xmlw::attr("Type") << "Cells"
      << xmlw::attr("Lambda") << lambda
      << xmlw::attr("Name") << globName;
}

void
CGCells::addCells(Iflt maxdiam, bool limitCells)
{
  cells.clear();
  partCellData.resize(Sim->lN); //Location data for particles

  NCells = 1;
  cellCount = CVector<long>(0);

  for (int iDim = 0; iDim < NDIM; iDim++)
    {
      cellCount[iDim] = static_cast<long>(Sim->aspectRatio[iDim] / maxdiam);
      
      if (cellCount[iDim] < 3)
	D_throw() << "Not enough cells in " << static_cast<char>('x'+iDim) << " dimension, need 3+";

      if (limitCells && cellCount[iDim] > 100)
	{
	  I_cout() << "Cell count was " << cellCount[iDim] 
		   << "\n Restricting to 100";
	  cellCount[iDim] = 100;
	}

      //Stop bad allocs!
      if (cellCount[iDim] > 500)
	I_cout() << "Cell count was " << cellCount[iDim] 
		 << "\n Restricting to " << (cellCount[iDim] = 500);
            
      NCells *= cellCount[iDim];
    }

  for (int iDim = 0; iDim < NDIM; iDim++)
    cellLatticeWidth[iDim] = Sim->aspectRatio[iDim] / cellCount[iDim];
  
  for (int iDim = 0; iDim < NDIM; iDim++)
    cellDimension[iDim] = cellLatticeWidth[iDim] 
      + (2.0 * cellLatticeWidth[iDim] - maxdiam - cellLatticeWidth[iDim]) 
      * lambda;
  
  I_cout() << "Cells <x,y,z>  " << cellCount[0] << ","
	   << cellCount[1] << "," << cellCount[2];

  I_cout() << "Cells dimension <x,y,z>  " 
	   << cellDimension[0] / Sim->Dynamics.units().unitLength()
	   << ","
	   << cellDimension[1] / Sim->Dynamics.units().unitLength()
	   << "," 
	   << cellDimension[2] / Sim->Dynamics.units().unitLength();

  I_cout() << "Lattice spacing <x,y,z>  " 
	   << cellLatticeWidth[0] / Sim->Dynamics.units().unitLength()
	   << ","
	   << cellLatticeWidth[1] / Sim->Dynamics.units().unitLength()
	   << "," 
	   << cellLatticeWidth[2] / Sim->Dynamics.units().unitLength();

  fflush(stdout);

  cells.resize(NCells); //Empty Cells created!

  for (size_t id = 0; id < NCells; id++)
    {
      cells[id].coords = getCoordsFromID(id);

      for (int iDim = 0; iDim < NDIM; iDim++)
	cells[id].origin[iDim] = cells[id].coords[iDim] 
	  * cellLatticeWidth[iDim] - 0.5 * Sim->aspectRatio[iDim];
      /*
      std::cerr << "\nID " << id << " " << "coords  " << cells[id].coords[0] 
		<< "," << cells[id].coords[1] << "," << cells[id].coords[2]
		<< " " << "Origin  " << cells[id].origin[0] 
		<< "," << cells[id].origin[1] << "," << cells[id].origin[2];*/
    }

  //Add the particles section
  //Required so particles find the right owning cell
  Sim->Dynamics.Liouvillean().updateAllParticles(); 

  ////initialise the data structures
  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    addToCell(part.getID(), getCellID(part.getPosition()));

  //Init the cell linking
  init_cells();
}

void
CGCells::init_cells()
{
  //Holds the displacement in each dimension, the unit is cells!
  CVector<long>displacement(-1);
  //This determines how many neighbours should exist
  long total = static_cast<long>(pow (3.0, static_cast<Iflt>(NDIM)) - 1) / 2;
  
  std::list<CVector<long> > neighbourVectors;
  //This loop iterates through each neighbour position
  for (long iter = 0; iter < total; iter++)
    {
      //Add the current vector to the list
      neighbourVectors.push_back (displacement);
      
      //Now update the displacement vector
      displacement[0] += 1;
      for (int iDim = 1; iDim < NDIM; iDim++)
	if (displacement[iDim - 1] == 2)
	  {
	    displacement[iDim - 1] = -1;
	    displacement[iDim] += 1;
	  }
    }


  for (size_t i = 0; i < NCells; ++i)
    cells[i].neighbours.push_back(i);
  
  //Tell the cells about their neighbours
  for (size_t i = 0; i < NCells; ++i)
    {
      //Tell about the neighbourhood
      BOOST_FOREACH(CVector<long> &neighbour, neighbourVectors)
	{
	  //Positive direction       
	  cells[i].neighbours.push_back(getCellID(cells[i].coords + neighbour));
	  //Negative direction
	  cells[getCellID(cells[i].coords + neighbour)].neighbours.push_back(i);
	}

      CVector<long> tmpvec;
      //Tell them about who's next to them
      for (int iDim = 0; iDim < NDIM; ++iDim)
	{
	  tmpvec = CVector<long>(0);
	  tmpvec[iDim] = 1;
	  cells[i].posCells[iDim] = getCellID(cells[i].coords + tmpvec);
	  cells[i].negCells[iDim] = getCellID(cells[i].coords - tmpvec);
	}
    }
}

void 
CGCells::addLocalEvents()
{
  BOOST_FOREACH(cellStruct& cell, cells)
    {
      cell.locals.clear();

      BOOST_FOREACH(const smrtPlugPtr<CLocal>& local, Sim->Dynamics.getLocals())
	if (local->isInCell(cell.origin, cellDimension))
	  cell.locals.push_back(local->getID());
    }
}

long 
CGCells::getCellID(CVector<long> coords) const
{
  //PBC for vectors
  for (int iDim = 0; iDim < NDIM; iDim++)
    {
      coords[iDim] %= cellCount[iDim];
      if (coords[iDim] < 0)
	coords[iDim] += cellCount[iDim];
    }
  
  return ((coords[0] + coords[1]*cellCount[0] + coords[2]*cellCount[0]*cellCount[1]) % NCells);
}

CVector<long> 
CGCells::getCoordsFromID(unsigned long i) const
{
  CVector<long> tmp;
  i = i % NCells; //PBC's for ID
  
  tmp[0] = i % cellCount[0];
  i /= cellCount[0];
  tmp[1] = i % cellCount[1];
  i /= cellCount[1];
  tmp[2] = i % cellCount[2];
  return tmp;
}

long 
CGCells::getCellID(CVector<> pos) const
{
  Sim->Dynamics.BCs().setPBC(pos);
  CVector<long> temp;
  
  for (int iDim = 0; iDim < NDIM; iDim++)
    temp[iDim] = (long) ( (pos[iDim] + 0.5 * Sim->aspectRatio[iDim]) / cellLatticeWidth[iDim]);
  
  return getCellID(temp);
}


void 
CGCells::getParticleNeighbourhood(const CParticle& part,
				  const nbhoodFunc& func) const
{
  BOOST_FOREACH(const int& nb, cells[partCellData[part.getID()].cell].neighbours)
    for (int next = cells[nb].list;
	 next != -1; next = partCellData[next].next)
      if (next != static_cast<int>(part.getID()))
	func(part, next);
}

void 
CGCells::getParticleLocalNeighbourhood(const CParticle& part, 
				       const nbhoodFunc& func) const
{
  BOOST_FOREACH(const size_t& id, cells[partCellData[part.getID()].cell].locals)
    func(part, id);
}
