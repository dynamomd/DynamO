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

CGCells::CGCells(const DYNAMO::SimData* nSim):
  CGlobal(new CRAll(nSim), nSim),
  cellCount(0),
  cellDimension(1.0),
  lambda(0.9), //Default to higher overlap
  NCells(0)

{
  globName = "Cells";
  I_cout() << "Cells Loaded";
}

CGCells::CGCells(const XMLNode &XML, const DYNAMO::SimData* ptrSim):
  CGlobal(ptrSim),
  cellCount(0),
  cellDimension(1.0),
  lambda(0.9), //Default to higher overlap
  NCells(0)
{
  globName = "Cells";
  operator<<(XML);

  I_cout() << "Cells Loaded";
}

void 
CGCells::operator<<(const XMLNode& XML)
{
  if (XML.isAttributeSet("lambda"))
    {
      try {
	
	lambda = boost::lexical_cast<Iflt>
	  (XML.getAttribute("lambda"));
	
      }
      catch(...)
	{
	  D_throw() << "Could not load the lambda value in cellular scheduler";
	}
      
      if (lambda < 0.0 || lambda > 1.0)
	D_throw() << "Lambda out of bounds [0,1), lambda = " << lambda;
    }
}

size_t 
CGCells::getLocalCellID(const CParticle& part) const
{
  return partCellData[part.getID()].cell;
}

CGlobEvent 
CGCells::getEvent(const CParticle& part) const
{
  //Update the wall collision

  return CGlobEvent(part, Sim->Dynamics.Liouvillean().
		    getSquareCellCollision2
		    (part, cells[partCellData[part.getID()].cell].origin, 
		     cellDimension), VIRTUAL, *this);
}

CNParticleData
CGCells::runEvent(const CGlobEvent& event) const
{
  const CParticle& part(event.getParticle());

  //Determine the cell transition direction, its saved
  size_t cellDirection(Sim->Dynamics.Liouvillean().
		    getSquareCellCollision3
		    (part, cells[partCellData[part.getID()].cell].origin, 
		     cellDimension));
  size_t endCell;
  size_t inPosition;

  if (std::signbit(part.getVelocity()[cellDirection])) 
    {
      endCell = cells[partCellData[part.getID()].cell].negCells[cellDirection];
      inPosition = cells[cells[endCell].negCells[cellDirection]].coords[cellDirection];
    }
  else
    {
      endCell = cells[partCellData[part.getID()].cell].posCells[cellDirection];
      inPosition = cells[cells[endCell].posCells[cellDirection]].coords[cellDirection];
    }
  
  //Debug section
#ifdef DYNAMO_WallCollDebug
  {      
    CVector<long> tmp = cells[partCellData[part.getID()].cell].coords;
    CVector<long> tmp2 = cells[endCell].coords;
    
    std::cerr << "\nCGWall sysdt " 
	      << (Sim->dSysTime + event.getdt())
	      << "  WALL ID "
	      << part.getID()
	      << "  dt " << event.getdt()
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
  Sim->ptrScheduler->popVirtualEvent();

  //Particle has just arrived into a new cell warn the scheduler about
  //its new neighbours so it can add them to the heap
  BOOST_FOREACH(const int& nb, cells[endCell].neighbours)
    if (cells[nb].coords[cellDirection] == inPosition)
      for (int next = cells[nb].list; next != -1; next = partCellData[next].next)
	Sim->ptrScheduler->virtualCellNewNeighbour(part, Sim->vParticleList[next]);

  //Push the next virtual event and update the heap
  Sim->ptrScheduler->pushAndUpdateVirtualEvent(part, CGCells::getEvent(part));

  return CNParticleData();
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
}

void
CGCells::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "Cells";
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
  I_cout() << "Cells dimension <x,y,z>  " << cellDimension[0] << ","
	   << cellDimension[1] << "," << cellDimension[2];
  I_cout() << "Lattice spacing <x,y,z>  " << cellLatticeWidth[0] << ","
	   << cellLatticeWidth[1] << "," << cellLatticeWidth[2];

  fflush(stdout);

  cells.resize(NCells); //Empty Cells created!

  for (size_t id = 0; id < NCells; id++)
    {
      cells[id].coords = getCoordsFromID(id);

      for (int iDim = 0; iDim < NDIM; iDim++)
	cells[id].origin[iDim] = cells[id].coords[iDim] * cellLatticeWidth[iDim] - 0.5 * Sim->aspectRatio[iDim];

      /*std::cerr << "\nID " << id << " " << "coords  " << cells[id].coords[0] 
		<< "," << cells[id].coords[1] << "," << cells[id].coords[2]
		<< " " << "Origin  " << cells[id].origin[0] 
		<< "," << cells[id].origin[1] << "," << cells[id].origin[2];*/
    }

  //Add the particles section
  //Required so particles find the right owning cell
  Sim->Dynamics.Liouvillean().updateAllParticles(); 

  ////initialise the data structures
  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    addToCell(part.getID(), getID(part.getPosition()));

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
	  cells[i].neighbours.push_back(getID(cells[i].coords + neighbour));
	  //Negative direction
	  cells[getID(cells[i].coords + neighbour)].neighbours.push_back(i);
	}

      CVector<long> tmpvec;
      //Tell them about who's next to them
      for (int iDim = 0; iDim < NDIM; ++iDim)
	{
	  tmpvec = CVector<long>(0);
	  tmpvec[iDim] = 1;
	  cells[i].posCells[iDim] = getID(cells[i].coords + tmpvec);
	  cells[i].negCells[iDim] = getID(cells[i].coords - tmpvec);
	}
    }
}

long 
CGCells::getID(CVector<long> coords) const
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
CGCells::getID(CVector<> pos) const
{
  Sim->Dynamics.BCs().setPBC(pos);
  CVector<long> temp;
  
  for (int iDim = 0; iDim < NDIM; iDim++)
    temp[iDim] = (long) ( (pos[iDim] + 0.5 * Sim->aspectRatio[iDim]) / cellLatticeWidth[iDim]);
  
  return getID(temp);
}
