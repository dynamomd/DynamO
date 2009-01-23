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

#include "gcells2.hpp"
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

CGCells2::CGCells2(DYNAMO::SimData* nSim, const std::string& name):
  CGNeighbourList(nSim, "GlobalCellularEvent2"),
  cellCount(0),
  cellDimension(1.0),
  lambda(0.9), //Default to higher overlap
  NCells(0)

{
  globName = name;
  I_cout() << "Cells Loaded";
}

CGCells2::CGCells2(const XMLNode &XML, DYNAMO::SimData* ptrSim):
  CGNeighbourList(ptrSim, "GlobalCellularEvent"),
  cellCount(0),
  cellDimension(1.0),
  lambda(0.9), //Default to higher overlap
  NCells(0)
{
  operator<<(XML);

  I_cout() << "Cells Loaded";
}

CGCells2::CGCells2(DYNAMO::SimData* ptrSim, const char* nom, void*):
  CGNeighbourList(ptrSim, nom),
  cellCount(0),
  cellDimension(1.0),
  lambda(0.9), //Default to higher overlap
  NCells(0)
{}

void 
CGCells2::operator<<(const XMLNode& XML)
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
      D_throw() << "Error loading CGCells2";
    }
  
  if (lambda < 0.0 || lambda > 1.0)
    D_throw() << "Lambda out of bounds [0,1), lambda = " << lambda;
}

void 
CGCells2::setLambda(const Iflt& nL)
{
  lambda = nL;
}

CGlobEvent 
CGCells2::getEvent(const CParticle& part) const
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
CGCells2::runEvent(const CParticle& part) const
{
  Sim->Dynamics.Liouvillean().updateParticle(part);
  
  size_t oldCell(partCellData[part.getID()].cell);

  //Determine the cell transition direction, its saved
  size_t cellDirection(Sim->Dynamics.Liouvillean().
		       getSquareCellCollision3
		       (part, cells[oldCell].origin, 
			cellDimension));

  CVector<int> newcoords = cells[oldCell].coords;

  newcoords[cellDirection]
    += ((part.getVelocity()[cellDirection] > 0) * 2) - 1;

  size_t endCell = getCellID(newcoords);

  newcoords[cellDirection]
    += ((part.getVelocity()[cellDirection] > 0) * 2) - 1;

  size_t inPosition 
    = newcoords[cellDirection] % cellCount[cellDirection]
    + (newcoords[cellDirection] < 0) * cellCount[cellDirection];

  //Debug section
#ifdef DYNAMO_WallCollDebug
  {      
    CVector<int> tmp = cells[partCellData[part.getID()].cell].coords;
    CVector<int> tmp2 = cells[endCell].coords;
    
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
  //Holds the displacement in each dimension, the unit is cells!
  CVector<int>displacement(-1);

  //This loop iterates through each neighbour position
  for (int iter = 0; iter < ctime_pow<3,NDIM>::result; ++iter)
    {
      //Add the current vector to the list
      const int nb = getCellID(cells[endCell].coords + displacement);

      if (size_t(cells[nb].coords[cellDirection]) == inPosition)
	for (int next = cells[nb].list; next != -1; 
	     next = partCellData[next].next)
	  BOOST_FOREACH(const nbHoodSlot& nbs, sigNewNeighbourNotify)
	    nbs.second(part, next);

      //Now update the displacement vector
      ++displacement[0];

      for (int iDim = 1; iDim < NDIM; ++iDim)
	if (displacement[iDim - 1] == 2)
	  {
	    displacement[iDim - 1] = -1;
	    displacement[iDim] += 1;
	  }
    }

  //Tell about the new locals
  BOOST_FOREACH(const size_t& lID, cells[endCell].locals)
    BOOST_FOREACH(const nbHoodSlot& nbs, sigNewLocalNotify)
    nbs.second(part, lID);
  
  //Push the next virtual event, this is the reason the scheduler
  //doesn't need a second callback
  Sim->ptrScheduler->pushEvent(part, getEvent(part));
  Sim->ptrScheduler->sort(part);

  BOOST_FOREACH(const nbHoodSlot& nbs, sigCellChangeNotify)
    nbs.second(part, oldCell);
  
  //This doesn't stream the system as its a virtual event
}

void 
CGCells2::initialise(size_t nID)
{
  ID=nID;
  
  reinitialise(Sim->Dynamics.getLongestInteraction());
}

void
CGCells2::reinitialise(const Iflt& maxdiam)
{
  I_cout() << "Reinitialising on collision " << Sim->lNColl;

  //Create the cells
  addCells(maxdiam, false);

  addLocalEvents();

  BOOST_FOREACH(const initSlot& nbs, sigReInitNotify)
    nbs.second();
}

void
CGCells2::outputXML(xmlw::XmlStream& XML) const
{
  //If you add anything here it also needs to go in gListAndCells.cpp too
  XML << xmlw::attr("Type") << "Cells2"
      << xmlw::attr("Lambda") << lambda
      << xmlw::attr("Name") << globName;
}

void
CGCells2::addCells(Iflt maxdiam, bool limitCells)
{
  cells.clear();
  partCellData.resize(Sim->lN); //Location data for particles

  NCells = 1;
  cellCount = CVector<int>(0);

  for (int iDim = 0; iDim < NDIM; iDim++)
    {
      cellCount[iDim] = int(Sim->aspectRatio[iDim] / maxdiam);
      
      if (cellCount[iDim] < 3)
	D_throw() << "Not enough cells in " << char('x'+iDim) << " dimension, need 3+";

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
}

void 
CGCells2::addLocalEvents()
{
  BOOST_FOREACH(cellStruct& cell, cells)
    {
      cell.locals.clear();

      BOOST_FOREACH(const smrtPlugPtr<CLocal>& local, Sim->Dynamics.getLocals())
	if (local->isInCell(cell.origin, cellDimension))
	  cell.locals.push_back(local->getID());
    }
}

size_t
CGCells2::getCellID(const CVector<int>& coords) const
{
  //PBC for vectors

  size_t id = 0;
  size_t pow = 1;

  for (size_t iDim = 0; iDim < NDIM - 1; ++iDim)
    {
      id += ((coords[iDim] % cellCount[iDim]) + (coords[iDim] < 0) * cellCount[iDim]) * pow;
      pow *= cellCount[iDim];
    }
  
  //Previous Loop runs to NDIM - 1 to stop an invalid memory read and
  //useless *=
  return id + pow * ((coords[NDIM - 1] % cellCount[NDIM-1]) 
		     + (coords[NDIM-1] < 0) * cellCount[NDIM-1]);
}

CVector<int> 
CGCells2::getCoordsFromID(size_t i) const
{
  CVector<int> tmp;
  i = i % NCells; //PBC's for ID
  
  tmp[0] = i % cellCount[0];
  i /= cellCount[0];
  tmp[1] = i % cellCount[1];
  i /= cellCount[1];
  tmp[2] = i % cellCount[2];
  return tmp;
}

size_t
CGCells2::getCellID(CVector<> pos) const
{
  Sim->Dynamics.BCs().setPBC(pos);
  CVector<int> temp;
  
  for (size_t iDim = 0; iDim < NDIM; iDim++)
    temp[iDim] = int((pos[iDim] + 0.5 * Sim->aspectRatio[iDim]) / cellLatticeWidth[iDim]);
  
  return getCellID(temp);
}


void 
CGCells2::getParticleNeighbourhood(const CParticle& part,
				  const nbHoodFunc& func) const
{
  CVector<int> coords
    (cells[partCellData[part.getID()].cell].coords - CVector<int>(1));

  const CVector<int> coordsorig(coords);

  //This loop iterates through each neighbour position
  BOOST_STATIC_ASSERT(NDIM==3);

  for (size_t iDim(0); iDim < 3; ++iDim)
    {
      for (size_t jDim(0); jDim < 3; ++jDim)
	{
	  for (size_t kDim(0); kDim < 3; ++kDim)
	    {
	      //Add the current vector to the list
	      const size_t& nb = getCellID(coords);
	    
	      for (int next = cells[nb].list;
		   next != -1; next = partCellData[next].next)
		if (next != int(part.getID()))
		  func(part, next);
	      ++coords[2];
	    }
	  ++coords[1];
	  coords[2] = coordsorig[2];
	}
      ++coords[0];
      coords[1] = coordsorig[1];
    }
}

void 
CGCells2::getParticleLocalNeighbourhood(const CParticle& part, 
				       const nbHoodFunc& func) const
{
  BOOST_FOREACH(const size_t& id, cells[partCellData[part.getID()].cell].locals)
    func(part, id);
}
