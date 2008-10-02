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

#include "cells.hpp"
#include "../dynamics/liouvillean/CompressionL.hpp"
#include "../dynamics/interactions/intEvent.hpp"
#include "../extcode/threadpool.hpp"
#include "../simulation/particle.hpp"
#include "../dynamics/dynamics.hpp"
#include "../dynamics/liouvillean/liouvillean.hpp"
#include "../dynamics/BC/BC.hpp"
#include "../base/is_simdata.hpp"
#include "../base/is_base.hpp"
#include "../dynamics/globals/globEvent.hpp"
#include "../dynamics/systems/system.hpp"
#include <cmath> //for huge val
#include <boost/lexical_cast.hpp>
#include "../extcode/xmlParser.h"

CSCells::CSCells(const DYNAMO::SimData* Sim, const char *nom):
  CScheduler(Sim, nom),
  cellCount(0),
  cellDimension(1.0),
  lambda(0.9), //Default to higher overlap
  NCells(0)
{
  I_cout() << "Cellular algorithm loaded";
}

void
CSCells::addCells(Iflt maxdiam, bool limitCells)
{
  cells.clear();
  partCellData.resize(Sim->lN); //Location data for particles

  NCells = 1;
  cellCount = CVector<long>(0);

  for (int iDim = 0; iDim < NDIM; iDim++)
    {
      cellCount[iDim] = static_cast<long>(Sim->aspectRatio[iDim] / maxdiam);
      
      if (cellCount[iDim] < 3)
	I_throw() << "Not enough cells in " << static_cast<char>('x'+iDim) << " dimension, need 3+";

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

  for (long id = 0; id < NCells; id++)
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
CSCells::init_cells()
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
  
  //Tell the cells about their neighbours
  for (int i = 0; i < NCells; ++i) 
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
  
  //This is horribly slow, but arbitrary dimension
  BOOST_FOREACH(int dim, unlinktasklist)
    {
      I_cout() << "Unlinking the cells in the " << dim << " dimension";
      std::list<int> list1;
      std::list<int> list2;
      
      BOOST_FOREACH(cellStruct& cell, cells)
	{
	  if (cell.coords[dim] == 0)
	    list1.push_back(getID(cell.coords));
	  
	  if (cell.coords[dim] == cellCount[dim]-1)
	    list2.push_back(getID(cell.coords));
	}
      
      //Perform the removing
      BOOST_FOREACH(int cell1, list1)
	BOOST_FOREACH(int cell2, list2)
	{
	  cells[cell1].neighbours.erase
	    (std::remove(cells[cell1].neighbours.begin(), 
			 cells[cell1].neighbours.end(), cell2), 
	     cells[cell1].neighbours.end());
	  
	  cells[cell2].neighbours.erase
	    (std::remove(cells[cell2].neighbours.begin(), 
			 cells[cell2].neighbours.end(), cell1), 
	     cells[cell2].neighbours.end());
	}
    }
}

void 
CSCells::link_LE_cells()
{
  //Boundaries are LE, modify the Cell list
  I_cout() << "Linking cells required for LE BC";
  
  //build a list of neighbours for the bottom cells
  //Holds the displacement in each dimension, the unit is cells!
  CVector<long> displacement;
  
  displacement[1] = -1;

  std::list<CVector<long> > neighbourVectors;
  
  for (displacement[0] = 0; displacement[0] < cellCount[0]; ++displacement[0])
    for (displacement[2] = -1; displacement[2] < 2; ++displacement[2])
      //Add the current vector to the list
      neighbourVectors.push_back (displacement);
  
  for (int i = 0; i < cellCount[0]; i++)
    for (int j = 0; j < cellCount[2]; j++)
      {
	CVector<long> currentCell;
	currentCell[0] = i;
	currentCell[1] = 0;
	currentCell[2] = j;
	
	BOOST_FOREACH(CVector<long>& disp, neighbourVectors)
	  {
	    //Positive direction
	    long currentID = getID(currentCell),
	      oppositeID = getID(currentCell + disp);

	    //Check its not already in the list
	    if (std::find(cells[currentID].neighbours.begin(), 
			  cells[currentID].neighbours.end(), 
			  oppositeID) == cells[currentID].neighbours.end())
	      cells[currentID].neighbours.push_back(oppositeID);


	    //Add the mirror neighbour
	    if (std::find(cells[oppositeID].neighbours.begin(), 
			  cells[oppositeID].neighbours.end(), 
			  currentID) == cells[oppositeID].neighbours.end())
	      cells[oppositeID].neighbours.push_back(currentID);
	  }
      }
}

long 
CSCells::getID(CVector<long> coords) const
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
CSCells::getCoordsFromID(unsigned long i) const
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
CSCells::getID(CVector<> pos) const
{
  Sim->Dynamics.BCs().setPBC(pos);
  CVector<long> temp;
  
  for (int iDim = 0; iDim < NDIM; iDim++)
    temp[iDim] = (long) ( (pos[iDim] + 0.5 * Sim->aspectRatio[iDim]) / cellLatticeWidth[iDim]);
  
  return getID(temp);
}

void 
CSCells::operator<<(const XMLNode& XML)
{
  if (XML.isAttributeSet("Unlink"))
    for (int i = 0; i < XML.nAttribute(); i++)
      {
	try {
	  const char * tmp = XML.getAttribute("Unlink",i);
	  unlinktasklist.push_back(boost::lexical_cast<int>(tmp));
	} catch (...)
	  {}
      }

    if (XML.isAttributeSet("lambda"))
      {
	try {
	  
	  lambda = boost::lexical_cast<Iflt>
	    (XML.getAttribute("lambda"));
	  
	}
	catch(...)
	  {
	    I_throw() << "Could not load the lambda value in cellular scheduler";
	  }

	if (lambda < 0.0 || lambda > 1.0)
	  I_throw() << "Lambda out of bounds [0,1), lambda = " << lambda;
      }
}

void 
CSCells::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("lambda") << lambda;
  
  BOOST_FOREACH(int a, unlinktasklist)
    XML << xmlw::attr("Unlink") << a;
}
