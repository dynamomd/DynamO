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

#ifndef CSCells_H
#define CSCells_H

#include "scheduler.hpp"
#include "../datatypes/vector.hpp"
#include <vector>
#include <list>
#include "../extcode/mathtemplates.hpp"

class CSCells : public CScheduler
{
 public:
  CSCells(const DYNAMO::SimData*, const char *);
  
  void init_cells();

  void link_LE_cells();
  
  virtual void operator<<(const XMLNode&);

  virtual void reinitialise(Iflt) = 0;

  CVector<> getCellDimensions() const 
  { return cellDimension; }

  inline void setLambda(const Iflt& nL) { lambda = nL; }
  inline Iflt getLambda() { return lambda; }

  inline void addUnlinkTask(const int& iDim) { unlinktasklist.push_back(iDim); }

 protected:
  virtual void outputXML(xmlw::XmlStream&) const;

  //Cell Numbering
  long getID(CVector<long>) const;
  CVector<long> getCoordsFromID(unsigned long) const; 
  long getID(CVector<>) const;

  //Cell initialisation
  //limit the cell count to 100?
  void addCells(Iflt, bool limitCells = true);

  //Variables
  CVector<long> cellCount;
  CVector<> cellDimension;
  CVector<> cellLatticeWidth;
  Iflt lambda;
  long NCells;
  

  std::list<int> unlinktasklist;  

  struct cellStruct
  {
    //Be smart about memory
    cellStruct():list(-1)
    { neighbours.reserve(ctime_pow<NDIM,NDIM>::result - 1); }

    std::vector<int> neighbours;
    int list;
    CVector<> origin;
    CVector<long> coords;
    size_t posCells[NDIM]; 
    size_t negCells[NDIM];
  };

  struct partCEntry
  {
    int prev;
    int next;
    int cell;
  };

  mutable std::vector<cellStruct> cells;
  mutable std::vector<partCEntry> partCellData;

  inline void addToCell(const int& ID, const int& cellID) const
  {
#ifdef DYNAMO_DEBUG
    if (cells.at(cellID).list != -1)
      partCellData.at(cells.at(cellID).list).prev = ID;
    
    partCellData.at(ID).next = cells.at(cellID).list;
    cells.at(cellID).list = ID;    
    partCellData.at(ID).prev = -1;
    partCellData.at(ID).cell = cellID;
# else
    if (cells[cellID].list != -1)
      partCellData[cells[cellID].list].prev = ID;
    
    partCellData[ID].next = cells[cellID].list;
    cells[cellID].list = ID;    
    partCellData[ID].prev = -1;
    partCellData[ID].cell = cellID;    
#endif
  }
  
  inline void removeFromCell(const int& ID) const
  {
    /* remove from linked list */    
    if (partCellData[ID].prev != -1)
      partCellData[partCellData[ID].prev].next = partCellData[ID].next;
    else
      cells[partCellData[ID].cell].list = partCellData[ID].next;
    
    if (partCellData[ID].next != -1)
      partCellData[partCellData[ID].next].prev = partCellData[ID].prev;

#ifdef DYNAMO_DEBUG
    partCellData[ID].cell = -1;
#endif
  }

};
#endif
