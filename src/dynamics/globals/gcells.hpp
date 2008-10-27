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

#ifndef CGCells_HPP
#define CGCells_HPP

#include "global.hpp"
#include "../../datatypes/vector.hpp"
#include <vector>
#include <list>
#include "../../extcode/mathtemplates.hpp"
#include "../../simulation/particle.hpp"

class CGCells: public CGlobal
{
public:
  CGCells(const XMLNode &, const DYNAMO::SimData*);

  CGCells(const DYNAMO::SimData*);

  virtual ~CGCells() {}

  virtual CGlobal* Clone() const { return new CGCells(*this); }

  virtual CGlobEvent getEvent(const CParticle &) const;

  virtual CNParticleData runEvent(const CGlobEvent&) const;

  virtual void initialise(size_t);

  virtual void reinitialise(const Iflt&);
  
  void setLambda(const Iflt&);
  inline const Iflt& getLambda() const { return lambda; }
  inline const size_t& getID() const { return ID; }
  
  virtual void operator<<(const XMLNode&);

  CVector<> getCellDimensions() const 
  { return cellDimension; }

  inline size_t getLocalCellID(const CParticle& part) const
  { return partCellData[part.getID()].cell; }

  inline const std::vector<size_t>& getCellNeighbourHood(const CParticle& part) const
  {
    return cells[partCellData[part.getID()].cell].neighbours;
  }

  inline const int& getCellLocalParticles(const size_t& id) const
  { return cells[id].list; }

  struct partCEntry
  {
    int prev;
    int next;
    int cell;
  };

  inline const partCEntry& getParticleData(const size_t& id) const
  { return partCellData[id]; }

  struct cellStruct
  {
    //Be smart about memory
    cellStruct():list(-1)
    { neighbours.reserve(ctime_pow<NDIM,NDIM>::result); }

    std::vector<size_t> neighbours;
    std::vector<size_t> locals;    
    int list;
    CVector<> origin;
    CVector<long> coords;
    size_t posCells[NDIM]; 
    size_t negCells[NDIM];
  };

  inline const cellStruct& getParticleCellData(const CParticle& part) const
  { return cells[partCellData[part.getID()].cell]; }

protected:
  virtual void outputXML(xmlw::XmlStream&) const;

private:
  //Cell Numbering
  long getID(CVector<long>) const;
  CVector<long> getCoordsFromID(unsigned long) const; 
  long getID(CVector<>) const;

  void addCells(Iflt, bool limitCells = true);
  void init_cells();
  void addLocalEvents();

  //Variables
  CVector<long> cellCount;
  CVector<> cellDimension;
  CVector<> cellLatticeWidth;
  Iflt lambda;
  size_t NCells;

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
