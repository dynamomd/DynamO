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

#pragma once
#include <dynamo/dynamics/globals/neighbourList.hpp>
#include <dynamo/simulation/particle.hpp>
#include <magnet/math/morton_number.hpp>
#include <vector>

class GCells: public GNeighbourList
{
public:
  GCells(const magnet::xml::Node&, dynamo::SimData*);
  GCells(dynamo::SimData*, const std::string&, size_t overlink = 1);

  virtual ~GCells() {}

  virtual Global* Clone() const { return new GCells(*this); }

  virtual GlobalEvent getEvent(const Particle &) const;

  virtual void runEvent(const Particle&, const double) const;

  virtual void initialise(size_t);

  virtual void reinitialise(const double&);

  virtual void getParticleNeighbourhood(const Particle&, 
					const nbHoodFunc&) const;

  virtual void getParticleLocalNeighbourhood(const Particle&, 
					     const nbHoodFunc&) const;
  
  virtual void operator<<(const magnet::xml::Node&);

  Vector getCellDimensions() const 
  { return cellDimension; }

  virtual double getMaxSupportedInteractionLength() const;

  virtual double getMaxInteractionLength() const;

protected:

  struct partCEntry
  {
    int prev;
    int next;
    int cell;
  };

  virtual void outputXML(magnet::xml::XmlStream&) const;
  void outputXML(magnet::xml::XmlStream&, const std::string&) const;

  magnet::math::MortonNumber<3> getCellID(Vector) const;

  void addCells(double);

  inline Vector calcPosition(const magnet::math::MortonNumber<3>& coords,
			     const Particle& part) const;

  Vector calcPosition(const magnet::math::MortonNumber<3>& coords) const;

  void addLocalEvents();

  //Variables
  size_t cellCount[3];
  magnet::math::DilatedInteger<3> dilatedCellMax[3];
  Vector cellDimension;
  Vector cellLatticeWidth;
  Vector cellOffset;

  double _oversizeCells;
  size_t NCells;
  size_t overlink;

  //! \brief The start of the list of particles in each cell.
  mutable std::vector<int> list;

  //! \brief The local events in each cell.
  mutable std::vector<std::vector<size_t> > cells;

  //! \brief The linked list of the particles in each cell.
  mutable std::vector<partCEntry> partCellData;

  inline void addToCell(const int& ID, const int& cellID) const
  {
#ifdef DYNAMO_DEBUG
    if (list.at(cellID) != -1)
      partCellData.at(list.at(cellID)).prev = ID;
    
    partCellData.at(ID).next = list.at(cellID);
    list.at(cellID) = ID;    
    partCellData.at(ID).prev = -1;
    partCellData.at(ID).cell = cellID;
# else
    if (list[cellID] != -1)
      partCellData[list[cellID]].prev = ID;
    
    partCellData[ID].next = list[cellID];
    list[cellID] = ID;    
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
      list[partCellData[ID].cell] = partCellData[ID].next;
    
    if (partCellData[ID].next != -1)
      partCellData[partCellData[ID].next].prev = partCellData[ID].prev;

#ifdef DYNAMO_DEBUG
    partCellData[ID].cell = -1;
#endif
  }

};
