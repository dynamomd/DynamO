/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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
#include <dynamo/globals/neighbourList.hpp>
#include <dynamo/particle.hpp>
#include <magnet/math/morton_number.hpp>
#include <boost/unordered_map.hpp>
#include <vector>

namespace dynamo {
  /*! \brief A regular cell neighbour list implementation.
    
    This neighbour list is the main neighbour list implemenetation for
    dynamo. It uses a regular grid of cells into which the particles
    are sorted to accelerate calculating the neighbourhood of a single
    particle.

    There are several "unusual" properties of this neighbour list
    which are used to optimise its behaviour.
    
    Although the neighbour list is a regular grid of cells, each cell
    overlaps with its neighbours. This means that if you cross from
    one cell into another, you enter the other cell some finite
    distance from the cells border. This helps remove "rattling"
    events where particles rapidly pass between two cells. It also
    helps prevent local event objects (like walls) falling numerically
    outside all cells when the wall lies on the cell border.

    The second property is that the contents of each cell is stored as
    a std::vector. In theory, a linked list is far more memory
    efficient however, the vector is much more cache friendly and can
    boost performance by 50% in cases where the cell has multiple
    particles inside of it.
   */
  class GCells: public GNeighbourList
  {
  public:
    GCells(const magnet::xml::Node&, dynamo::SimData*);
    GCells(SimData*, const std::string&, size_t overlink = 1);

    virtual ~GCells() {}

    virtual GlobalEvent getEvent(const Particle &) const;

    virtual void runEvent(Particle&, const double) const;

    virtual void initialise(size_t);

    virtual void reinitialise();

    virtual void getParticleNeighbourhood(const Particle&, 
					  const nbHoodFunc&) const;

    virtual void getParticleNeighbourhood(const Vector&, 
					  const nbHoodFunc2&) const;

    virtual void getLocalNeighbourhood(const Particle&, 
				       const nbHoodFunc&) const;
    
    virtual void operator<<(const magnet::xml::Node&);

    Vector getCellDimensions() const 
    { return cellDimension; }

    virtual double getMaxSupportedInteractionLength() const;

  protected:
    size_t cellCount[3];
    magnet::math::DilatedInteger<3> dilatedCellMax[3];
    Vector cellDimension;
    Vector cellLatticeWidth;
    Vector cellOffset;

    double _oversizeCells;
    size_t NCells;
    size_t overlink;

    boost::signals2::scoped_connection _particleAdded;
    boost::signals2::scoped_connection _particleRemoved;

    //! \brief The start of the list of particles in each cell.
    mutable std::vector<std::vector<size_t> > list;

    //! \brief The local events in each cell.
    mutable std::vector<std::vector<size_t> > cells;

    /*! \brief The cell for a given particle.
      
      This container is an unordered map, so we only store the linked
      list for the particles actually inserted into this neighborlist.
     */
    mutable boost::unordered_map<size_t, size_t> partCellData;

    GCells(const GCells&);

    virtual void outputXML(magnet::xml::XmlStream&) const;
    void outputXML(magnet::xml::XmlStream&, const std::string&) const;

    magnet::math::MortonNumber<3> getCellID(Vector) const;

    void addCells(double);

    inline Vector calcPosition(const magnet::math::MortonNumber<3>& coords,
			       const Particle& part) const;

    Vector calcPosition(const magnet::math::MortonNumber<3>& coords) const;

    void addLocalEvents();

    inline void addToCell(size_t ID)
    { addToCell(ID, getCellID(Sim->particleList[ID].getPosition()).getMortonNum()); }

    inline void addToCell(size_t ID, size_t cellID) const
    {
      list[cellID].push_back(ID);
      partCellData[ID] = cellID;
    }
  
    inline void removeFromCell(size_t ID) const
    { removeFromCellwIt(ID, partCellData.find(ID)); }

    inline void removeFromCellwIt(size_t ID, boost::unordered_map<size_t, size_t>::iterator it) const
    {
#ifdef DYNAMO_DEBUG
    if (it == partCellData.end())
      M_throw() << "Could not find the particle's cell data";
#endif
      const size_t cellID = it->second;
      //Erase the cell data
      partCellData.quick_erase(it);

      std::vector<std::vector<size_t> >::iterator listIt = list.begin() + cellID;
#ifdef DYNAMO_DEBUG
    if (listIt == list.end())
      M_throw() << "Could not find the particle's cell data";
#endif

      std::vector<size_t>::iterator pit = std::find(listIt->begin(), listIt->end(), ID);      

#ifdef DYNAMO_DEBUG
      if (pit == listIt->end())
	M_throw() << "Removing a particle (ID=" << ID << ") which is not in a cell";
#endif

      *pit = listIt->back();
      listIt->pop_back();
    }
  };
}
