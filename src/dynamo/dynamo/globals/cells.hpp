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
#ifdef DYNAMO_JUDY
# include <magnet/containers/judy.hpp>
#endif
#include <magnet/containers/vector_set.hpp>
#include <magnet/containers/multimaps.hpp>
#include <magnet/containers/ordering.hpp>
#include <unordered_map>
#include <vector>

namespace dynamo {
  namespace detail {
    /*! \brief A container for storing the cell contents (and which
        particle is in which cell).
	
	\tparam CellList A type which gives a multimap-style container
	mapping from cell ids to particle ids. Examples include
	SetCellList<JudySet<uint64_t>>, and
	VectorSetCellList<JudySet<size_t>> although
	Vector_Multimap<VectorSet<size_t>> appears to be the most
	performant.

	\tparam Map A map container which links particle IDs to cell
	IDs. Examples include std::unordered_map<size_t, size_t> but
	JudyMap<size_t, size_t> appears to be the best.
     */
    template<typename CellList, typename Map>
    class CellParticleList {
      CellList _cellcontents;
      Map _particleCell;
     
    public:
      void add(size_t cell, size_t particle) {
	_cellcontents.insert(cell, particle);
	_particleCell[particle] = cell;
      }
      
      void remove(size_t cell, size_t particle) {
	_cellcontents.erase(cell, particle);
	_particleCell.erase(particle);
      }

      void moveTo(size_t oldcell, size_t newcell, size_t particle) {
	_cellcontents.erase(oldcell, particle);
	_cellcontents.insert(newcell, particle);
	_particleCell[particle] = newcell;
      }

      typename CellList::RangeType getCellContents(const size_t cellID) {
	return _cellcontents.getKeyContents(cellID);
      }

      size_t getCellID(const size_t particle) const {
#ifdef MAGNET_DEBUG
	if (_particleCell.find(particle) == _particleCell.end())
	  M_throw() << "Could not find the cell for particle " << particle << " during cell look-up";
#endif
	return _particleCell.find(particle)->second;
      }

      void resize(size_t cellcount, size_t N) { 
	_cellcontents.resize(cellcount); 
      }

      size_t size() const { return _particleCell.size(); }
      void clear() { _particleCell.clear(); _cellcontents.clear(); }
    };
  }

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
    events where particles rapidly pass between two cells.

    The second property is that the contents of each cell is stored as
    a std::vector. In theory, a linked list is far more memory
    efficient however, the vector is much more cache friendly and can
    boost performance by 50% in cases where the cell has multiple
    particles inside of it.
   */
  class GCells: public GNeighbourList
  {
  public:
    GCells(const magnet::xml::Node&, dynamo::Simulation*);
    GCells(Simulation*, const std::string&);

    virtual ~GCells() {}

    virtual GlobalEvent getEvent(const Particle &) const;

    virtual void runEvent(Particle&, const double) const;

    virtual void initialise(size_t);

    virtual void reinitialise();

    void getParticleNeighbours(const Particle&, std::vector<size_t>&) const;
    void getParticleNeighbours(const Vector&, std::vector<size_t>&) const;
    
    virtual void operator<<(const magnet::xml::Node&);

    Vector getCellDimensions() const 
    { return _cellDimension; }

    virtual double getMaxSupportedInteractionLength() const;

    void setConfigOutput(bool val) { _inConfig = val; }

  protected:
    virtual void getParticleNeighbours(const std::array<size_t, 3>&, std::vector<size_t>&) const;

    typedef magnet::containers::RowMajorOrdering<3> Ordering;
    Ordering _ordering;

    Vector _cellDimension;
    Vector _cellLatticeWidth;
    Vector _cellOffset;

    bool _inConfig;
    double _oversizeCells;
    size_t overlink;

#ifdef DYNAMO_JUDY
    mutable detail::CellParticleList<magnet::containers::Vector_Multimap<magnet::containers::VectorSet<size_t>>, 
				     magnet::containers::JudyMap<size_t, size_t>> _cellData;
#else
    mutable detail::CellParticleList<magnet::containers::Vector_Multimap<magnet::containers::VectorSet<size_t>>, 
				     std::unordered_map<size_t, size_t> > _cellData;
#endif
    GCells(const GCells&);

    virtual void outputXML(magnet::xml::XmlStream&) const;

    std::array<size_t, 3> getCellCoords(Vector) const;

    void addCells(double);
    void buildCells();

    Vector calcPosition(const size_t cellIndex, const Particle& part) const { return calcPosition(_ordering.toCoord(cellIndex), part);}
    Vector calcPosition(const std::array<size_t, 3>& coords, const Particle& part) const ;
    Vector calcPosition(const size_t cellIndex) const { return calcPosition(_ordering.toCoord(cellIndex));}
    Vector calcPosition(const std::array<size_t, 3>& coords) const;
  };
}
