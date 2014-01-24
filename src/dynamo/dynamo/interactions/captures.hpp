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

#include <dynamo/particle.hpp>
#include <dynamo/interactions/interaction.hpp>
#include <magnet/exception.hpp>
#include <magnet/containers/judy.hpp>
#include <map>

namespace dynamo {
  namespace detail {
    namespace {
      ::std::size_t
      hash_combine(const ::std::size_t hash1, const ::std::size_t hash2)
      {
	return hash1 ^ (hash2 + 0x9e3779b9 + (hash1 << 6) + (hash1 >> 2));
      }
    }
    
    static_assert(sizeof(size_t) == 8, "DynamO now requires 64bit systems as it uses Judy arrays for capture maps");
    /*! \brief A key used to represent a pair of two particles.
      
      This key sorts the particle ID's into ascending order. This way
      the keys can be compared and symmetric keys will compare equal.
      \code assert(cMapKey(a,b) == cMapKey(b,a)); \endcode
    */
    struct PairKey {
      inline PairKey(const PairKey& o) { _key = o._key; }
      inline PairKey(uint64_t val) { _key = val; }
      inline PairKey(const size_t p1, size_t p2):
	first(std::min(p1, p2)), second(std::max(p1, p2))
      {
#ifdef DYNAMO_DEBUG
	if (first == second) M_throw() << "Particle ID's should not be equal!";
#endif
      }
      
      operator uint64_t() const { return _key; }      
      union {
	struct {uint32_t first; uint32_t second; };
	uint64_t _key;
      };
    };

    /*!\brief This is a container that stores a single size_t
      identified by a pair of particles.
       
      To efficiently store the state of all possible particle
      pairings, a map is used and entries are only stored if the
      state is non-zero. Although an unordered_map may be slightly
      faster, a map is used to allow comparison of CaptureMaps and a
      rapid hashing if CaptureMaps are to be used as an index of the
      simulation state.
       
      To facilitate the storage only if non-zero behaviour, the array
      access operator is overloaded to automatically return a size_t
      0 for any entry which is missing. It also returns a proxy which
      deletes entries when they are set to 0.
    */
    class CaptureMap: public magnet::containers::JudyMap<PairKey, size_t>
    {
      typedef magnet::containers::JudyMap<PairKey, size_t> Container;
    public:
      /*!\brief This proxy is used to double check if an assignment of
	zero is done, and delete the entry if it is. */
      struct EntryProxy {
      public:
	EntryProxy(Container& container, const PairKey& key):
	  _container(container), _key(key) {}

	operator const size_t() const {
	  const auto it (_container.find(_key));
	  return (it == _container.end()) ? 0 : (it->second);
	}
	
	EntryProxy& operator=(size_t newval) {
	  if (newval == 0)
	    _container.erase(_key);
	  else
	    _container[_key] = newval;

	  return *this;
	}
	
      private:
	Container& _container;
	const PairKey _key;
      };
      
      /*! \brief This non-const array access operator uses EntryProxy
	to check if any values assigned are zero so they may be deleted. */
      EntryProxy operator[](const PairKey& key) {
	return EntryProxy(*this, key); 
      }

      /*! \brief A simple const array access operator which returns 0
	if the entry is missing. */
      size_t operator[](const PairKey& key) const {
	Container::const_iterator it = Container::find(key);
	return (it == Container::end()) ? 0 : (it->second);
      }
    };

    struct CaptureMapKey: public std::vector<CaptureMap::value_type>
    {
      typedef std::vector<CaptureMap::value_type> Container;
      CaptureMapKey(const CaptureMap& map):
      Container(map.begin(), map.end()) {}

      std::size_t hash() const {
	std::size_t hash(0);
	for (const Container::value_type& val : *this)
	  hash = hash_combine(hash, hash_combine(size_t(val.first), size_t(val.second)));
	return hash;
      }
    };

    /*! \brief A functor to allow the storage of CaptureMapKey types
      in unordered containers. */
    struct CaptureMapKeyHash {
      std::size_t operator() (const CaptureMapKey& map) const 
      { return map.hash(); }
    };
  }

  /*! \brief A general interface for \ref Interaction classes with
    states for the particle pairs.
   
    This class is a general interface to Interaction classes that
    allow particles to "capture" each other and store some state. The
    state might be the internal energy between particle pairs (e.g.,
    \ref ISquareWell), or it might be used to track if the particles
    are within each others bounding sphere (e.g., \ref ILines).
  */
  class ICapture: public Interaction, public detail::CaptureMap
  {
    typedef detail::CaptureMap Map;

  public:
    ICapture(dynamo::Simulation* sim, IDPairRange* range): Interaction(sim, range), _mapUninitialised(true) {}

    //! \brief A test if two particles are captured
    size_t isCaptured(const Particle& p1, const Particle& p2) const {
      return Map::operator[](detail::PairKey(p1, p2));
    }

    //! \brief A test if two particles are captured
    size_t isCaptured(const size_t p1, const size_t p2) const {
      Map::const_iterator it = Map::find(detail::PairKey(p1, p2));
      return (it == Map::end()) ? 0 : it->second;
    }

    /*! \brief This function tells an uninitialised capture map to
        forget the data loaded from the xml file.
     */
    void forgetMap() { _mapUninitialised = true; }

    void initCaptureMap();

    virtual size_t captureTest(const Particle&, const Particle&) const = 0;

  protected:  
    bool _mapUninitialised;

    void loadCaptureMap(const magnet::xml::Node&);

    void outputCaptureMap(magnet::xml::XmlStream&) const;

    virtual size_t validateState(bool textoutput = true, size_t max_reports = std::numeric_limits<size_t>::max()) const;

    virtual void testAddToCaptureMap(const Particle& p1, const size_t& p2);

    //! \brief Add a pair of particles to the capture map.
    void add(const Particle& p1, const Particle& p2) {
#ifdef DYNAMO_DEBUG
      if (Map::operator[](detail::PairKey(p1, p2)) != 0)
	M_throw() << "Adding a particle while its already added!";
#endif
      Map::operator[](detail::PairKey(p1, p2)) = 1;
    }
  
    //! \brief Remove a pair of particles to the capture map.
    void remove(const Particle& p1, const Particle& p2)
    {
#ifdef DYNAMO_DEBUG
      if (Map::operator[](detail::PairKey(p1, p2)) == 0)
	M_throw() << "Deleting a particle while its already gone!";
#endif
      Map::operator[](detail::PairKey(p1, p2)) = 0;
    } 
  };
}
