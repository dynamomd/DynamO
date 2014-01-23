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
#include <map>

#include <Judy.h>

namespace dynamo {
  namespace detail {
    namespace {
      ::std::size_t
      hash_combine(const ::std::size_t hash1, const ::std::size_t hash2)
      {
       return hash1 ^ (hash2 + 0x9e3779b9 + (hash1 << 6) + (hash1 >> 2));
      }
    }

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
class CaptureMap
{
public:
  typedef std::pair<size_t, size_t> key_type;
  typedef Word_t mapped_type;
  typedef std::pair<key_type, mapped_type> value_type;

  class const_iterator : public std::iterator<std::input_iterator_tag, value_type> {
  public:
    const_iterator(const CaptureMap& container, size_t idx):_container(container), _idx(idx)  {}
    const_iterator& operator++() { ++_idx; return *this; }
    value_type operator*() const { return _container[_idx]; }
    bool operator!=(const const_iterator& o) const { return !(*this == o); }
    bool operator==(const const_iterator& o) const { return _idx == o._idx;}

  private:
    const CaptureMap& _container;
    size_t _idx;
  };

  const_iterator begin() const { return const_iterator(*this, 0); }
  const_iterator end() const { return const_iterator(*this, size()); }

private:
  static const size_t half_shift = sizeof(Word_t) * 4;
  static const size_t mask = (size_t(1) << half_shift) - 1;
      
  Pvoid_t _array = (Pvoid_t) NULL;
  size_t _count = 0;
      
  static inline Word_t keyToID(key_type key) {
    return Word_t((std::min(key.first, key.second) << half_shift) | std::max(key.first, key.second)); 
  }
  
  static inline key_type IDtoKey(Word_t key) {
    return key_type((size_t(key) >> half_shift) & mask,  size_t(key) & mask);
  }

public:

  ~CaptureMap() { clear(); }

  void clear() { JudyLFreeArray(&_array, NULL); }
  /*!\brief This proxy is used to double check if an assignment of
    zero is done, and delete the entry if it is. */
  struct EntryProxy {
  public:
    EntryProxy(CaptureMap& container, const key_type& key):
      _container(container), _ID(keyToID(key)) {}

    operator const mapped_type() const {
      mapped_type* PValue = NULL;
      PValue = (mapped_type*)JudyLGet(_container._array, _ID, NULL);
      return (PValue == NULL) ? 0 : *PValue;
    }
	
    EntryProxy& operator=(mapped_type newval) {
      if (newval == 0)
	_container._count -= JudyLDel(&(_container._array), _ID, NULL);
      else
	{
	  mapped_type& Value = *(mapped_type*)JudyLIns(&(_container._array), _ID, NULL);
	  if (Value == 0) ++_container._count;
	  Value = newval;
	}
      return *this;
    }
	
  private:
    CaptureMap& _container;
    const Word_t _ID;
  };

  size_t size() const { return _count; }
  
  value_type operator[](const size_t i) const {
    Word_t ID;
    mapped_type val = *(mapped_type*)JudyLByCount(_array, i+1, &ID, NULL);
    return value_type(IDtoKey(ID), val);
  }
      
  /*! \brief This non-const array access operator uses EntryProxy
    to check if any values assigned are zero so they may be deleted. */
  EntryProxy operator[](const key_type& key) {
    return EntryProxy(*this, key); 
  }

  mapped_type operator[](const key_type& key) const {
    mapped_type* PValue = NULL;
    PValue = (mapped_type*)JudyLGet(_array, keyToID(key), NULL);
    return (PValue == NULL) ? 0 : *PValue;
  }
};

    struct CaptureMapKey: public std::vector<CaptureMap::value_type>
    {
      typedef std::vector<CaptureMap::value_type> Container;
      CaptureMapKey(const CaptureMap& map):
       Container(map.begin(), map.end())
      {}

      std::size_t hash() const {
	std::size_t hash(0);
	for (const Container::value_type& val : *this)
	  hash = hash_combine(hash, hash_combine(val.first.first, hash_combine(val.first.second, val.second)));
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
      return Map::operator[](Map::key_type(p1, p2));
    }

    //! \brief A test if two particles are captured
    size_t isCaptured(const size_t p1, const size_t p2) const {
      return Map::operator[](Map::key_type(p1, p2));
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
      if (Map::operator[](Map::key_type(p1.getID(), p2.getID())) != 0)
	M_throw() << "Adding a particle while its already added!";
#endif
      Map::operator[](Map::key_type(p1.getID(), p2.getID())) = 1;
    }
  
    //! \brief Remove a pair of particles to the capture map.
    void remove(const Particle& p1, const Particle& p2)
    {
#ifdef DYNAMO_DEBUG
      if (Map::operator[](Map::key_type(p1.getID(), p2.getID())) == 0)
	M_throw() << "Deleting a particle while its already gone!";
#endif
      Map::operator[](Map::key_type(p1.getID(), p2.getID())) = 0;
    } 
  };
}
