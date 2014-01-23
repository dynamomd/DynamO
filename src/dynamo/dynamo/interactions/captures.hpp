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

    typedef magnet::containers::JudyPairMap CaptureMap;

    struct CaptureMapKey: public std::vector<CaptureMap::value_type>
    {
      typedef std::vector<CaptureMap::value_type> Container;
      CaptureMapKey(const CaptureMap& map):
       Container(map.begin(), map.end()) {}

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
