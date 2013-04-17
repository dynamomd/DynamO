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
#include <boost/functional/hash.hpp>
#include <tr1/unordered_set>
#include <tr1/unordered_map>
#include <vector>

namespace dynamo {
  /*! \brief A general interface for \ref Interaction classes with
     states for the particle pairs.
   
    This class is a general interface to Interaction classes that
    allow particles to "capture" each other and store some state. The
    state might be the internal energy between particle pairs (e.g.,
    \ref ISquareWell), or it might be used to track if the particles
    are within each others bounding sphere (e.g., \ref ILines).
  */
  class ICapture: public Interaction
  {
  protected:
    /*! \brief A key used to represent two particles.
     
      This key sorts the particle ID's into ascending order. This way
      the keys can be compared and symmetric keys will compare equal.
      \code assert(cMapKey(a,b) == cMapKey(b,a)); \endcode
     */
    struct cMapKey: public std::pair<size_t,size_t>
    {
      inline cMapKey() {}
    
      //The ContactMap output plugin relies that the ID pair is sorted
      //min,max
      inline cMapKey(const size_t& a, const size_t& b):
	std::pair<size_t,size_t>(std::min(a, b), std::max(a, b))
      {
#ifdef DYNAMO_DEBUG
	if (a == b) M_throw() << "Particle ID's should not be equal!";
#endif
      }
    };
    
  public:
    typedef std::tr1::unordered_map<cMapKey, size_t, boost::hash<cMapKey> > captureMapType;
    typedef captureMapType::key_type key_type;
    typedef captureMapType::mapped_type mapped_type;
    typedef captureMapType::value_type value_type;
    typedef captureMapType::iterator cmap_it;
    typedef captureMapType::const_iterator const_cmap_it;

    ICapture(dynamo::Simulation* sim, IDPairRange* range): Interaction(sim, range), noXmlLoad(true) {}

    //! \brief Returns the number of particles that are captured in some way
    size_t getTotalCaptureCount() const { return captureMap.size(); }
  
    //! \brief A test if two particles are captured
    size_t isCaptured(const Particle& p1, const Particle& p2) const {
      return isCaptured(p1.getID(), p2.getID()); 
    }

    //! \brief A test if two particles are captured
    size_t isCaptured(const size_t p1, const size_t p2) const {
      cmap_it it = captureMap.find(cMapKey(p1, p2));
      return (it == captureMap.end()) ? 0 : (it->second); 
    }

    //! \brief Returns the total internal energy stored in this Interaction.
    virtual double getInternalEnergy() const = 0;

    /*! \brief This function tells an uninitialised capture map to
        forget the data loaded from the xml file.
     */
    void forgetXMLCaptureMap() { noXmlLoad = true; }

    void initCaptureMap();

    virtual size_t captureTest(const Particle&, const Particle&) const = 0;

    void clear() const { captureMap.clear(); }
    
    const captureMapType& getMap() const {
      return captureMap;
    }
    
  protected:  
    bool noXmlLoad;
    mutable captureMapType captureMap;

    void loadCaptureMap(const magnet::xml::Node&);

    void outputCaptureMap(magnet::xml::XmlStream&) const;

    virtual size_t validateState(bool textoutput = true, size_t max_reports = std::numeric_limits<size_t>::max()) const;

    virtual void testAddToCaptureMap(const Particle& p1, const size_t& p2) const;

    //! \brief Add a pair of particles to the capture map.
    void addToCaptureMap(const Particle& p1, const Particle& p2) const {
#ifdef DYNAMO_DEBUG
      if (captureMap.find(cMapKey(p1.getID(), p2.getID())) != captureMap.end())
	M_throw() << "Adding a particle while its already added!";
#endif
      captureMap[cMapKey(p1.getID(), p2.getID())] = 1;
     }
  
    //! \brief Remove a pair of particles to the capture map.
    void removeFromCaptureMap(const Particle& p1, const Particle& p2) const
    {
#ifdef DYNAMO_DEBUG
      if (captureMap.find(cMapKey(p1.getID(), p2.getID())) == captureMap.end())
	M_throw() << "Deleting a particle while its already gone!";
#endif
      captureMap.erase(cMapKey(p1.getID(), p2.getID()));
    } 

    inline void delFromCaptureMap(const Particle& p1, const Particle& p2) const
    {
#ifdef DYNAMO_DEBUG
      if (captureMap.find(cMapKey(p1.getID(), p2.getID())) == captureMap.end())
	M_throw() << "Deleting a particle while its already gone!";
#endif 
      captureMap.erase(cMapKey(p1.getID(), p2.getID()));
    }

    inline cmap_it getCMap_it(const Particle& p1, const Particle& p2) const
    { return captureMap.find(cMapKey(p1.getID(), p2.getID())); }
  };
}
