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
#include <boost/tr1/unordered_set.hpp>
#include <boost/tr1/unordered_map.hpp>
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
  public:
    ICapture(dynamo::Simulation* sim, C2Range* range): Interaction(sim, range), noXmlLoad(true) {}

    //! \brief Returns the number of particles that are captured in some way
    virtual size_t getTotalCaptureCount() const = 0;
  
    //! \brief A test if two particles are captured
    virtual bool isCaptured(const Particle&, const Particle&) const = 0;

    //! \brief Returns the total internal energy stored in this Interaction.
    virtual double getInternalEnergy() const = 0;

    /*! \brief This function tells an uninitialised capture map to
        forget the data loaded from the xml file.
     */
    void forgetXMLCaptureMap() { noXmlLoad = true; }

    //! \brief Add a pair of particles to the capture map.
    virtual void addToCaptureMap(const Particle& p1, const Particle& p2) const = 0;

    void initCaptureMap();

    virtual void clear() const = 0;
    
  protected:
    bool noXmlLoad;

    virtual void testAddToCaptureMap(const Particle& p1, const size_t& p2) const = 0;

    /*! \brief A key used to represent two particles.
     
      This key sorts the particle ID's into ascending order. This way
      the keys can be compared and symmetric keys will compare equal.
      \code assert(cMapKey(a,b) == cMapKey(b,a)); \endcode
     */
    struct cMapKey: public std::pair<size_t,size_t>
    {
      inline cMapKey() {}
    
      inline cMapKey(const size_t& a, const size_t& b):
	std::pair<size_t,size_t>(std::min(a, b), std::max(a, b))
      {
#ifdef DYNAMO_DEBUG
	if (a == b) M_throw() << "Particle ID's should not be equal!";
#endif
      }
    };
  };

  /*! \brief This base class is for Interaction classes which only
    "capture" particle pairs in one state.
   
    There is only one state a pair of particles can be in, either
    captured or not.  This can be contrasted with IMultiCapture where
    a pair of particles may be in a range of captured states.
   */
  class ISingleCapture: public ICapture
  {
  public:
    ISingleCapture(dynamo::Simulation* sim, C2Range* range): ICapture(sim, range) {}

    size_t getTotalCaptureCount() const { return captureMap.size(); }
  
    virtual bool isCaptured(const Particle& p1, const Particle& p2) const
    { return captureMap.count(cMapKey(p1.getID(), p2.getID())); }

    virtual void clear() const { captureMap.clear(); }

  protected:

    mutable std::tr1::unordered_set<cMapKey > captureMap;

    /*! \brief Test if two particles should be "captured".
    
      This function should provide a test of the particles current
      position and velocity to determine if they're captured. Used in
      rebuilding the captureMap.
    */
    virtual bool captureTest(const Particle&, const Particle&) const = 0;

    /*! \brief Function to load the capture map. 
     
      Should be called by the derived classes
      Interaction::operator<<(const magnet::xml::Node&) function.
     */
    void loadCaptureMap(const magnet::xml::Node&);

    /*! \brief Function to write out the capture map. 
     
      Should be called by the derived classes Interaction::outputXML()
      function.
     */
    void outputCaptureMap(magnet::xml::XmlStream&) const;

    virtual void testAddToCaptureMap(const Particle& p1, const size_t& p2) const;

    //! \brief Add a pair of particles to the capture map
    void addToCaptureMap(const Particle& p1, const Particle& p2) const
    {
#ifdef DYNAMO_DEBUG
      if (captureMap.count(cMapKey(p1.getID(), p2.getID())))
	M_throw() << "Insert found " << std::min(p1.getID(), p2.getID())
		  << " and " << std::max(p1.getID(), p2.getID()) << " in the capture map";
#endif
    
      captureMap.insert(cMapKey(p1.getID(), p2.getID()));
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

  };

  /*! \brief This base class is for Interaction classes which "capture"
   * particle pairs in multiple states.
   */
  class IMultiCapture: public ICapture
  {
  public:
    IMultiCapture(dynamo::Simulation* sim, C2Range* range): ICapture(sim, range) {}

    size_t getTotalCaptureCount() const { return captureMap.size(); }
  
    virtual bool isCaptured(const Particle& p1, const Particle& p2) const
    { return captureMap.count(cMapKey(p1.getID(), p2.getID())); }

    virtual void clear() const { captureMap.clear(); }

  protected:
  
    typedef std::tr1::unordered_map<cMapKey, int, boost::hash<cMapKey> > captureMapType;
    typedef captureMapType::iterator cmap_it;
    typedef captureMapType::const_iterator const_cmap_it;

    mutable captureMapType captureMap;

    virtual int captureTest(const Particle&, const Particle&) const = 0;

    void loadCaptureMap(const magnet::xml::Node&);

    void outputCaptureMap(magnet::xml::XmlStream&) const;

    inline cmap_it getCMap_it(const Particle& p1, const Particle& p2) const
    { return captureMap.find(cMapKey(p1.getID(), p2.getID())); }

    inline void addToCaptureMap(const Particle& p1, const Particle& p2) const
    {
#ifdef DYNAMO_DEBUG
      if (captureMap.find(cMapKey(p1.getID(), p2.getID())) != captureMap.end())
	M_throw() << "Adding a particle while its already added!";
#endif
    
      captureMap[cMapKey(p1.getID(), p2.getID())] = 1;
    }

    //! \brief Add a pair of particles to the capture map
    void addToCaptureMap(const Particle& p1, const size_t& p2) const;

    virtual void testAddToCaptureMap(const Particle& p1, const size_t& p2) const;

    inline void delFromCaptureMap(const Particle& p1, const Particle& p2) const
    {
#ifdef DYNAMO_DEBUG
      if (captureMap.find(cMapKey(p1.getID(), p2.getID())) == captureMap.end())
	M_throw() << "Deleting a particle while its already gone!";
#endif 
      captureMap.erase(cMapKey(p1.getID(), p2.getID()));
    }
  };
}
