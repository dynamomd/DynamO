/*  DYNAMO:- Event driven molecular dynamics simulator 
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

#include "interaction.hpp"
#include "../../simulation/particle.hpp"
#include <boost/tr1/unordered_set.hpp>
#include <boost/tr1/unordered_map.hpp>
#include <set>

//Expose the boost TR1
namespace std {
  using namespace tr1;
}

//! This class is a general interface to Interaction classes that
//! allow particles to capture each other and store internal energy.
class ICapture
{
public:
  //! Returns the number of particles that are captured in some way
  virtual size_t getTotalCaptureCount() const = 0;
  
  //! A test if two particles are captured
  virtual bool isCaptured(const Particle&, const Particle&) const = 0;

  //! Returns the total internal energy stored in this Interaction.
  virtual double getInternalEnergy() const = 0;
  
protected:

};

//! This base class is for Interaction classes which only capture in
//! one state. 
//!
//! There is only one state a pair of particles can be in, either
//! captured or not.  This can be contrasted with IMultiCapture where
//! a pair of particles may be in a range of captured states.
class ISingleCapture: public ICapture, virtual public Interaction
{
public:
  ISingleCapture():noXmlLoad(true) {}

  //! \sa ICapture::getTotalCaptureCount()
  size_t getTotalCaptureCount() const { return captureMap.size(); }
  
  //! \sa ICapture::isCaptured()
  virtual bool isCaptured(const Particle& p1, const Particle& p2) const
  {
#ifdef DYNAMO_DEBUG
    if (p1.getID() == p2.getID())
      M_throw() << "Particle is testing if it captured itself";
#endif 
    
    return (p1.getID() < p2.getID())
      ? captureMap.count(std::pair<size_t, size_t>(p1.getID(), p2.getID()))
      : captureMap.count(std::pair<size_t, size_t>(p2.getID(), p1.getID()));
  }

protected:

  mutable std::unordered_set<std::pair<size_t, size_t> > captureMap;

  //! This function should provide a test of the particles current
  //! position and velocity to determine if they're captured. Used in
  //! rebuilding the captureMap.
  virtual bool captureTest(const Particle&, const Particle&) const = 0;

  bool noXmlLoad;

  void initCaptureMap();

  //! Function to load the capture map. Should be called by the
  //! derived classes Interaction::operator<<(const magnet::xml::Node&) function.
  void loadCaptureMap(const magnet::xml::Node&);

  //! Function to write out the capture map. Should be called by the
  //! derived classes Interaction::outputXML() function.
  void outputCaptureMap(xml::XmlStream&) const;

  //! Add a pair of particles to the capture map.
  void addToCaptureMap(const Particle&, const Particle&) const;
  
  //! Remove a pair of particles to the capture map.
  void removeFromCaptureMap(const Particle&, const Particle&) const;
  
};

//! This base class is for Interaction classes which only capture in
//! one state. 
//!
//! There is only one state a pair of particles can be in, either
//! captured or not.  This can be contrasted with IMultiCapture where
//! a pair of particles may be in a range of captured states.
//!
class IMultiCapture: public ICapture, virtual public Interaction
{
public:
  IMultiCapture():noXmlLoad(true) {}

  size_t getTotalCaptureCount() const { return captureMap.size(); }
  
  virtual bool isCaptured(const Particle&, const Particle&) const;

protected:

  struct cMapKey: public std::pair<size_t,size_t>
  {
    inline cMapKey(const size_t&a, const size_t&b):
      std::pair<size_t,size_t>(std::min(a,b), std::max(a,b))
    {
      //size_t t1(a),t2(b), t3(a);
      //asm ("PMINUW %1,%0" : "=i" (t1) : "i" (t2));
      //asm ("PMAXUW %1,%0" : "=i" (t2) : "i" (t3));
      //first = t1;
      //second = t2;
    }
  };
  
  struct keyHash
  {
    inline size_t operator()(const cMapKey& key) const
    {
      //return key.first * key.second + (key.first ^ key.second);
      //return (key.first > key.second) ?  ((key.first << 16) ^ key.second) : ((key.second << 16) ^ key.first);
      return (key.first << (sizeof(size_t) * 4)) ^ key.second;
    }
  };

  typedef std::unordered_map<cMapKey, int, keyHash> captureMapType;
  typedef captureMapType::iterator cmap_it;
  typedef captureMapType::const_iterator const_cmap_it;

  mutable captureMapType captureMap;

  virtual int captureTest(const Particle&, const Particle&) const = 0;

  bool noXmlLoad;

  void initCaptureMap();

  void loadCaptureMap(const magnet::xml::Node&);

  void outputCaptureMap(xml::XmlStream&) const;

  inline cmap_it getCMap_it(const Particle& p1, const Particle& p2) const
  {
    return captureMap.find(cMapKey(p1.getID(), p2.getID()));
  }

  inline void addToCaptureMap(const Particle& p1, const Particle& p2) const
  {
#ifdef DYNAMO_DEBUG
    if (captureMap.find(cMapKey(p1.getID(), p2.getID())) != captureMap.end())
      M_throw() << "Adding a particle while its already added!";
#endif
    
    captureMap[cMapKey(p1.getID(), p2.getID())] = 1;
  }

  inline void delFromCaptureMap(const Particle& p1, const Particle& p2) const
  {
#ifdef DYNAMO_DEBUG
    if (captureMap.find(cMapKey(p1.getID(), p2.getID())) == captureMap.end())
      M_throw() << "Deleting a particle while its already gone!";
#endif 
    captureMap.erase(cMapKey(p1.getID(), p2.getID()));
  }
};
