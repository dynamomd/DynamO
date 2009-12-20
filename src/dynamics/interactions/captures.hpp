/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef CICapture_H
#define CICapture_H

#include "interaction.hpp"
#include <set>
#include <vector>
#include <boost/tr1/unordered_set.hpp>
#include <boost/tr1/unordered_map.hpp>
#include "../../simulation/particle.hpp"

//Expose the boost TR1
namespace std {
  using namespace tr1;
}

class CICapture: public CInteraction
{
public:
  CICapture(DYNAMO::SimData*, C2Range*);

  virtual size_t getTotalCaptureCount() const = 0;
  
  virtual bool isCaptured(const CParticle&, const CParticle&) const = 0;

  virtual Iflt getInternalEnergy() const = 0;
  
protected:

};


class CISingleCapture: public CICapture
{
public:
  CISingleCapture(DYNAMO::SimData* Sim, C2Range* r):
    CICapture(Sim,r),
    noXmlLoad(true)
  {}

  size_t getTotalCaptureCount() const { return captureMap.size(); }
  
  virtual bool isCaptured(const CParticle& p1, const CParticle& p2) const
  {
#ifdef DYNAMO_DEBUG
    if (p1.getID() == p2.getID())
      D_throw() << "Particle is testing if it captured itself";
#endif 
    
    return (p1.getID() < p2.getID())
      ? captureMap.count(std::pair<size_t, size_t>(p1.getID(), p2.getID()))
      : captureMap.count(std::pair<size_t, size_t>(p2.getID(), p1.getID()));
  }

protected:

  mutable std::unordered_set<std::pair<size_t, size_t> > captureMap;

  virtual bool captureTest(const CParticle&, const CParticle&) const = 0;

  bool noXmlLoad;

  void initCaptureMap();

  void loadCaptureMap(const XMLNode&);

  void outputCaptureMap(xmlw::XmlStream&) const;

  void addToCaptureMap(const CParticle&, const CParticle&) const;
  
  void removeFromCaptureMap(const CParticle&, const CParticle&) const;
  
};

class CIMultiCapture: public CICapture
{
public:
  CIMultiCapture(DYNAMO::SimData* Sim, C2Range* r):
    CICapture(Sim,r),
    noXmlLoad(true)
  {}

  size_t getTotalCaptureCount() const { return captureMap.size(); }
  
  virtual bool isCaptured(const CParticle&, const CParticle&) const;

protected:
  
  //typedef std::pair<size_t, size_t> cMapKey;

  
  struct cMapKey: public std::pair<size_t,size_t>
  {
    inline cMapKey(const size_t&a, const size_t&b):
      std::pair<size_t,size_t>(a,b) {}

    inline bool operator==(const cMapKey& o) const
    {
      return (((first == o.first) && (second == o.second))
	      || ((second == o.first) && (first == o.second)));
    }
  };
  
  //  typedef boost::hash<cMapKey> keyHash;
  
  struct keyHash
  {
    inline size_t operator()(const cMapKey& key) const
    {
      return key.first * key.second + (key.first ^ key.second);
      //return (key.first > key.second) ?  ((key.first << 16) ^ key.second) : ((key.second << 16) ^ key.first);
    }
  };
  
  typedef std::unordered_map<cMapKey, int, keyHash> captureMapType;
  typedef captureMapType::iterator cmap_it;
  typedef captureMapType::const_iterator const_cmap_it;

  mutable captureMapType captureMap;

  virtual int captureTest(const CParticle&, const CParticle&) const = 0;

  bool noXmlLoad;

  void initCaptureMap();

  void loadCaptureMap(const XMLNode&);

  void outputCaptureMap(xmlw::XmlStream&) const;

  inline cmap_it getCMap_it(const CParticle& p1, const CParticle& p2) const
  {
    return captureMap.find(cMapKey(p1.getID(), p2.getID()));
  }

  inline void addToCaptureMap(const CParticle& p1, const CParticle& p2) const
  {
#ifdef DYNAMO_DEBUG
    if (captureMap.find(cMapKey(p1.getID(), p2.getID())) != captureMap.end())
      D_throw() << "Adding a particle while its already added!";
#endif
    
    captureMap[cMapKey(p1.getID(), p2.getID())] = 1;
  }

  inline void delFromCaptureMap(const CParticle& p1, const CParticle& p2) const
  {
#ifdef DYNAMO_DEBUG
    if (captureMap.find(cMapKey(p1.getID(), p2.getID())) == captureMap.end())
      D_throw() << "Deleting a particle while its already gone!";
#endif 
    captureMap.erase(cMapKey(p1.getID(), p2.getID()));
  }
};

#endif
