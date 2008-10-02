/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef DYNAMICS_H
#define DYNAMICS_H

#include "../datatypes/pluginpointer.hpp"
#include "../base/is_base.hpp"
#include <vector>
#include "interactions/interaction.hpp"
#include "interactions/intEvent.hpp"
#include <boost/foreach.hpp>

class CBC;
class CUnits;
class CSpecies;
class CGlobEvent;
class CGlobal;
class CSystem;
class CTopology;
class CParticle;
template<class T>
class CVector;
class CLiouvillean;
class XMLNode;
class CNParticleData;
class C2ParticleData;
class C1ParticleData;
namespace xmlw
{
  class XmlStream;
}

class CDynamics: public DYNAMO::SimBase
{
public:
  //Constructors
  CDynamics(DYNAMO::SimData*);

  CDynamics(const XMLNode&, DYNAMO::SimData*);

  explicit CDynamics(const CDynamics &);

  ~CDynamics();
  
  void setUnits(CUnits*);

  void setLiouvillean(CLiouvillean*);

  CInteraction* addInteraction(CInteraction*);

  void addSpecies(CSpecies);
  
  void addGlobal(CGlobal*);

  void addSystem(CSystem*);

  void addStructure(CTopology*);

  const CSpecies& getSpecies(const CParticle&) const;
  
  void setAspectRatio(const CVector<> &);

  const smrtPlugPtr<CInteraction>& 
    getInteraction(const CParticle&, const CParticle&) const; 
  
  void stream(const Iflt&);
  
  C2ParticleData runEvent(const CIntEvent&);
  
  CNParticleData runEvent(const CGlobEvent&);

  inline CIntEvent getEvent(const CParticle &p1, const CParticle &p2) const
  {
    BOOST_FOREACH(const smrtPlugPtr<CInteraction>& ptr, interactions)
      if (ptr->isInteraction(p1,p2))
#ifdef DYNAMO_UpdateCollDebug
	{
	  std::cerr << "\nGOT INTERACTION P1 = " << p1.getID() << " P2 = " 
		    << p2.getID() << " NAME = " << typeid(*ptr.get_ptr()).name();
	  return ptr->getCollision(p1,p2);
	}
# else
        return ptr->getCollision(p1,p2);
#endif
    
    I_throw() << "Could not find the right interaction to test for";
  }

  CGlobEvent getEvent(const CParticle&) const;
  
  void operator<<(const XMLNode&);


  void initialise();
  
  Iflt getLongestInteraction() const;
  
  void zeroMomentum(std::vector<CParticle> &);

  void rescaleLengths(Iflt);

  void SystemOverlapTest();
  
  CVector<> getLabVelocity(const CParticle &) const;
  
  Iflt calcInternalEnergy() const;

  CDynamics* Clone() const { return new CDynamics(*this); }

  std::vector<smrtPlugPtr<CInteraction> >& getInteractions();
  const std::vector<smrtPlugPtr<CInteraction> >& getInteractions() const;
  smrtPlugPtr<CInteraction>& getInteraction(std::string);
  const smrtPlugPtr<CInteraction>& getInteraction(std::string) const;

  const std::vector<smrtPlugPtr<CGlobal> >& getGlobals() const;
  smrtPlugPtr<CGlobal>& getGlobal(std::string);
  const smrtPlugPtr<CGlobal>& getGlobal(std::string) const;

  const std::vector<CSpecies>& getSpecies() const;
  const CSpecies& getSpecies(std::string) const;

  std::vector<smrtPlugPtr<CTopology> >& getTopology();
  const std::vector<smrtPlugPtr<CTopology> >& getTopology() const;
  smrtPlugPtr<CTopology>& getTopology(std::string);
  const smrtPlugPtr<CTopology>& getTopology(std::string) const;

  std::vector<smrtPlugPtr<CSystem> >& getSystemEvents();
  const std::vector<smrtPlugPtr<CSystem> >& getSystemEvents() const;
  const smrtPlugPtr<CSystem>& getSystem(std::string) const;
  smrtPlugPtr<CSystem>& getSystem(std::string);

  void addSystemTicker();
  
  friend xmlw::XmlStream& operator<<(xmlw::XmlStream&, const CDynamics&);

  //Inlines
  inline const CUnits& units() const 
  { return *p_units; }

  inline CUnits& units()
  { return *p_units; }
  
  inline const CBC& BCs() const 
  { return *p_BC; }

  inline const CLiouvillean& Liouvillean() const
  { return *p_liouvillean; }

  inline  CLiouvillean& Liouvillean()
  { return *p_liouvillean; }

  template<class T>
  inline bool liouvilleanTypeTest() const
  { return dynamic_cast<const T*>(p_liouvillean.get_ptr()) != NULL; }

  template<class T>
  inline bool BCTypeTest() const
  { return dynamic_cast<const T*>(p_BC.get_ptr()) != NULL; }

  template<class T>
  inline bool unitTypeTest() const
  { return dynamic_cast<const T*>(p_units.get_ptr()) != NULL; }

  //templates
  template<class T> void setPBC()
    {
      if (p_BC != NULL)
	I_cout() << "Warning, resetting the BC's";
      
      p_BC.set_ptr(new T(Sim));
    }

  Iflt getKineticEnergy() const;
  
  CVector<> getVecKineticEnergy() const;

  Iflt getkT() const;
  
  CVector<> getVeckT() const;

  Iflt getNumberDensity() const;
  
  Iflt getPackingFraction() const;

  Iflt getParticleEnergy(const CParticle&) const;

 protected:
  void outputXML(xmlw::XmlStream &) const;

  std::vector<smrtPlugPtr<CInteraction> > interactions;
  std::vector<smrtPlugPtr<CGlobal> > globals;
  std::vector<smrtPlugPtr<CSystem> > systems;
  std::vector<smrtPlugPtr<CTopology> > topology;
  std::vector<CSpecies> species;
  smrtPlugPtr<CBC> p_BC;
  smrtPlugPtr<CLiouvillean> p_liouvillean;
  mutable smrtPlugPtr<CUnits> p_units;
};

#endif
