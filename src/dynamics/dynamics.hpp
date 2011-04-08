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

#include <magnet/cloneptr.hpp>
#include "../base/is_base.hpp"
#include <vector>
#include "interactions/interaction.hpp"
#include "interactions/intEvent.hpp"
#include <boost/foreach.hpp>

class BoundaryCondition;
class Units;
class Species;
class GlobalEvent;
class Global;
class Local;
class LocalEvent;
class System;
class Topology;
class Particle;
template<class T> class CVector;
class Liouvillean;
class NEventData;
class PairEventData;
class ParticleEventData;
namespace magnet {
  namespace xml
  { class Node; }
}

class Dynamics: public DYNAMO::SimBase
{
public:
  //Constructors
  Dynamics(DYNAMO::SimData*);

  Dynamics(const magnet::xml::Node& XML, DYNAMO::SimData*);

  explicit Dynamics(const Dynamics &);

  ~Dynamics();
  
  void setUnits(Units*);

  void setLiouvillean(Liouvillean*);

  Interaction* addInteraction(Interaction*);

  void addSpecies(const magnet::ClonePtr<Species>&);
  
  void addGlobal(Global*);

  void addLocal(Local*);

  void addSystem(System*);

  void addStructure(Topology*);

  const Species& getSpecies(const Particle&) const;
  
  const magnet::ClonePtr<Interaction>& 
    getInteraction(const Particle&, const Particle&) const; 
  
  void stream(const double&);
  
  inline IntEvent getEvent(const Particle& p1, const Particle& p2) const
  {
    BOOST_FOREACH(const magnet::ClonePtr<Interaction>& ptr, interactions)
      if (ptr->isInteraction(p1,p2))
	{
#ifdef DYNAMO_UpdateCollDebug
	  std::cerr << "\nGOT INTERACTION P1 = " << p1.getID() << " P2 = " 
		    << p2.getID() << " NAME = " << typeid(*ptr.get_ptr()).name();
#endif
	  return ptr->getEvent(p1,p2);
	}
    
    M_throw() << "Could not find the right interaction to test for";
  }

  void operator<<(const magnet::xml::Node&);

  void initialise();
  
  double getLongestInteraction() const;
  
  /*! \brief Sets the Centre of Mass (COM) velocity of the system 
   * 
   *  The COM momentum of the system is
   * \f[ \bm{P}_{system} = \sum_i m_i \bm{v}_i \f]
   * 
   * We want to first remove any motion of the system, so we subtract
   * the COM momentum based on the mass of each particle (E.g. \f$ m_i
   * / \sum_j m_j\f$). This has two nice effects, first, particles
   * store their velocities, not their momentums so we convert by
   * dividing by \f$m_i\f$ which gives 
   *
   * \f[ \bm{v}_i \to \bm{v}_i -
   * (\sum_i m_i \bm{v}_i) / \sum_i m_i \f] 
   *
   * So relative velocities are preserved as the subtraction is a
   * constant for all particles. Also we can now just add the offset to give
   *
   * \f[ \bm{v}_i \to \bm{v}_i -(\sum_i m_i \bm{v}_i) / \sum_i m_i  + \bm{V}_{COM}\f]  
   *
   * \param COMVelocity The target velocity for the COM of the system.
   */  
  void setCOMVelocity(const Vector COMVelocity = Vector(0,0,0));

  void rescaleLengths(double);

  void SystemOverlapTest();
  
  double calcInternalEnergy() const;

  Dynamics* Clone() const { return new Dynamics(*this); }

  std::vector<magnet::ClonePtr<Interaction> >& getInteractions() { return interactions; }
  const std::vector<magnet::ClonePtr<Interaction> >& getInteractions() const { return interactions; }

  magnet::ClonePtr<Interaction>& getInteraction(std::string);
  const magnet::ClonePtr<Interaction>& getInteraction(std::string) const;

  const std::vector<magnet::ClonePtr<Global> >& getGlobals() const { return globals; }
  std::vector<magnet::ClonePtr<Global> >& getGlobals() { return globals; }
  magnet::ClonePtr<Global>& getGlobal(std::string);
  const magnet::ClonePtr<Global>& getGlobal(std::string) const;

  std::vector<magnet::ClonePtr<Local> >& getLocals() { return locals; }
  const std::vector<magnet::ClonePtr<Local> >& getLocals() const { return locals; }
  magnet::ClonePtr<Local>& getLocal(std::string);
  const magnet::ClonePtr<Local>& getLocal(std::string) const;

  const std::vector<magnet::ClonePtr<Species> >& getSpecies() const { return species; }
  const Species& getSpecies(std::string) const;
  Species& getSpecies(std::string);

  std::vector<magnet::ClonePtr<Topology> >& getTopology() { return topology; }
  const std::vector<magnet::ClonePtr<Topology> >& getTopology() const { return topology; }

  magnet::ClonePtr<Topology>& getTopology(std::string);

  const magnet::ClonePtr<Topology>& getTopology(std::string) const;

  std::vector<magnet::ClonePtr<System> >& getSystemEvents() { return systems; }
  const std::vector<magnet::ClonePtr<System> >& getSystemEvents() const { return systems; }
  const magnet::ClonePtr<System>& getSystem(std::string) const;
  magnet::ClonePtr<System>& getSystem(std::string);

  void addSystemTicker();
  
  friend xml::XmlStream& operator<<(xml::XmlStream&, const Dynamics&);

  //Inlines
  inline const Units& units() const 
  { return *p_units; }

  inline Units& units()
  { return *p_units; }
  
  inline const BoundaryCondition& BCs() const 
  { return *p_BC; }

  inline const Liouvillean& getLiouvillean() const
  { return *p_liouvillean; }

  inline  Liouvillean& getLiouvillean()
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
  template<class T> void applyBC()
    {
      if (p_BC.empty())
	I_cout() << "Warning, resetting the BC's";
      
      p_BC.set_ptr(new T(Sim));
    }

  double getNumberDensity() const;
  
  double getPackingFraction() const;

 protected:
  void outputXML(xml::XmlStream &) const;

  std::vector<magnet::ClonePtr<Interaction> > interactions;
  std::vector<magnet::ClonePtr<Global> > globals;
  std::vector<magnet::ClonePtr<Local> > locals;
  std::vector<magnet::ClonePtr<System> > systems;
  std::vector<magnet::ClonePtr<Topology> > topology;
  std::vector<magnet::ClonePtr<Species> > species;
  magnet::ClonePtr<BoundaryCondition> p_BC;
  magnet::ClonePtr<Liouvillean> p_liouvillean;
  mutable magnet::ClonePtr<Units> p_units;
};
