/*  dynamo:- Event driven molecular dynamics simulator 
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

#include <dynamo/dynamics/interactions/interaction.hpp>
#include <dynamo/dynamics/interactions/intEvent.hpp>
#include <dynamo/dynamics/units/units.hpp>
#include <tr1/memory>
#include <boost/foreach.hpp>
#include <vector>
#include <tr1/memory>

namespace magnet {    
  namespace xml { 
    class Node; 
  }
}

namespace dynamo {
  class BoundaryCondition;
  class Particle;
  class Species;
  class GlobalEvent;
  class Global;
  class Local;
  class LocalEvent;
  class System;
  class Topology;
  class Particle;
  class Liouvillean;
  class NEventData;
  class PairEventData;
  class ParticleEventData;

  class Dynamics: public dynamo::SimBase
  {
  public:
    //Constructors
    Dynamics(dynamo::SimData*);

    Dynamics(const magnet::xml::Node& XML, dynamo::SimData*);

    explicit Dynamics(const Dynamics &);

    ~Dynamics();
  
    inline void setLiouvillean(std::tr1::shared_ptr<Liouvillean> ptr)
    { p_liouvillean = ptr; }

    Interaction* addInteraction(Interaction*);

    void addSpecies(const std::tr1::shared_ptr<Species>&);
  
    void addGlobal(Global*);

    void addLocal(Local*);

    void addSystem(System*);

    void addStructure(Topology*);

    const Species& getSpecies(const Particle&) const;
  
    const std::tr1::shared_ptr<Interaction>& 
    getInteraction(const Particle&, const Particle&) const; 
  
    void stream(const double&);
  
    inline IntEvent getEvent(const Particle& p1, const Particle& p2) const
    {
      BOOST_FOREACH(const std::tr1::shared_ptr<Interaction>& ptr, interactions)
	if (ptr->isInteraction(p1,p2))
	  {
#ifdef dynamo_UpdateCollDebug
	    std::cerr << "\nGOT INTERACTION P1 = " << p1.getID() << " P2 = " 
		      << p2.getID() << " NAME = " << typeid(*(ptr.get())).name();
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

    void SystemOverlapTest();
  
    double calcInternalEnergy() const;

    std::vector<std::tr1::shared_ptr<Interaction> >& getInteractions() { return interactions; }
    const std::vector<std::tr1::shared_ptr<Interaction> >& getInteractions() const { return interactions; }

    std::tr1::shared_ptr<Interaction>& getInteraction(std::string);
    const std::tr1::shared_ptr<Interaction>& getInteraction(std::string) const;

    const std::vector<std::tr1::shared_ptr<Global> >& getGlobals() const { return globals; }
    std::vector<std::tr1::shared_ptr<Global> >& getGlobals() { return globals; }
    std::tr1::shared_ptr<Global>& getGlobal(std::string);
    const std::tr1::shared_ptr<Global>& getGlobal(std::string) const;

    std::vector<std::tr1::shared_ptr<Local> >& getLocals() { return locals; }
    const std::vector<std::tr1::shared_ptr<Local> >& getLocals() const { return locals; }
    std::tr1::shared_ptr<Local>& getLocal(std::string);
    const std::tr1::shared_ptr<Local>& getLocal(std::string) const;

    const std::vector<std::tr1::shared_ptr<Species> >& getSpecies() const { return species; }
    const Species& getSpecies(std::string) const;
    Species& getSpecies(std::string);

    std::vector<std::tr1::shared_ptr<Topology> >& getTopology() { return topology; }
    const std::vector<std::tr1::shared_ptr<Topology> >& getTopology() const { return topology; }

    std::tr1::shared_ptr<Topology>& getTopology(std::string);

    const std::tr1::shared_ptr<Topology>& getTopology(std::string) const;

    std::vector<std::tr1::shared_ptr<System> >& getSystemEvents() { return systems; }
    const std::vector<std::tr1::shared_ptr<System> >& getSystemEvents() const { return systems; }
    const std::tr1::shared_ptr<System>& getSystem(std::string) const;
    std::tr1::shared_ptr<System>& getSystem(std::string);

    void addSystemTicker();
  
    friend magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream&, const Dynamics&);

    inline const Units& units() const { return _units; }

    inline Units& units() { return _units; }
  
    inline const BoundaryCondition& BCs() const 
    { return *p_BC; }

    inline const Liouvillean& getLiouvillean() const
    { return *p_liouvillean; }

    inline  Liouvillean& getLiouvillean()
    { return *p_liouvillean; }

    inline  std::tr1::shared_ptr<Liouvillean>& getLiouvilleanPtr()
    { return p_liouvillean; }

    template<class T>
    inline bool liouvilleanTypeTest() const
    { return std::tr1::dynamic_pointer_cast<T>(p_liouvillean); }

    template<class T>
    inline bool BCTypeTest() const
    { return std::tr1::dynamic_pointer_cast<T>(p_BC); }

    //templates
    template<class T> void applyBC()
    {
      if (p_BC)
	dout << "Warning, resetting the BC's" << std::endl;
      
      p_BC = std::tr1::shared_ptr<BoundaryCondition>(new T(Sim));
    }

    double getSimVolume() const;

    double getNumberDensity() const;
  
    double getPackingFraction() const;

  protected:
    void outputXML(magnet::xml::XmlStream &) const;

    std::vector<std::tr1::shared_ptr<Interaction> > interactions;
    std::vector<std::tr1::shared_ptr<Global> > globals;
    std::vector<std::tr1::shared_ptr<Local> > locals;
    std::vector<std::tr1::shared_ptr<System> > systems;
    std::vector<std::tr1::shared_ptr<Topology> > topology;
    std::vector<std::tr1::shared_ptr<Species> > species;
    std::tr1::shared_ptr<BoundaryCondition> p_BC;
    std::tr1::shared_ptr<Liouvillean> p_liouvillean;
    Units _units;
  };
}
