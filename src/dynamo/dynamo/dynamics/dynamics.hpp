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

#include <dynamo/dynamics/interactions/intEvent.hpp>
#include <dynamo/dynamics/units/units.hpp>
#include <boost/foreach.hpp>
#include <vector>

namespace magnet {    
  namespace xml { 
    class Node; 
  }
}

namespace dynamo {
  class Particle;
  class BoundaryCondition;
  class GlobalEvent;
  class Global;
  class Local;
  class LocalEvent;
  class System;
  class Particle;
  class NEventData;
  class PairEventData;
  class ParticleEventData;

  class Dynamics: public dynamo::SimBase
  {
  public:
    //Constructors
    Dynamics(dynamo::SimData*);

    Dynamics(const magnet::xml::Node& XML, dynamo::SimData*);

    ~Dynamics();
  
    void addGlobal(shared_ptr<Global> ptr);
    void addLocal(shared_ptr<Local> ptr);
    void addSystem(shared_ptr<System> ptr);
    
    void stream(const double&);
  
    void operator<<(const magnet::xml::Node&);

    void initialise();
    
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

    shared_ptr<Interaction>& getInteraction(std::string);
    const shared_ptr<Interaction>& getInteraction(std::string) const;

    const std::vector<shared_ptr<Global> >& getGlobals() const { return globals; }
    std::vector<shared_ptr<Global> >& getGlobals() { return globals; }
    shared_ptr<Global>& getGlobal(std::string);
    const shared_ptr<Global>& getGlobal(std::string) const;

    std::vector<shared_ptr<Local> >& getLocals() { return locals; }
    const std::vector<shared_ptr<Local> >& getLocals() const { return locals; }
    shared_ptr<Local>& getLocal(std::string);
    const shared_ptr<Local>& getLocal(std::string) const;

    std::vector<shared_ptr<System> >& getSystemEvents() { return systems; }
    const std::vector<shared_ptr<System> >& getSystemEvents() const { return systems; }
    const shared_ptr<System>& getSystem(std::string) const;
    shared_ptr<System>& getSystem(std::string);

    void addSystemTicker();
  
    friend magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream&, const Dynamics&);

    inline const Units& units() const { return _units; }

    inline Units& units() { return _units; }
  
    double getSimVolume() const;

    double getNumberDensity() const;
  
    double getPackingFraction() const;

  protected:
    Dynamics(const Dynamics &dyn);

    void outputXML(magnet::xml::XmlStream &) const;

    std::vector<shared_ptr<Global> > globals;
    std::vector<shared_ptr<Local> > locals;
    std::vector<shared_ptr<System> > systems;
    Units _units;
  };
}
