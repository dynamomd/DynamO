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

#include <cfloat>
#include "../eventtypes.hpp"
#include "../../simulation/particle.hpp"

class XMLNode;
namespace xml
{
  class XmlStream;
}
namespace DYNAMO {
  class SimData;
}
class IntEvent;
class Global;

class GlobalEvent
{
public:  
  GlobalEvent (const Particle&, const double&, 
	      EEventType, const Global&);

  inline bool operator== (const Particle &partx) const 
    { return (*particle_ == partx); }
  
  bool areInvolved(const IntEvent&) const; 
  
  inline void invalidate() 
  { 
    dt = DBL_MAX; 
    CType = NONE; 
  }

  inline bool operator< (const GlobalEvent & C2) const 
  { return dt < C2.dt;}
  
  inline bool operator> (const GlobalEvent & C2) const 
    { return dt > C2.dt;}

  inline void incrementTime(const double& deltat) {dt -= deltat; }

  inline void addTime(const double& deltat) {dt += deltat; }

  inline const Particle& getParticle() const { return *particle_; }

  inline const double& getdt() const { return dt; }
  inline void setdt(double nt) { dt = nt; }

  inline EEventType getType() const
    { return CType; }
  
  friend xml::XmlStream& operator<<(xml::XmlStream&, const GlobalEvent&);

  std::string stringData(const DYNAMO::SimData*) const;

  const size_t& getGlobalID() const { return globalID; } 

  inline void scaleTime(const double& scale)
  { dt *= scale; }

protected:
  const Particle*  particle_;
  double dt;
  EEventType CType;
  size_t globalID;
};
