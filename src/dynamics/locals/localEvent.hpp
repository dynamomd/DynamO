/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef CLocalEvent_H
#define CLocalEvent_H

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
class Local;

class LocalEvent
{
public:  
  LocalEvent(const Particle&, const double&,
	      EEventType, const Local&);

  inline bool operator== (const Particle &partx) const 
    { return (*particle_ == partx); }
  
  bool areInvolved(const IntEvent&) const; 
  
  inline void invalidate() 
  { 
    dt = DBL_MAX; 
    CType = NONE; 
  }

  inline bool operator< (const LocalEvent & C2) const 
  { return dt < C2.dt;}
  
  inline bool operator> (const LocalEvent & C2) const 
    { return dt > C2.dt;}

  inline void incrementTime(const double& deltat) {dt -= deltat; }

  inline void addTime(const double& deltat) {dt += deltat; }

  inline const Particle& getParticle() const { return *particle_; }

  inline const double& getdt() const { return dt; }

  inline EEventType getType() const
  { return CType; }
  
  friend xml::XmlStream& operator<<(xml::XmlStream&, const LocalEvent&);

  std::string stringData(const DYNAMO::SimData*) const;

  const size_t& getLocalID() const { return localID; } 

  inline void scaleTime(const double& scale)
  { dt *= scale; }

protected:
  const Particle*  particle_;
  double dt;
  mutable EEventType CType;
  const size_t localID;
};

#endif
