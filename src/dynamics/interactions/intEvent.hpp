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

#ifndef CIntEvent_H
#define CIntEvent_H

#include <cfloat>
#include "../eventtypes.hpp"
#include "../../simulation/particle.hpp"


class XMLNode;
class CInteraction;
namespace DYNAMO {
  class SimData;
}
namespace xmlw
{
  class XmlStream;
}

class CIntEvent
{
public:  
  //A way to recover the collision name from a type at compile time
  static const char* getCollEnumName(EEventType);

  friend struct intPart;

  inline CIntEvent ():
  particle1(NULL), particle2(NULL), 
  dt(HUGE_VAL), CType(NONE),
  ptrInteraction(NULL) {}

  inline CIntEvent(const CParticle& part1, const CParticle& part2, 
		     const Iflt &delt, EEventType nType, 
		     const CInteraction& pI):
  particle1(&part1), particle2(&part2), dt(delt), 
  CType(nType), ptrInteraction(&pI) {}

  inline CIntEvent (const CParticle& part1):
  particle1(&part1), particle2(NULL), 
  dt(HUGE_VAL), CType(NONE),
  ptrInteraction(NULL) {}

  inline CIntEvent(const CParticle& part1, const Iflt& dt, 
		     EEventType etype):
  particle1(&part1), particle2(NULL), 
  dt(dt), CType(etype),
  ptrInteraction(NULL) {}

  inline bool operator== (const CParticle &partx) const 
    { return ((*particle1 == partx) || (*particle2 == partx)); }
  
  inline bool areInvolved(const CIntEvent &coll) const 
    { return ((coll == *particle1) || (coll == *particle2)); }
  
  inline void invalidate() 
  { 
    dt = DBL_MAX; 
    CType = NONE; 
  }

  inline bool operator< (const CIntEvent & C2) const 
  { return dt < C2.dt;}
  
  inline bool operator> (const CIntEvent & C2) const 
    { return dt > C2.dt;}

  inline void incrementTime(const Iflt deltat) {dt -= deltat; }

  inline void addTime(const Iflt deltat) {dt += deltat; }

  inline const CParticle& getParticle1() const { return *particle1; }

  inline const CParticle& getParticle2() const { return *particle2; }

  inline const CParticle* getPtrParticle2() const { return particle2; }

  inline bool hasParticle2() const { return particle2 != NULL; }

  inline const Iflt& getdt() const { return dt; }

  inline EEventType getType() const
    { return CType; }
  
  friend xmlw::XmlStream& operator<<(xmlw::XmlStream&, const CIntEvent&);

  std::string stringData(const DYNAMO::SimData*) const;

  inline void setType(EEventType a) const
  { 
    //A bit of nastiness required by SquareWells
    CType = a;
  }

  inline void scaleTime(const Iflt& scale)
  { dt *= scale; }

  const CInteraction& getInteraction() const 
  { 
    return *ptrInteraction; 
  }

private:
  const CParticle*  particle1;
  const CParticle*  particle2;
  Iflt dt;
  mutable EEventType CType;
  const CInteraction* ptrInteraction;
};

#endif
