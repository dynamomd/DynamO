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

#ifndef PARTICLE_H
#define PARTICLE_H

#include "../datatypes/vector.hpp"

class XMLNode;
namespace xmlw
{
  class XmlStream;
}

class Particle
{
public:
  friend xmlw::XmlStream& operator<<(xmlw::XmlStream&, const Particle&);
  
  inline Particle (const Vector  &position, 
		    const Vector  &velocity,
		    const unsigned long& nID):
    posVector(position), velVector(velocity), 
    ID(nID), pecTime(0.0)
  {}
  
  Particle(const XMLNode&, unsigned long);
  
  inline bool operator==(const Particle &p) const { return (ID == p.ID); }
  inline bool operator!=(const Particle &p) const { return (ID != p.ID); }
  
  inline const Vector  &getPosition() const { return posVector; }
  inline const Vector  &getVelocity() const { return velVector; }
  
  inline Vector  &getPosition() { return posVector; }
  inline Vector  &getVelocity() { return velVector; }
  
  inline const unsigned long &getID() const { return ID; };
  inline const Iflt& getPecTime() const { return pecTime; }
  inline Iflt& getPecTime() { return pecTime; }
  
  inline void scaleVelocity(const Iflt& vs) { velVector *= vs; }
  inline void scalePosition(const Iflt& vs) { posVector *= vs; }
  
private:
  //NOTE, changing these members type must be reflected in the liouvillean
  //where binary data is written
  Vector posVector;
  Vector velVector;
  unsigned long ID;
  Iflt pecTime;
};

#endif
