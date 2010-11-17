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

#pragma once

#include "../datatypes/vector.hpp"

class XMLNode;
namespace xml
{
  class XmlStream;
}

struct ParticleState
{
};

class Particle
{
public:
  friend xml::XmlStream& operator<<(xml::XmlStream&, const Particle&);
  
  inline Particle (const Vector  &position, 
		    const Vector  &velocity,
		    const unsigned long& nID):
    _pos(position), _vel(velocity), 
    _ID(nID), _peculiarTime(0.0),
    _state(DEFAULT)
  {}
  
  Particle(const XMLNode&, unsigned long);
  
  inline bool operator==(const Particle &p) const { return (_ID == p._ID); }
  inline bool operator!=(const Particle &p) const { return (_ID != p._ID); }
  
  inline const Vector  &getPosition() const { return _pos; }
  inline const Vector  &getVelocity() const { return _vel; }
  
  inline Vector  &getPosition() { return _pos; }
  inline Vector  &getVelocity() { return _vel; }
  
  inline const unsigned long &getID() const { return _ID; };
  inline const double& getPecTime() const { return _peculiarTime; }
  inline double& getPecTime() { return _peculiarTime; }
  
  inline void scaleVelocity(const double& vs) { _vel *= vs; }
  inline void scalePosition(const double& vs) { _pos *= vs; }
  

  typedef enum {
    DEFAULT = 0x01 | 0x02,
    DYNAMIC = 0x01, //Will the free streaming part of the
    //Liouvillean be applied to his particle
    //and will events be generated for this
    //particle.                            
    ALIVE = 0x02 //Is this particle in the simulation
  } State;

  inline bool testState(State teststate) const { return _state & teststate; }
  inline void setState(State nState) { _state |= nState; }
  inline void clearState(State nState) { _state &= (~nState); }
  
private:
  //NOTE, changing these members type must be reflected in the liouvillean
  //where binary data is written
  Vector _pos;
  Vector _vel;
  unsigned long _ID;
  double _peculiarTime;

  int _state;
};
