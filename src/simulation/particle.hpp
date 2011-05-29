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

#include <magnet/math/vector.hpp>

namespace magnet { namespace xml { class Node; } }
namespace xml { class XmlStream; }

//! \brief The fundamental data structure for a Particle.
//!
//! This class holds only the very fundamental information on a
//! particle, such as its position, velocity, ID, and state
//! flags. Other data is "attached" to this particle using
//! Property classes stored in the PropertyStore.
class Particle
{
public:
  //! \brief Operator to write out an XML representation of a Particle.
  friend xml::XmlStream& operator<<(xml::XmlStream&, const Particle&);
  
  //! \brief Constructor to build a particle from passed values.
  inline Particle (const Vector  &position, 
		   const Vector  &velocity,
		   const unsigned long& nID):
    _pos(position), _vel(velocity), 
    _ID(nID), _peculiarTime(0.0),
    _state(DEFAULT)
  {}
  
  //! \brief Constructor to build a particle from an XML node.
  Particle(const magnet::xml::Node&, unsigned long);

  //! \brief Equal to comparison operator.
  //! This comparison operator only compares the ID's of the Particle
  //! classes. 
  inline bool operator==(const Particle &p) const { return (_ID == p._ID); }
  //! \brief Not equal to comparison operator.
  //! This comparison operator only compares the ID's of the Particle
  //! classes. 
  inline bool operator!=(const Particle &p) const { return (_ID != p._ID); }
  
  //! \brief Const position accessor function.
  inline const Vector  &getPosition() const { return _pos; }
  //! \brief Const velocity accessor function.
  inline const Vector  &getVelocity() const { return _vel; }
  
  //! \brief Position accessor function.
  inline Vector& getPosition() { return _pos; }
  //! \brief Velocity accessor function.
  inline Vector& getVelocity() { return _vel; }
  
  //! \brief ID accessor function.
  //! This ID is a unique value for each Particle in the Simulation
  //! and so it can also be used as a reference to a particle.
  inline const unsigned long &getID() const { return _ID; };

  //! \brief Const peculiar time accessor function.
  //! This value is used in the "delayed states" or "Time warp" algorithm.
  inline const double& getPecTime() const { return _peculiarTime; }
  //! \brief Peculiar time accessor function.
  //! This value is used in the "delayed states" or "Time warp" algorithm.
  inline double& getPecTime() { return _peculiarTime; }
  
  //! \brief The possible State flags of the Particle, these states may be combined.
  typedef enum {
    DEFAULT = 0x01 | 0x02,//!< The default flags for the Particle's State.
    DYNAMIC = 0x01, //!< For the LNewtonianGravity Liouvillean it Enables/Disables the gravity force for acting on this Particle.
    ALIVE = 0x02 //!< Flags if the particle is actually in the Simulation.
  } State;
  
  //! \brief Used to test if the Particle has a State flag set.
  //! \param teststate The State flag to test.
  inline bool testState(State teststate) const { return _state & teststate; }

  //! \brief Sets a State flag of the Particle
  //! \param nState The State flag to set.
  inline void setState(State nState) { _state |= nState; }

  //! \brief Clears a State flag of the Particle
  //! \param nState The State flag to clear.
  inline void clearState(State nState) { _state &= (~nState); }  

private:
  Vector _pos;
  Vector _vel;
  unsigned long _ID;
  double _peculiarTime;
  int _state;
};
