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
#include <dynamo/base.hpp>
#include <dynamo/eventtypes.hpp>

namespace magnet {
namespace xml {
class Node;
class XmlStream;
} // namespace xml
} // namespace magnet

namespace dynamo {
/*! \brief Future Event Lists (FEL) sort the Particle Event Lists
    (PEL) to determine the next event to occur.

    Classes Derived from this base class provide a mechanism to sort
    \ref Event s. These events are first pre-sorted using a Particle
    Event List before being sorted by these classes.
 */
class FEL {
public:
  virtual ~FEL() {}

  virtual void clear() = 0;

  // An implementation of size() is impossible, as events may be
  // discarded by the PELs and replaced by RECALCULATE events.
  // virtual size_t size() const = 0;

  /*! \brief Test if the queue is out of events.

    This function requires the FEL to be in a sorted state. (all
    pushes/pops/invalidates have been followed by updated_particle()
    or updated_all() calls).

    \param ID The particle ID, whose events have become invalid.
   */
  virtual bool empty() = 0;

  /*! \brief Initialise the FEL and prepare it for a maximum of N
      particles.

    \param N The maximum particle ID that will be pushed to the
    queue.
   */
  virtual void init(const size_t N) = 0;

  /*! \brief Invalidate all events involving a particle.

    This function may place the FEL in a "unsorted" state until
    updated_particle() is called with the same particle ID (or
    updated_all() is called). This is to allow calls to invalidate()
    and push() to be combined for efficiency.

    \param ID The particle ID, whose events have become invalid.
   */
  virtual void invalidate(const size_t ID) = 0;

  /*! \brief Remove the next event in the queue.

    This function may place the FEL in a "unsorted" state until
    updated_particle() is called with the same particle ID as the
    event removed (or updated_all() is called). This is to allow
    calls to invalidate(), pop(), and push() to be combined for
    efficiency.
   */
  virtual void pop() = 0;

  /*! \brief Add an event to the FEL.

     This function may place the FEL in a "unsorted" state until
     updated_particle() is called with the same particle ID as the
     event's _particle1ID (or updated_all() is called). This is to
     allow calls to invalidate() and push() to be combined for
     efficiency.

     \param event The new event to push.
    */
  virtual void push(Event event) = 0;

  virtual void rescaleTimes(const double) = 0;
  virtual void stream(const double) = 0;

  virtual Event top() = 0;

  static shared_ptr<FEL> getClass(const magnet::xml::Node &);
  friend ::magnet::xml::XmlStream &operator<<(::magnet::xml::XmlStream &,
                                              const FEL &);

private:
  virtual void outputXML(magnet::xml::XmlStream &) const = 0;
};
} // namespace dynamo
