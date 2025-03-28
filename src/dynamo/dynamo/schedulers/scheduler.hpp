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
#include <dynamo/ranges/IDRange.hpp>
#include <dynamo/schedulers/sorters/FEL.hpp>
#include <magnet/function/delegate.hpp>
#include <magnet/math/vector.hpp>
#include <memory>
#include <vector>

namespace magnet {
namespace xml {
class Node;
}
} // namespace magnet

namespace dynamo {
class Particle;
class Event;

class Scheduler : public dynamo::SimBase {
public:
  Scheduler(dynamo::Simulation *const, const char *, FEL *);

  virtual ~Scheduler() = 0;

  virtual void initialise();
  virtual void initialiseNBlist() = 0;

  void rebuildList();

  /*! \brief Retest for events for a single particle.
   */
  inline void fullUpdate(Particle &part) {
    invalidateEvents(part);
    addEvents(part);
  }

  /*! \brief Retest for events for two particles.

    Even though we would have less invalid events in the queue if we
    interleaved the following updates, we only want one valid event
    for the (p1,p2) interaction. So we still split the p1 and p2
    interactions.

    We want only one valid p1,p2 interaction to help prevent loops in
    the event recalculation code. So if we try to exectue one p1,p2
    interaction, but find the p2,p1 interaction is sooner by a
    numerically insignificant amount caused by being pushed into the
    sorter, we will enter a loop which has to be broken by the
    _interactionRejectionCounter logic.
  */
  inline void fullUpdate(Particle &p1, Particle &p2) {
    fullUpdate(p1);
    fullUpdate(p2);
  }

  void invalidateEvents(const Particle &);

  void addEvents(Particle &);

  void popNextEvent();

  void pushEvent(const Event &);

  void stream(const double dt) { sorter->stream(dt); }

  void runNextEvent();

  friend magnet::xml::XmlStream &operator<<(magnet::xml::XmlStream &,
                                            const Scheduler &);

  static shared_ptr<Scheduler> getClass(const magnet::xml::Node &,
                                        dynamo::Simulation *const);

  virtual void operator<<(const magnet::xml::Node &);

  void rescaleTimes(const double &scale) { sorter->rescaleTimes(scale); }

  const shared_ptr<FEL> &getSorter() const { return sorter; }

  void rebuildSystemEvents() const;

  void addInteractionEvent(const Particle &, const size_t &) const;

  void addLocalEvent(const Particle &, const size_t &) const;

  virtual double getNeighbourhoodDistance() const = 0;
  virtual std::unique_ptr<IDRange>
  getParticleNeighbours(const Particle &) const = 0;
  virtual std::unique_ptr<IDRange>
  getParticleNeighbours(const Vector &) const = 0;
  virtual std::unique_ptr<IDRange>
  getParticleLocals(const Particle &) const = 0;

protected:
  mutable shared_ptr<FEL> sorter;

  size_t _interactionRejectionCounter;
  size_t _localRejectionCounter;

  virtual void outputXML(magnet::xml::XmlStream &) const = 0;
};
} // namespace dynamo
