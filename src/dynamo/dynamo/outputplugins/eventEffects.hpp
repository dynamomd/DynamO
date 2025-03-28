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
#include <dynamo/eventtypes.hpp>
#include <dynamo/outputplugins/eventtypetracking.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/math/vector.hpp>
#include <map>

namespace dynamo {
class Particle;

using namespace EventTypeTracking;

class OPEventEffects : public OutputPlugin {
public:
  OPEventEffects(const dynamo::Simulation *, const magnet::xml::Node &);
  ~OPEventEffects();

  virtual void initialise();

  virtual void eventUpdate(const Event &, const NEventData &);

  void output(magnet::xml::XmlStream &);

  // This is fine to replica exchange as the interaction, global and
  // system lookups are done using id's
  virtual void replicaExchange(OutputPlugin &plug) {
    std::swap(counters, static_cast<OPEventEffects &>(plug).counters);
  }

protected:
  void newEvent(const EEventType &, const EventSourceKey &, const double &,
                const Vector &);

  struct counterData {
    counterData() : count(0), energyLoss(0), momentumChange({0, 0, 0}) {}
    unsigned long count;
    double energyLoss;
    Vector momentumChange;
  };

  std::map<EventKey, counterData> counters;
};
} // namespace dynamo
