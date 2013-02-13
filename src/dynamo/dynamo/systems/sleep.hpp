/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>
    Copyright (C) 2011  Sebastian Gonzalez <tsuresuregusa@gmail.com>

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

#include <dynamo/systems/system.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/ranges/IDRange.hpp>
#include <map>

namespace dynamo {
  class SSleep: public System
  {
  public:
    SSleep(const magnet::xml::Node& XML, dynamo::Simulation*);

    SSleep(dynamo::Simulation*, std::string, IDRange*, double);
  
    virtual void runEvent() const;

    virtual void initialise(size_t);

    virtual void operator<<(const magnet::xml::Node&);

  protected:
    virtual void outputXML(magnet::xml::XmlStream&) const;

    void particlesUpdated(const NEventData&);

    void recalculateTime();

    bool sleepCondition(const Particle& part, const Vector& g, const Vector& vel = Vector(0,0,0));

    shared_ptr<IDRange> _range;
    double _sleepDistance;
    double _sleepTime;
    double _sleepVelocity;

    mutable std::map<size_t, Vector> stateChange;

    std::vector<std::pair<Vector, long double> > _lastData;
  };
}
