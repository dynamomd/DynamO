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
#include <dynamo/systems/system.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/ranges/IDRange.hpp>
#include <dynamo/interactions/potentials/potential.hpp>
#include <map>

namespace dynamo {
  class SysUmbrella: public System
  {
  public:
    SysUmbrella(const magnet::xml::Node& XML, dynamo::Simulation*);
  
    virtual void runEvent();

    virtual void initialise(size_t);

    virtual void operator<<(const magnet::xml::Node&);

    virtual void outputData(magnet::xml::XmlStream&) const;

  protected:
    virtual void outputXML(magnet::xml::XmlStream&) const;

    void particlesUpdated(const NEventData&);

    void recalculateTime();

    std::size_t _stepID;
    shared_ptr<Potential> _potential;
    shared_ptr<IDRange> range1;
    shared_ptr<IDRange> range2;    
    double _energyScale;
    double _lengthScale;
    std::map<size_t, double> _histogram;
    long double _lastSystemTime;
  };
}
