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

namespace dynamo {
  class SysAndersen: public System
  {
  public:
    SysAndersen(const magnet::xml::Node& XML, dynamo::Simulation*);

    SysAndersen(dynamo::Simulation*, double, double, std::string);
  
    virtual void runEvent();

    virtual void initialise(size_t);

    virtual void operator<<(const magnet::xml::Node&);

    double getTemperature() const { return Temp; }
    double getReducedTemperature() const;
    void setTemperature(double nT) { Temp = nT; sqrtTemp = std::sqrt(Temp); }
    void setReducedTemperature(double nT);

    virtual void replicaExchange(System& os) { 
      SysAndersen& s = static_cast<SysAndersen&>(os);
      std::swap(dt, s.dt);
      std::swap(meanFreeTime, s.meanFreeTime);
      std::swap(Temp, s.Temp);
      std::swap(sqrtTemp, s.sqrtTemp);
      std::swap(tune, s.tune);
      std::swap(dimensions, s.dimensions);
      std::swap(setPoint, s.setPoint);
      std::swap(eventCount, s.eventCount);
      std::swap(lastlNColl, s.lastlNColl);
      std::swap(setFrequency, s.setFrequency);
    }
  
  protected:
    virtual void outputXML(magnet::xml::XmlStream&) const;
    double meanFreeTime;
    double Temp, sqrtTemp;
    bool tune;
    size_t dimensions;
    double setPoint;
    size_t eventCount;
    size_t lastlNColl;
    size_t setFrequency;

    double getGhostt() const;
  
    shared_ptr<IDRange> range;
  };
}
