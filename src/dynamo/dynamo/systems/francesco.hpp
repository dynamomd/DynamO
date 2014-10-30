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
  class SysFrancesco: public System
  {
  public:
    SysFrancesco(const magnet::xml::Node& XML, dynamo::Simulation*);

    SysFrancesco(dynamo::Simulation*, double, double, std::string);
  
    virtual NEventData runEvent();

    virtual void initialise(size_t);

    virtual void operator<<(const magnet::xml::Node&);

    double getTemperature() const { return Temp; }
    double getReducedTemperature() const;
    void setTemperature(double nT) { Temp = nT; sqrtTemp = std::sqrt(Temp); }
    void setReducedTemperature(double nT);

  protected:
    virtual void outputXML(magnet::xml::XmlStream&) const;
    double meanFreeTime;
    double Temp, sqrtTemp;
    double _R;
    size_t dimensions;
    size_t eventCount;
    size_t lastlNColl;

    double getGhostt() const;
  
    shared_ptr<IDRange> range;
  };
}
