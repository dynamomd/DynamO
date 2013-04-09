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

namespace dynamo {
  //! \brief A System Event which periodically saves the state of the system.
  class SSnapshot: public System
  {
  public:
    SSnapshot(dynamo::Simulation*, double, std::string, std::string, bool);
  
    virtual void runEvent() const;

    virtual void initialise(size_t);

    virtual void operator<<(const magnet::xml::Node&) {}

    void setdt(double);

    void increasedt(double);

    void setTickerPeriod(const double&);

    const double& getPeriod() const { return _period; }
  protected:
    virtual void outputXML(magnet::xml::XmlStream&) const {}

    double _period;
    bool _applyBC;
    std::string _format;
    mutable size_t _saveCounter;
  };
}
