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
class SysSnapshot : public System {
public:
  SysSnapshot(dynamo::Simulation *, double, std::string, std::string, bool);
  SysSnapshot(dynamo::Simulation *, size_t, std::string, std::string, bool);

  virtual NEventData runEvent();

  virtual void initialise(size_t);

  virtual void operator<<(const magnet::xml::Node &) {}

  void setdt(double);

  void increasedt(double);

  const double &getPeriod() const { return _period; }

  virtual void replicaExchange(System &os) {
    SysSnapshot &s = static_cast<SysSnapshot &>(os);
    std::swap(dt, s.dt);
    std::swap(_period, s._period);
    std::swap(_applyBC, s._applyBC);
    std::swap(_format, s._format);
    std::swap(_saveCounter, s._saveCounter);
  }

  void setTickerPeriod(const double &);

protected:
  void eventCallback(const NEventData &);
  virtual void outputXML(magnet::xml::XmlStream &) const {}

  double _period;
  bool _applyBC;
  std::string _format;
  size_t _saveCounter;
  size_t _eventPeriod;
  size_t _lastEventCount;
};
} // namespace dynamo
