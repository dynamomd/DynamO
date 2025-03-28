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
class NEventData;

class System : public dynamo::SimBase {
public:
  System(dynamo::Simulation *);

  virtual ~System() {}

  inline void stream(const double &ndt) { dt -= ndt; }

  virtual NEventData runEvent() = 0;

  virtual void initialise(size_t) = 0;

  virtual void operator<<(const magnet::xml::Node &) = 0;

  Event getEvent() const { return Event(ID, dt, SYSTEM, type, ID); }

  friend magnet::xml::XmlStream &operator<<(magnet::xml::XmlStream &,
                                            const System &);

  static shared_ptr<System> getClass(const magnet::xml::Node &,
                                     dynamo::Simulation *);

  void setName(const std::string &tmp) { sysName = tmp; }

  const std::string &getName() const { return sysName; }

  inline double getdt() const { return dt; }

  inline const size_t &getID() const { return ID; }

  virtual void replicaExchange(System &s) {
    M_throw() << "The System \"" << getName() << "\"Not replica exchange safe";
  }

  virtual void outputData(magnet::xml::XmlStream &) const {}

protected:
  virtual void outputXML(magnet::xml::XmlStream &) const = 0;

  std::string sysName;
  double dt;
  EEventType type;
  size_t ID;
};
} // namespace dynamo
