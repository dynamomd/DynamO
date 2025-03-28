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
#include <dynamo/outputplugins/outputplugin.hpp>
#include <fstream>

namespace dynamo {
class OPTrajectory : public OutputPlugin {
public:
  OPTrajectory(const dynamo::Simulation *, const magnet::xml::Node &);

  OPTrajectory(const OPTrajectory &);

  ~OPTrajectory() {}

  void eventUpdate(const Event &, const NEventData &);

  virtual void replicaExchange(OutputPlugin &) {
    M_throw()
        << "This output plugin hasn't been prepared for changes of system";
  }

  virtual void initialise();

  virtual void output(magnet::xml::XmlStream &);

private:
  mutable std::ofstream logfile;
};
} // namespace dynamo
