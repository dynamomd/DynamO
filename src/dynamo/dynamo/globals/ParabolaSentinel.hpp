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
#include <dynamo/globals/global.hpp>
#include <vector>

namespace dynamo {
class GParabolaSentinel : public Global {
public:
  GParabolaSentinel(dynamo::Simulation *, const std::string &);

  virtual ~GParabolaSentinel() {}

  virtual Event getEvent(const Particle &) const;

  virtual void runEvent(Particle &, const double);

  virtual void operator<<(const magnet::xml::Node &) {}

protected:
  virtual void outputXML(magnet::xml::XmlStream &) const {}
};
} // namespace dynamo
