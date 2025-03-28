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
#include <dynamo/outputplugins/tickerproperty/ticker.hpp>
#include <list>
#include <vector>

namespace dynamo {
class OPMSD;
class Topology;

class OPPeriodicMSD : public OPTicker {
public:
  OPPeriodicMSD(const dynamo::Simulation *, const magnet::xml::Node &);

  virtual void initialise();

  void output(magnet::xml::XmlStream &);

protected:
  virtual void stream(double) {}
  virtual void ticker();

  typedef std::pair<double, double> localpair;

  std::list<localpair> results;

  typedef std::pair<const Topology *, std::list<localpair>> localpair2;

  std::vector<localpair2> structResults;
  std::vector<std::vector<localpair>> speciesData;

  shared_ptr<const OPMSD> ptrOPMSD;
};
} // namespace dynamo
