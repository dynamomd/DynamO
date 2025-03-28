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
#include <boost/circular_buffer.hpp>
#include <dynamo/outputplugins/tickerproperty/ticker.hpp>
#include <magnet/math/vector.hpp>
#include <vector>

namespace dynamo {
class OPMSDOrientationalCorrelator : public OPTicker {
public:
  OPMSDOrientationalCorrelator(const dynamo::Simulation *,
                               const magnet::xml::Node &);

  virtual void initialise();

  void output(magnet::xml::XmlStream &);

  virtual void operator<<(const magnet::xml::Node &);

  typedef std::pair<Vector, Vector> RUpair;

protected:
  virtual void stream(double) {}
  virtual void ticker();

  void accPass();

  std::vector<boost::circular_buffer<RUpair>> historicalData;
  std::vector<double> stepped_data_parallel, stepped_data_perpendicular,
      stepped_data_rotational_legendre1, stepped_data_rotational_legendre2;

  size_t length;
  size_t currCorrLength;
  size_t ticksTaken;
  bool notReady;
};
} // namespace dynamo
