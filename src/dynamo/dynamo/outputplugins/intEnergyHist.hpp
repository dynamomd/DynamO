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
#include <magnet/math/histogram.hpp>
#include <unordered_map>

namespace dynamo {
class OPMisc;

class OPIntEnergyHist : public OutputPlugin {
public:
  OPIntEnergyHist(const dynamo::Simulation *, const magnet::xml::Node &);

  virtual void initialise();

  virtual void eventUpdate(const Event &, const NEventData &);

  virtual void output(magnet::xml::XmlStream &);

  virtual void replicaExchange(OutputPlugin &);

  void operator<<(const magnet::xml::Node &);

  std::unordered_map<int, double> getImprovedW() const;
  inline double getBinWidth() const { return intEnergyHist.getBinWidth(); }

protected:
  magnet::math::HistogramWeighted<> intEnergyHist;
  shared_ptr<const OPMisc> _ptrOPMisc;
  double binwidth;
};
} // namespace dynamo
