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
#include <dynamo/eventtypes.hpp>
#include <dynamo/outputplugins/eventtypetracking.hpp>
#include <magnet/math/vector.hpp>
#include <map>
#include <magnet/math/histogram.hpp>

namespace dynamo {
  class Particle;

  using namespace EventTypeTracking;

  class OPBrenner: public OutputPlugin
  {
  public:
    OPBrenner(const dynamo::Simulation*, const magnet::xml::Node&);
    ~OPBrenner();

    virtual void initialise();

    virtual void eventUpdate(const Event&, const NEventData&);

    void output(magnet::xml::XmlStream &);

    virtual void replicaExchange(OutputPlugin& plug) { 
      M_throw() << "Not implemented";
    }
  
  protected:
    magnet::math::HistogramWeighted<> _sysmomentum_hist[3];
    Vector _sysMomentum;
  };
}
