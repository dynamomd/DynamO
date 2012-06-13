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
#include <dynamo/dynamics/eventtypes.hpp>
#include <dynamo/outputplugins/1partproperty/1partproperty.hpp>
#include <magnet/math/histogram.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <vector>
#include <map>

namespace dynamo {
  class OPCollEnergyChange: public OP1PP
  {
  public:
    OPCollEnergyChange(const dynamo::SimData*, const magnet::xml::Node&);

    void A1ParticleChange(const ParticleEventData&);

    void A2ParticleChange(const PairEventData&);

    void stream(const double&) {}

    void output(magnet::xml::XmlStream &); 

    void periodicOutput() {}

    virtual void initialise();

    void operator<<(const magnet::xml::Node&);

  protected:  
    double binWidth; 

    static double KEBinWidth;

    struct histogram: public magnet::math::Histogram<>
    {
      histogram(): magnet::math::Histogram<>(OPCollEnergyChange::KEBinWidth) {}
    };

    typedef boost::tuple<size_t, size_t, EEventType> mapkey;

    std::map<mapkey, histogram> collisionKE;

    std::vector<magnet::math::Histogram<> > data;
    magnet::math::Histogram<> specialhist;
  };
}
