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
#include <dynamo/outputplugins/1partproperty/1partproperty.hpp>
#include <magnet/math/histogram.hpp>
#include <boost/circular_buffer.hpp>
#include <vector>

namespace dynamo {
  class OPMFT: public OP1PP
  {
  public:
    OPMFT(const dynamo::Simulation*, const magnet::xml::Node&);

    void A1ParticleChange(const ParticleEventData&);

    void stream(const double&) {}

    void output(magnet::xml::XmlStream &); 

    void periodicOutput() {}

    virtual void initialise();

    virtual void operator<<(const magnet::xml::Node&);
  
  protected:
    size_t collisionHistoryLength;
  
    double binwidth;

    //! Each particles last collision times
    std::vector<boost::circular_buffer<double> > lastTime;

    //! A histogram for each species
    std::vector<std::vector<magnet::math::Histogram<> > > data;
  };
}
