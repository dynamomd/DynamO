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

#include <dynamo/outputplugins/2partproperty/2partproperty.hpp>
#include <dynamo/include.hpp>
#include <boost/foreach.hpp>

namespace dynamo {
  OP2PP::OP2PP(const dynamo::Simulation* t1,const char *t2):
    OutputPlugin(t1,t2)
  {}

  void 
  OP2PP::eventUpdate(const IntEvent &event, 
		     const PairEventData &SDat) 
  {
    stream(event.getdt());
    A2ParticleChange(SDat);
  }

  void 
  OP2PP::eventUpdate(const GlobalEvent &event, const NEventData& SDat) 
  {
    stream(event.getdt());

    BOOST_FOREACH(const PairEventData& pData, SDat.L2partChanges)
      A2ParticleChange(pData);
  }

  void 
  OP2PP::eventUpdate(const LocalEvent &event, const NEventData& SDat) 
  {
    stream(event.getdt());

    BOOST_FOREACH(const PairEventData& pData, SDat.L2partChanges)
      A2ParticleChange(pData);
  }

  void 
  OP2PP::eventUpdate(const System&, const NEventData& SDat, const double& dt)
  {
    stream(dt);

    BOOST_FOREACH(const PairEventData& pData, SDat.L2partChanges)
      A2ParticleChange(pData);
  }
}
