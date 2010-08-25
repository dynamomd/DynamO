/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "2partproperty.hpp"
#include <boost/foreach.hpp>
#include "../../dynamics/include.hpp"

OP2PP::OP2PP(const DYNAMO::SimData* t1,const char *t2):
  OutputPlugin(t1,t2)
{}

void 
OP2PP::eventUpdate(const IntEvent &event, 
		    const C2ParticleData &SDat) 
{
  stream(event.getdt());
  A2ParticleChange(SDat);
}

void 
OP2PP::eventUpdate(const CGlobEvent &event, const CNParticleData& SDat) 
{
  stream(event.getdt());

  BOOST_FOREACH(const C2ParticleData& pData, SDat.L2partChanges)
    A2ParticleChange(pData);
}

void 
OP2PP::eventUpdate(const CLocalEvent &event, const CNParticleData& SDat) 
{
  stream(event.getdt());

  BOOST_FOREACH(const C2ParticleData& pData, SDat.L2partChanges)
    A2ParticleChange(pData);
}

void 
OP2PP::eventUpdate(const CSystem&, const CNParticleData& SDat, const Iflt& dt)
{
  stream(dt);

  BOOST_FOREACH(const C2ParticleData& pData, SDat.L2partChanges)
    A2ParticleChange(pData);
}
