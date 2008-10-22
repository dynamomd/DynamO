/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "sentinel.hpp"
#include "globEvent.hpp"
#include "../NparticleEventData.hpp"
#include "../../extcode/xmlParser.h"
#include "../../extcode/xmlwriter.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../units/units.hpp"
#include "../ranges/1RAll.hpp"

CGSentinel::CGSentinel(const XMLNode &XML, const DYNAMO::SimData* ptrSim):
  CGlobal(ptrSim)
{
  globName = "CollisionSentinel";
  operator<<(XML);

  I_cout() << "Sentinel Loaded";
}

CGSentinel::CGSentinel(const DYNAMO::SimData* nSim):
  CGlobal(new CRAll(nSim), nSim)
{
  globName = "CollisionSentinel";
  I_cout() << "Sentinel Loaded";
}

CGlobEvent 
CGSentinel::getEvent(const CParticle& part) const
{
  Sim->Dynamics.Liouvillean().updateParticle(part);

  return CGlobEvent(part, Sim->Dynamics.Liouvillean().getHalfBoxTraversalTime(part), VIRTUAL, *this);
}

CNParticleData
CGSentinel::runEvent(const CGlobEvent& event) const
{
  return CNParticleData();
}

void 
CGSentinel::initialise(size_t nID)
{
  ID=nID;
}

void 
CGSentinel::operator<<(const XMLNode& XML)
{}

void
CGSentinel::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "Sentinel";
}
