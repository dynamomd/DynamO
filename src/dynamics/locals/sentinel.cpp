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
#include "../liouvillean/liouvillean.hpp"
#include "localEvent.hpp"
#include "../NparticleEventData.hpp"

CLSentinel::CLSentinel(const XMLNode&, const DYNAMO::SimData* tmp):
  CLocal(tmp, "GlobalSentinel")
{}

CLocalEvent 
CLSentinel::getEvent(const CParticle& part) const
{
  Sim->Dynamics.Liouvillean().updateParticle(part);
  
  return CLocalEvent(part, HUGE_VAL, NONE, *this);
}

CNParticleData 
CLSentinel::runEvent(const CLocalEvent&) const
{ return CNParticleData(); }

bool 
CLSentinel::isInCell(const CVector<>& Origin, const CVector<>& CellDim) const
{
  return false;
}

void 
CLSentinel::initialise(size_t nID)
{
  ID = nID;
}

void 
CLSentinel::operator<<(const XMLNode& XML)
{}

void 
CLSentinel::outputXML(xmlw::XmlStream& XML) const
{
    XML << xmlw::attr("Type") << "Sentinel";
}
