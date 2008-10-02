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

#include "rdotv.hpp"
#include "../../dynamics/include.hpp"
#include <boost/foreach.hpp>
#include "../../base/is_simdata.hpp"
#include "../0partproperty/collMatrix.hpp"

COPRdotV::COPRdotV(const DYNAMO::SimData* tmp):
  COutputPlugin(tmp, "RdotV"),
  collMatrixPlug(NULL)
{}

void 
COPRdotV::initialise()
{
  collMatrixPlug = Sim->getOutputPlugin<COPCollMatrix>();
}

void 
COPRdotV::eventUpdate(const CIntEvent& iEvent, const C2ParticleData& pDat)
{
  rvdotacc[std::make_pair(iEvent.getType(),collMatrixPlug->getID(iEvent.getInteraction()))].addVal(pDat.rij % pDat.particle1_.getDeltaP());
}

void 
COPRdotV::eventUpdate(const CGlobEvent& globEvent, const CNParticleData& SDat)
{
  BOOST_FOREACH(const C2ParticleData& pDat, SDat.L2partChanges)
    {
      rvdotacc[std::make_pair(globEvent.getType(), collMatrixPlug->getID(globEvent.getGlobal()))].addVal(pDat.rij % pDat.particle1_.getDeltaP());
    }
}

void
COPRdotV::eventUpdate(const CSystem& sysEvent, const CNParticleData& SDat, const Iflt&)
{
  BOOST_FOREACH(const C2ParticleData& pDat, SDat.L2partChanges)
    {
      rvdotacc[std::make_pair(sysEvent.getType(), collMatrixPlug->getID(sysEvent))].addVal(pDat.rij % pDat.particle1_.getDeltaP());
    } 
}

void
COPRdotV::output(xmlw::XmlStream &XML)
{
  XML << xmlw::tag("rdotV");

  BOOST_FOREACH(mappair& pair1, rvdotacc)
    {
      XML << xmlw::tag("Element")
	  << xmlw::attr("Type") 
	  << CIntEvent::getCollEnumName(pair1.first.first)
	  << xmlw::attr("EventName") 
	  << collMatrixPlug->getName(pair1.first.second)
	  << xmlw::attr("Val") << pair1.second.getAvg() 
	/ (Sim->Dynamics.units().unitVelocity() 
	   * Sim->Dynamics.units().unitLength()
	   * Sim->Dynamics.units().unitMass())
	  << xmlw::endtag("Element");
    }
  
  
  XML << xmlw::endtag("rdotV");
}
