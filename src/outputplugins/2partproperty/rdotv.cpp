/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

COPRdotV::COPRdotV(const DYNAMO::SimData* tmp, const XMLNode&):
  COutputPlugin(tmp, "RdotV")
{}

void 
COPRdotV::initialise()
{}

void 
COPRdotV::eventUpdate(const CIntEvent& iEvent, const C2ParticleData& pDat)
{
  mapdata& ref = rvdotacc[mapKey(iEvent.getType(), getClassKey(iEvent))];

  ref.addVal(pDat.rij | pDat.particle1_.getDeltaP());
  ref.costheta.addVal(pDat.rij | pDat.vijold 
		      / (pDat.rij.nrm() * pDat.vijold.nrm()));
}

void 
COPRdotV::eventUpdate(const CGlobEvent& globEvent, const CNParticleData& SDat)
{
  BOOST_FOREACH(const C2ParticleData& pDat, SDat.L2partChanges)
    {
      mapdata& ref = rvdotacc[mapKey(globEvent.getType(), 
				     getClassKey(globEvent))];
      
      ref.addVal(pDat.rij | pDat.particle1_.getDeltaP());

      ref.costheta.addVal(pDat.rij | pDat.vijold / (pDat.rij.nrm() * pDat.vijold.nrm()));
    }
}

void 
COPRdotV::eventUpdate(const CLocalEvent& localEvent, const CNParticleData& SDat)
{
  BOOST_FOREACH(const C2ParticleData& pDat, SDat.L2partChanges)
    {
      mapdata& ref = rvdotacc[mapKey(localEvent.getType(), 
				     getClassKey(localEvent))];
      
      ref.addVal(pDat.rij | pDat.particle1_.getDeltaP());

      ref.costheta.addVal(pDat.rij | pDat.vijold / (pDat.rij.nrm() * pDat.vijold.nrm()));
    }
}

void
COPRdotV::eventUpdate(const CSystem& sysEvent, const CNParticleData& SDat, const Iflt&)
{
  BOOST_FOREACH(const C2ParticleData& pDat, SDat.L2partChanges)
    {
      mapdata& ref = rvdotacc[mapKey(sysEvent.getType(), 
				     getClassKey(sysEvent))];
      
      ref.addVal(pDat.rij | pDat.particle1_.getDeltaP());

      ref.costheta.addVal(pDat.rij | pDat.vijold / (pDat.rij.nrm() * pDat.vijold.nrm()));
    } 
}

void
COPRdotV::output(xmlw::XmlStream &XML)
{
  XML << xmlw::tag("RdotV");
  
  typedef std::pair<const mapKey, mapdata> mappair;

  BOOST_FOREACH(const mappair& pair1, rvdotacc)
    {
      XML << xmlw::tag("Element")
	  << xmlw::attr("Type") 
	  << CIntEvent::getCollEnumName(pair1.first.first)
	  << xmlw::attr("EventName") 
	  << getName(pair1.first.second, Sim)
	  << xmlw::attr("RijdotDeltaMomentum") << pair1.second.getAvg()
	/ (Sim->Dynamics.units().unitVelocity() 
	   * Sim->Dynamics.units().unitLength()
	   * Sim->Dynamics.units().unitMass());
      
      pair1.second.costheta.outputHistogram(XML, 1.0);
      
      XML << xmlw::endtag("Element");
    }

    XML << xmlw::endtag("RdotV");
}
