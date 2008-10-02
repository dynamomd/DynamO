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

#include "momentum.hpp"
#include <boost/foreach.hpp>
#include <cmath>
#include "../../dynamics/interactions/captures.hpp"
#include "../../dynamics/include.hpp"
#include "../../dynamics/interactions/intEvent.hpp"
#include "../../base/is_simdata.hpp"
#include "../../datatypes/vector.xml.hpp"
#include "../../datatypes/pluginpointer.hpp"

COPMomentum::COPMomentum(const DYNAMO::SimData* tmp):
  COP1PP(tmp,"Momentum")
{}

void
COPMomentum::initialise()
{
  CVector<> tmp(0.0);

  BOOST_FOREACH(const CSpecies& spec, Sim->Dynamics.getSpecies())
    {
      CVector<> tmp(0.0);
      BOOST_FOREACH(unsigned int ID, *spec.getRange())
	  {
	    const CParticle& part = Sim->vParticleList[ID];
	    tmp += part.getVelocity() * spec.getMass();
	  }
      MomentumVal[&spec] = tmp;
    }
}
  
void 
COPMomentum::A1ParticleChange(const C1ParticleData& PDat)
{
  MomentumVal[&(Sim->Dynamics.getSpecies(PDat.getParticle()))] 
    += PDat.getDeltaP();
}

void 
COPMomentum::stream(Iflt dt)
{
  typedef std::pair<const CSpecies* const, CVector<> > mypair;
  BOOST_FOREACH(mypair& dat, MomentumVal)
    MomentumAcc[dat.first] += dt * dat.second;
}

void
COPMomentum::output(xmlw::XmlStream &XML)
{

  XML << xmlw::tag("Momentum");

  typedef std::pair<const CSpecies* const, CVector<> > mypair;
  BOOST_FOREACH(mypair& dat, MomentumAcc)
    XML << xmlw::tag("Species")
	<< xmlw::attr("Name") << dat.first->getName()
	<< (dat.second / (Sim->dSysTime * Sim->Dynamics.units().unitMomentum()))
	<< xmlw::endtag("Species");
  
  XML << xmlw::endtag("Momentum");
}

void
COPMomentum::periodicOutput()
{
//  Iflt unitPwrloss = Sim->Dynamics.units().unitLength()
//    * pow(Sim->Dynamics.units().unitTime(),3) / Sim->Dynamics.units().unitMass();
//   
//  I_Pcout() << "T " << getTheta() << ", U " << intECurrent << ", <T> " 
//	    << getAvgTheta() << ", <PwrLoss> " 
//	    << sumPowerLoss * unitPwrloss
//    /(Sim->Dynamics.units().simVolume()*Sim->dSysTime)
//	    << ", ";
//
}
