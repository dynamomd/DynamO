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

#include "coll_vel_dist.hpp"
#include <cmath>
#include <boost/foreach.hpp>
#include "../simulation/particle.hpp"
#include "../dynamics/dynamics.hpp"
#include "../dynamics/interactions/intEvent.hpp"
#include "../dynamics/interactions/intEventData.hpp"
#include "../dynamics/units/units.hpp"
#include "ke.hpp"
#include "../base/is_simdata.hpp"

//When to collect the velocity distribution
#define collectFreq 100

Iflt CollVelHist::initVal = 1.0;

COPCollVelDist::COPCollVelDist(DYNAMO::SimData* tmp):
  COutputPlugin(tmp,"CollVelDistribution"),
  ptrKE(NULL)
{
  CollVelHist::initVal = sqrt(getkT())*0.02;
  
  BOOST_FOREACH(const smrtPlugPtr<COutputPlugin>& plugin, Sim->outputPlugins)
    if (dynamic_cast<const COPKE*>(plugin.get_ptr()) != NULL)
      {
	ptrKE = dynamic_cast<const COPKE*>(plugin.get_ptr());
	return;
      }
  
  I_throw() << "The Velocity distribution plugin(s) require COPKE <kinetic energy> plugin to be loaded, so it can normalise correctly";
}

void 
COPCollVelDist::collisionUpdate(const CIntEvent &collision, const CIntEventData &preColl)
{
  vf[collision.getType()].addVal(preColl.oldVelVec2.length());
  
  if (collision.getParticle1().getID() != collision.getParticle2().getID())
    vf[collision.getType()].addVal(preColl.oldVelVec2.length());
}

void
COPCollVelDist::output(xmlw::XmlStream &XML)
{
  XML << xmlw::tag("Coll_Vel_Dist");
  
  typedef std::pair<const EEventType, CollVelHist> mapPair;
  BOOST_FOREACH(mapPair& tmp, vf)
    {
      XML << xmlw::tag(CIntEvent::getCollEnumName(tmp.first));
      tmp.second.outputHistogram(XML, 1.0/sqrt(ptrKE->getAvgkT()));
      XML << xmlw::endtag(CIntEvent::getCollEnumName(tmp.first));
    }

  XML << xmlw::endtag("Coll_Vel_Dist");
}

