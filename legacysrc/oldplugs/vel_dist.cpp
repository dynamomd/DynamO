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

#include "vel_dist.hpp"
#include <cmath>
#include <boost/foreach.hpp>
#include "../simulation/particle.hpp"
#include "../dynamics/dynamics.hpp"
#include "../dynamics/units/units.hpp"
#include "ke.hpp"
#include "../base/is_simdata.hpp"

//When to collect the velocity distribution
#define collectFreq 100

COPVelDist::COPVelDist(DYNAMO::SimData* tmp):
  COutputPlugin(tmp,"VelDistribution"),
  ptrKE(NULL)
{
  for (int i=0; i<NDIM; i++)
    vf[i] = C1DHistogram(sqrt(getkT(i))*0.02);
  
  BOOST_FOREACH(const smrtPlugPtr<COutputPlugin>& plugin, Sim->outputPlugins)
    if (dynamic_cast<const COPKE*>(plugin.get_ptr()) != NULL)
      {
	ptrKE = dynamic_cast<const COPKE*>(plugin.get_ptr());
	return;
      }
  
  I_throw() << "The Velocity distribution plugin requires COPKE <kinetic energy> plugin to be loaded to normalise correctly";
}

COPVelDist::COPVelDist(const COPVelDist &C2):
  COutputPlugin(C2)
{
  for (int i=0; i<NDIM; i++)
    vf[i] = C2.vf[i];
}

void 
COPVelDist::collisionUpdate(const CIntEvent &collision, const CIntEventData &preColl)
{
  if (!(Sim->lNColl % collectFreq))
    BOOST_FOREACH( const CParticle & Part, Sim->vParticleList)
      for (int i = 0; i < NDIM; i++)
	vf[i].addVal(Part.getVelocity()[i]);
}

void
COPVelDist::output(xmlw::XmlStream &XML)
{
  XML << xmlw::tag("Vel_Dist");

  for (int i = 0; i < NDIM; i++)
    {
      char name[6] = "x-dim";
      name[0] += i;
      XML << xmlw::tag(name);
      vf[i].outputHistogram(XML, 1.0/(Sim->Dynamics.units().unitLength()*sqrt(ptrKE->getAvgTheta(i))));     
      XML << xmlw::endtag(name);
    }

  XML << xmlw::endtag("Vel_Dist");
}

