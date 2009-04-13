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

COPMomentum::COPMomentum(const DYNAMO::SimData* tmp, const XMLNode&):
  COP1PP(tmp,"Momentum", 250),
  accMom(0,0,0), accMomsq(0,0,0), sysMom(0,0,0)
{}

void
COPMomentum::initialise()
{  
  accMom = Vector (0,0,0);
  accMomsq = Vector (0,0,0);
  sysMom = Vector (0,0,0);

  BOOST_FOREACH(const smrtPlugPtr<CSpecies>& spec, Sim->Dynamics.getSpecies())
    BOOST_FOREACH(const size_t& ID, *spec->getRange())
    sysMom += spec->getMass() * Sim->vParticleList[ID].getVelocity();
}

void 
COPMomentum::A1ParticleChange(const C1ParticleData& PDat)
{
  sysMom += PDat.getDeltaP();
}

void 
COPMomentum::stream(const Iflt& dt)
{
  Vector  tmp(sysMom * dt);
  accMom += tmp;
  for (size_t i(0); i < NDIM; ++i)
    accMomsq[i] += sysMom[i] * tmp[i];
}

void
COPMomentum::output(xmlw::XmlStream &XML)
{
  XML << xmlw::tag("Momentum")
      << xmlw::tag("Current")
      << sysMom / Sim->Dynamics.units().unitMomentum()
      << xmlw::endtag("Current")
      << xmlw::tag("Avg") 
      << accMom / (Sim->dSysTime * Sim->Dynamics.units().unitMomentum())
      << xmlw::endtag("Avg")    
      << xmlw::tag("SqAvg") 
      << accMomsq / (Sim->dSysTime * Sim->Dynamics.units().unitMomentum() * Sim->Dynamics.units().unitMomentum())
      << xmlw::endtag("SqAvg")
      << xmlw::endtag("Momentum");
}
