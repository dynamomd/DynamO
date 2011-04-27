/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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

#include "momentum.hpp"
#include <boost/foreach.hpp>
#include <cmath>
#include "../../dynamics/interactions/captures.hpp"
#include "../../dynamics/include.hpp"
#include "../../dynamics/interactions/intEvent.hpp"
#include "../../base/is_simdata.hpp"
#include "../../datatypes/vector.xml.hpp"

OPMomentum::OPMomentum(const DYNAMO::SimData* tmp, const magnet::xml::Node&):
  OP1PP(tmp,"Momentum", 250),
  accMom(0,0,0), accMomsq(0,0,0), sysMom(0,0,0)
{}

void
OPMomentum::initialise()
{  
  accMom = Vector (0,0,0);
  accMomsq = Vector (0,0,0);
  sysMom = Vector (0,0,0);

  BOOST_FOREACH(const magnet::ClonePtr<Species>& spec, Sim->dynamics.getSpecies())
    BOOST_FOREACH(const size_t& ID, *spec->getRange())
    sysMom += spec->getMass(ID) * Sim->particleList[ID].getVelocity();
}

void 
OPMomentum::A1ParticleChange(const ParticleEventData& PDat)
{
  sysMom += PDat.getDeltaP();
}

void 
OPMomentum::stream(const double& dt)
{
  Vector  tmp(sysMom * dt);
  accMom += tmp;
  for (size_t i(0); i < NDIM; ++i)
    accMomsq[i] += sysMom[i] * tmp[i];
}

void
OPMomentum::output(xml::XmlStream &XML)
{
  XML << xml::tag("Momentum")
      << xml::tag("Current")
      << sysMom / Sim->dynamics.units().unitMomentum()
      << xml::endtag("Current")
      << xml::tag("Avg") 
      << accMom / (Sim->dSysTime * Sim->dynamics.units().unitMomentum())
      << xml::endtag("Avg")    
      << xml::tag("SqAvg") 
      << accMomsq / (Sim->dSysTime * Sim->dynamics.units().unitMomentum() * Sim->dynamics.units().unitMomentum())
      << xml::endtag("SqAvg")
      << xml::endtag("Momentum");
}
