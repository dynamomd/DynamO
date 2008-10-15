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

#include "collenergychange.hpp"
#include <boost/foreach.hpp>
#include "../../base/is_simdata.hpp"
#include "../../dynamics/dynamics.hpp"
#include "../../dynamics/species/species.hpp"
#include "../../dynamics/1particleEventData.hpp"
#include "../../dynamics/units/units.hpp"

COPCollEnergyChange::COPCollEnergyChange(const DYNAMO::SimData* tmp):
  COP1PP(tmp,"MeanFreeLength", 250)
{}

void
COPCollEnergyChange::initialise()
{
  data.resize(Sim->Dynamics.getSpecies().size(), 
	      C1DHistogram(Sim->Dynamics.units().unitEnergy() * 0.05));
}

void 
COPCollEnergyChange::A1ParticleChange(const C1ParticleData& PDat)
{
  data[PDat.getSpecies().getID()]
    .addVal(PDat.getDeltae());
}

void
COPCollEnergyChange::output(xmlw::XmlStream &XML)
{
  XML << xmlw::tag("CollEnergyChange");
  
  for (size_t id = 0; id < data.size(); ++id)
    {
      XML << xmlw::tag("Species")
	  << xmlw::attr("Name")
	  << Sim->Dynamics.getSpecies()[id].getName();

      data[id].outputHistogram(XML, 1.0 / Sim->Dynamics.units().unitEnergy());

      XML << xmlw::endtag("Species");
    }

  XML << xmlw::endtag("CollEnergyChange");
}

