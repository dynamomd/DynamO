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

#include "MFT.hpp"
#include <boost/foreach.hpp>
#include "../../base/is_simdata.hpp"
#include "../../dynamics/dynamics.hpp"
#include "../../dynamics/species/species.hpp"
#include "../../dynamics/1particleEventData.hpp"
#include "../../dynamics/units/units.hpp"

COPMFT::COPMFT(const DYNAMO::SimData* tmp):
  COP1PP(tmp,"MeanFreeLength", 250),
  collisionHistoryLength(10)
{}

void
COPMFT::initialise()
{
  lastTime.resize(Sim->lN, boost::circular_buffer<Iflt>(collisionHistoryLength, 0.0));
  
  std::vector<C1DHistogram> vecTemp;
  
  vecTemp.resize(collisionHistoryLength, 
		 C1DHistogram(Sim->Dynamics.units().unitLength() * 0.05));
  
  data.resize(Sim->lN, vecTemp);
}

void 
COPMFT::A1ParticleChange(const C1ParticleData& PDat)
{
  //We ignore stuff that hasn't had an event yet

  for (size_t collN = 0; collN < collisionHistoryLength; ++collN)
    if (lastTime[PDat.getParticle().getID()][collN] != 0.0)
      {
	data[PDat.getSpecies().getID()][collN]
	  .addVal(Sim->dSysTime 
		  - lastTime[PDat.getParticle().getID()][collN]);
      }
  
  lastTime[PDat.getParticle().getID()].push_front(Sim->dSysTime);
}

void
COPMFT::output(xmlw::XmlStream &XML)
{
  XML << xmlw::tag("MFT");
  
  for (size_t id = 0; id < data.size(); ++id)
    {
      XML << xmlw::tag("Species")
	  << xmlw::attr("Name")
	  << Sim->Dynamics.getSpecies()[id].getName();
      
      for (size_t collN = 0; collN < collisionHistoryLength; ++collN)
	{
	  XML << xmlw::tag("Collisions")
	      << xmlw::attr("val") << collN + 1;

	  data[id].outputHistogram(XML, 1.0 / Sim->Dynamics.units().unitLength());

	  XML << xmlw::tag("Collisions");
	}
	
      XML << xmlw::tag("Species");
    }
  
  XML << xmlw::endtag("MFT");
}

