/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

OPMFT::OPMFT(const DYNAMO::SimData* tmp, const XMLNode& XML):
  OP1PP(tmp,"MeanFreeLength", 250),
  collisionHistoryLength(10),
  binwidth(0.01)
{
  operator<<(XML);
}

void 
OPMFT::operator<<(const XMLNode& XML)
{
  try 
    {
      if (XML.isAttributeSet("binwidth"))
	binwidth = boost::lexical_cast<Iflt>(XML.getAttribute("binwidth"));
      
      if (XML.isAttributeSet("length"))
	collisionHistoryLength 
	  = boost::lexical_cast<size_t>(XML.getAttribute("length"));
    }
  catch (boost::bad_lexical_cast&)
    {
      D_throw() << "Failed a lexical cast in OPMFL";
    }
}

void
OPMFT::initialise()
{
  lastTime.resize(Sim->N, 
		  boost::circular_buffer<Iflt>(collisionHistoryLength, 0.0));
  
  std::vector<C1DHistogram> vecTemp;
  
  vecTemp.resize(collisionHistoryLength, 
		 C1DHistogram(Sim->dynamics.units().unitTime() * binwidth));
  
  data.resize(Sim->dynamics.getSpecies().size(), vecTemp);
}

void 
OPMFT::A1ParticleChange(const ParticleEventData& PDat)
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
OPMFT::output(xml::XmlStream &XML)
{
  XML << xml::tag("MFT");
  
  for (size_t id = 0; id < data.size(); ++id)
    {
      XML << xml::tag("Species")
	  << xml::attr("Name")
	  << Sim->dynamics.getSpecies()[id]->getName();
      
      for (size_t collN = 0; collN < collisionHistoryLength; ++collN)
	{
	  XML << xml::tag("Collisions")
	      << xml::attr("val") << collN + 1;
	  
	  data[id][collN].outputHistogram
	    (XML, 1.0 / Sim->dynamics.units().unitTime());
	  
	  XML << xml::endtag("Collisions");
	}
	
      XML << xml::endtag("Species");
    }
  
  XML << xml::endtag("MFT");
}

