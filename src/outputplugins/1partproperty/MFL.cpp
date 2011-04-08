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

#include "MFL.hpp"
#include "../../base/is_simdata.hpp"
#include "../../dynamics/dynamics.hpp"
#include "../../dynamics/species/species.hpp"
#include "../../dynamics/1particleEventData.hpp"
#include "../../dynamics/units/units.hpp"
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>

OPMFL::OPMFL(const DYNAMO::SimData* tmp, const magnet::xml::Node& XML):
  OP1PP(tmp,"MeanFreeLength", 250),
  binwidth(0.01)
{
  operator<<(XML);
}

void 
OPMFL::operator<<(const magnet::xml::Node& XML)
{
  try 
    {
      binwidth = XML.getAttribute("binwidth").as<double>(0.01);
    }
  catch (boost::bad_lexical_cast&)
    {
      M_throw() << "Failed a lexical cast in OPMFL";
    }
}

void
OPMFL::initialise()
{
  lastTime.resize(Sim->N, 0.0);
  data.resize(Sim->dynamics.getSpecies().size(), 
	      C1DHistogram(Sim->dynamics.units().unitLength() * binwidth));
}

void 
OPMFL::A1ParticleChange(const ParticleEventData& PDat)
{
  //We ignore stuff that hasn't had an event yet
  if (lastTime[PDat.getParticle().getID()] != 0.0)
    { 
     data[PDat.getSpecies().getID()]
	.addVal(PDat.getParticle().getVelocity().nrm() 
	     * (Sim->dSysTime 
		- lastTime[PDat.getParticle().getID()]));
    }
  
  lastTime[PDat.getParticle().getID()] = Sim->dSysTime;
}

void
OPMFL::output(xml::XmlStream &XML)
{
  XML << xml::tag("MFL");
  
  for (size_t id = 0; id < data.size(); ++id)
    {
      XML << xml::tag("Species")
	  << xml::attr("Name")
	  << Sim->dynamics.getSpecies()[id]->getName();

      data[id].outputHistogram(XML, 1.0 / Sim->dynamics.units().unitLength());
      
      XML << xml::endtag("Species");
    }

  XML << xml::endtag("MFL");
}

